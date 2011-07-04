#!/bin/sh -e

DATAFILE=befor-hyp.txt

NODEP=n
if [ "$1" = "-n" ]; then
	NODEP=y
	shift
fi

if [ -z "$*" ]; then
    echo "usage: `basename $0` [-n] targets" >&2
    echo "options:" >&2
    echo "    -n    do not build target dependencies, do not use .targets folder" >&2
    exit 1
fi

log() {
    printf "%s $$ [%s] %s\n" "`basename $0`" "`date '+%Y-%m-%d %H:%M:%S'`" "$*"
}

require() {
	if [ "$1" = "-f" ]; then
		$2
	else
		if [ "$NODEP" = "y" ]; then
			return
		fi
        mkdir -p .targets
		if [ ! -r .targets/$1 ]; then
			$1
		fi
	fi
}

completed() {
	if [ "$NODEP" = "n" ]; then
		mkdir -p .targets
		touch .targets/$1
	fi
}

simulate() {
    log "init"
    init.debug >init.log 2>init.errlog
    log "simulate"
    rf-plasma.release >simulate.log 2>simulate.errlog <1.stdin || true
    completed simulate
}

gen_times() {
    log "gen_times"
    ls *.out | egrep -o "N[0-9]+" | sed "s/N0*//" | N-to-T.sh $DATAFILE > times.out
}

outdata() {
    log "outdata"
    if [ ! -e out.cfg ]; then
        echo "file does not exist: out.cfg" >&2
        exit 1
    fi
    rm -f *.out
    scatter-by-N.sh `head -1 out.cfg` < $DATAFILE 
    gen_times
    completed outdata
}

lastdata() {
    log "lastdata"
    rm -f *.out
    scatter-by-N.sh `parse-results.sh ' ' ' ' ' ' 'print N;' < $DATAFILE` < $DATAFILE
    gen_times
    completed lastdata
}

prelastdata() {
    log "prelastdata"
    rm -f *.out
    scatter-by-N.sh `parse-results.sh ' ' 'print N' < $DATAFILE | tail -3 | head -1` < $DATAFILE
    gen_times
    completed prelastdata
}

filter() {
    log "filter"
    for f in `ls *.out | egrep "N[0-9]+"`; do
        local MK=`tail -1 $f | cut -f 1`
        awk -v mk=$MK 'BEGIN{l=-1;}int(NR*1000/mk)>l{l=int(NR*1000/mk);print;}' < $f > $f.tmp && mv $f.tmp $f
    done
}

plot() {
    log "plot"
    local LEN=`head -3 parskhgs.cfg | tail -1 | perl -ne 's/d/e/i; /\s*([0-9.\-+ed]*)/i; print $1;'`
    local MK=`head -1 parskhgs.cfg | perl -ne 's/d/e/i; /\s*([0-9.\-+ed]*)/i; print $1;'`
    local DT=`head -24 bconst.cfg | tail -1 | perl -ne 's/d/e/i; /\s*([0-9.\-+ed]*)/i; print $1;'`
    local EXTRA=' '
    local CS=3e8
    local NSTEP=`head -3 1.stdin | tail -1`
    ./plot.sh "$LEN" "$MK" "$EXTRA" "$CS" "$DT" "$NSTEP"
    completed plot
}

all() {
    require simulate
    require outdata
    require plot
}

last() {
    require simulate
    require lastdata
    require plot
}

clean() {
    rm -f *.out StateDump befor-hyp.txt Data1D.in simulate.log simulate.errlog init.log init.errlog
}

targets() {
    if ! [ -d .targets ]; then
        echo "There is no targets built in current folder"
        exit
    fi
    ls -lt .targets
}

for t in $*; do
    require -f $t
done

log "Complete"
