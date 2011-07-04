#/bin/sh -e

N="$1";
shift || true

parse-results.sh '
if(N=='$N') {print m "\t" ne "\t" vz "\t" vx "\t" Et "\t" Ex "\t" Hy "\t" Ez;}
' | amline.pl "$@"
