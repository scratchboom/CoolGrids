BUILDDIR=../build
BINDIR=../bin
CONFIG=debug
OPTS=-ffixed-line-length-132 -g -O0 -D"LINUX"
#CONFIG=release
#OPTS=-ffixed-line-length-132 -O3 -D"LINUX"
FC=gfortran
PP=fpp

ODIR=$(BUILDDIR)/$(CONFIG)
RFPLASMA=$(BINDIR)/rf-plasma.$(CONFIG)
INIT=$(BINDIR)/init.$(CONFIG)

.PHONY: all dirs clean

all: dirs $(RFPLASMA) $(INIT)

$(RFPLASMA): $(ODIR)/rf-plasma
	ln -f $< $@

$(INIT): $(ODIR)/init
	ln -f $< $@

$(ODIR)/rf-plasma: $(ODIR)/Hyperbolic_1D6.i.for $(ODIR)/rhs-calc.i.for
	$(FC) $(OPTS) -o $@ $^ 

$(ODIR)/init: $(ODIR)/Set_IC_1D6.i.for
	$(FC) $(OPTS) -o $@ $^ 

$(ODIR)/%.i.for: %.for
	$(PP) $< | grep -v '^#' > $@

$(ODIR)/Hyperbolic_1D6.i.for: linux.conf

$(ODIR)/rhs-calc.i.for: linux.conf

dirs:
	mkdir -p $(BINDIR)
	mkdir -p $(ODIR)

clean:
	rm -rf $(ODIR)
