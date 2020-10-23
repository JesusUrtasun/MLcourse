# makefile for ggHiggs
# LO, NLO, NNLO and N3LO inclusive gg->H cross-section

# Ball, Bonvini, Forte, Marzani, Ridolfi,  arxiv:1303.3590
# Ball, Bonvini, Forte, Marzani, Ridolfi,  arxiv:1404.3204
# Bonvini, Marzani, Muselli, Rottoli,  arXiv:1603.08000 (version 3.0+)
# Ahmed, Bonvini, Kumar, Mathews, Rana, Ravindran, Rottoli,  arXiv:1606.00837 (for pseudo-scalar Higgs production, version 3.3+)
# Bonvini, Marzani, arXiv:1802.07758 (version 4.0+)
# Bonvini, arXiv:1805:08xxx (new exact N3LO, new small-x, small-x resummation, version 4.0+)

# author: Marco Bonvini, 2013

# path to CUBA 4.x library, if not already in LD_LIBRARY_PATH (example: /usr/local/Cuba-4.x)
CUBA_PATH = 

# path to gfortran, if not already in LD_LIBRARY_PATH (example: /usr/local/)
GFORTRAN_PATH = 

# path to LHAPDF
LHAPDF_PATH = $(shell lhapdf-config --prefix)

# path to GSL
GSL_PATH = $(shell gsl-config --prefix)

# path to HELL (if you want to use small-x resummation, download and make HELL-x from https://www.ge.infn.it/~bonvini/hell/ and specify the path here)
HELL_PATH = /home/jesus/Desktop/N3PDF/Tools/Bonvini_gg_Higgs/ggHiggs.v4.0/ggHiggs_Nspace/HELLx

######## do not edit below this line #######

HELLx = hell-x
HELLdata=$(HELL_PATH)/data/
ifneq ($(wildcard $(HELL_PATH)),)
withHELL= -DwithHELL -DHELLdata=\"$(HELLdata)\"
HELLlib = $(HELL_PATH)/lib$(HELLx).a
HELLhh  = $(HELL_PATH)/include/hell-x.hh
$(info HELL-x found. Compiling with small-x resummation feature available.)
else
withHELL=
HELLlib =
HELLhh  =
$(info HELL-x not found. Small-x resummation not available.)
endif
#$(info $$withHELL is [${withHELL}])

CC = g++
FF = gfortran

CFLAGS = -Wall -O2
CLIBS = -lgsl -lgslcblas -lgfortran -lcuba -lLHAPDF

INCLUDE_PATH = 
LIBS_PATH = 
#
ifneq "$(CUBA_PATH)" "" 
INCLUDE_PATH += -I$(CUBA_PATH)
LIBS_PATH += -L$(CUBA_PATH)
endif
#
ifneq "$(GFORTRAN_PATH)" "" 
LIBS_PATH += -L$(GFORTRAN_PATH)/lib
endif
#
ifneq "$(LHAPDF_PATH)" "" 
INCLUDE_PATH += -I$(LHAPDF_PATH)/include
LIBS_PATH += -L$(LHAPDF_PATH)/lib
endif
#
ifneq "$(GSL_PATH)" ""
INCLUDE_PATH += -I$(GSL_PATH)/include
LIBS_PATH += -L$(GSL_PATH)/lib
endif
#
ifneq "$(withHELL)" ""
INCLUDE_PATH += -I$(HELL_PATH)/include
#LIBS_PATH += -L$(HELL_PATH)
endif
#

CFLAGS += $(INCLUDE_PATH)

deps    = src/parameters.hh src/NLO_finite-mt/bdv.hh src/CFs_pointlike/CFs_pointlike.hh src/math/complex.hh src/small-x/finite-mt_small-x_coeffs.hh src/ggHiggs.hh src/math/mellin_Dk.hh src/deltaEW.cc

oNLOmt  = objects/HPL_full.o objects/bdv.o \
	objects/harlander.o objects/coefficients.o objects/HPL2.o objects/HPL3.o objects/HPL4.o

objects = $(oNLOmt) objects/nnlo_finite-mt.o objects/nlo_finite-mt.o \
	objects/NNLOpl.o objects/NNLOpl_F.o objects/N3LOpl.o objects/HSum.o objects/Sigma2.o \
	objects/complex.o objects/mellin_Dk.o  objects/integration.o objects/special_functions.o\
	objects/finite-mt_small-x_coeffs.o objects/hHiggs_new_var.o \
	objects/chebyshev.o objects/PDFs.o objects/msbar_mass.o

all: ggHiggs

ggH%: ggH%.cc libggHiggs.a $(HELLlib)
	$(CC) $(withHELL) -o $@ $< $(CFLAGS) libggHiggs.a $(HELLlib) $(LIBS_PATH) $(CLIBS) -lgfortran

libggHiggs.a: objects/ggHiggs.o $(objects)
	ar -r $@ $^

ifneq "$(withHELL)" ""
$(HELL_PATH)/lib$(HELLx).a: $(HELL_PATH)/include/hell-x.hh $(HELL_PATH)/src/hell-x.cc
	make lib$(HELLx).a -C $(HELL_PATH)
endif

objects/ggHiggs.o: src/ggHiggs.cc $(deps) $(HELLhh)
	$(CC) $(withHELL) $(CFLAGS) -o $@ -c $<
objects/%.o: src/math/%.cc src/math/%.hh
	$(CC) $(CFLAGS) -o $@ -c $<
objects/%.o: src/small-x/%.cc src/small-x/%.hh
	$(CC) $(CFLAGS) -o $@ -c $<
objects/%.o: src/%.f
	$(FF) -o $@ -c $<
objects/%_F.o: src/%.f
	$(FF) -o $@ -c $<
objects/%.o: src/%.cc src/%.hh
	$(CC) $(CFLAGS) -o $@ -c $<

# NLO finite m_top
objects/%.o: src/NLO_finite-mt/%.cc src/NLO_finite-mt/%.hh
	$(CC) $(CFLAGS) -o $@ -c $<
objects/%.o: src/NLO_finite-mt/%.F
	$(FF) -o $@ -c $<
objects/%.o: src/NLO_finite-mt/harlander/%.F
	$(FF) -o $@ -c $<
objects/%.o: src/NLO_finite-mt/harlander/polylogs/%.f
	$(FF) -o $@ -c $<

# NNLO finite m_top
objects/nnlo%.o: src/NNLO_finite-mt/nnlo%.cc src/NNLO_finite-mt/nnlo%.hh src/NNLO_finite-mt/nnlo_pars.hh src/NNLO_finite-mt/nnlo_mt_Nspace.cc src/NNLO_finite-mt/nnlo_mt_zspace.cc
	$(CC) $(CFLAGS) -c -o $@ $<
objects/nlo%.o: src/NNLO_finite-mt/nlo%.cc src/NNLO_finite-mt/nlo%.hh src/NNLO_finite-mt/nlo_pars.hh src/NNLO_finite-mt/nlo_mt_Nspace.cc src/NNLO_finite-mt/nlo_mt_zspace.cc
	$(CC) $(CFLAGS) -c -o $@ $<

# NNLO pointlike
objects/%_F.o: src/CFs_pointlike/%.f
	$(FF) -o $@ -c $<
objects/%.o: src/CFs_pointlike/%.cc src/CFs_pointlike/CFs_pointlike.hh
	$(CC) $(CFLAGS) -c -o $@ $<
objects/%.o: src/CFs_pointlike/NNLO_Exact_Mellin/%.cc src/CFs_pointlike/NNLO_Exact_Mellin/%.hh
	$(CC) $(CFLAGS) -c -o $@ $<

# NNLO anom dim
objects/%_F.o: src/anomalous_dimensions/%.f
	$(FF) -o $@ -c $<
objects/%.o: src/anomalous_dimensions/%.cc src/anomalous_dimensions/%.hh
	$(CC) $(CFLAGS) -o $@ -c $<

check:

clean: clean-tilde
	rm -f ggHiggs
	rm -f ggHiggs-verbose
	rm -f objects/*.o
	rm -f libggHiggs.a
	rm -f *.txt

clean-tilde:
	rm -f *~
	rm -f src/*~
	rm -f src/NLO_finite-mt/*~
	rm -f src/NNLO_finite-mt/*~
	rm -f src/math/*~

.PHONY: clean clean-tilde all check