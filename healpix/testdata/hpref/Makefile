# healpix_cxx truth harness (see hpref.cc, ../PROVENANCE.md).
#   make        build ./hpref (gitignored build artifact)
#   make clean  remove it
# HEALPIX_CXX is a healpix_cxx 8.3.0 source tree (GPL-2.0, not copied
# into this repository); the default is the copy Siril vendors.
HEALPIX_CXX ?= $(error set HEALPIX_CXX to a healpix_cxx 8.3.0 source directory, e.g. Siril 1.4.0 subprojects/healpix_cxx)
CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++11 -ffp-contract=off
SRCS := healpix_base.cc healpix_tables.cc error_handling.cc geom_utils.cc pointing.cc string_utils.cc

hpref: hpref.cc $(addprefix $(HEALPIX_CXX)/,$(SRCS))
	$(CXX) $(CXXFLAGS) -I$(HEALPIX_CXX) -o hpref hpref.cc $(addprefix $(HEALPIX_CXX)/,$(SRCS)) -lm

.PHONY: clean sha
clean:
	rm -f hpref

# sha256 of the healpix_cxx sources used (recorded in ../PROVENANCE.md)
sha:
	@cd $(HEALPIX_CXX) && sha256sum $(SRCS) healpix_base.h math_utils.h
