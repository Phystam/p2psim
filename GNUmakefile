#
# Makefile for Geant4 progrum
#

name := CATANA_LaBr
G4TARGET := $(name)
G4EXLIB := true
G4WORKDIR := .

ifndef G4INSTALL
#   G4INSTALL = /home/yamada/g/geant4.9.6.p04
   G4INSTALL = /usr/local/src/geant4.9.6.p04
endif


##-----------------------------------------------------------
## environment setting for ROOT output
##-----------------------------------------------------------
#
## enable ROOT output option (see src/T00EventAction.cc etc.)
CPPFLAGS += -D__OUTPUT_ROOTFILE__ -g3 -O0
#
## add compile flags and extra libs for ROOT output
CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS = $(shell $(ROOTSYS)/bin/root-config --glibs)
EXTRALIBS += $(ROOTLIBS)
#
##-----------------------------------------------------------

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
