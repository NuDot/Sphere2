name := sphere2
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
#  G4INSTALL = ../../..
  G4INSTALL=/home/taritree/software/geant4.10.03/release/share/Geant4-10.3.0/geant4make
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

ifdef ROOTSYS
  CPPFLAGS += -I$(ROOTSYS)/include
  LDFLAGS += -L$(ROOTSYS)/lib
  LOADLIBS += $(shell root-config --libs) 
endif                                      

G4BINDIR = .
