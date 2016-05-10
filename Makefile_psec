name := sphere2
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
#  G4INSTALL = ../../..
  G4INSTALL=/code/geant4.9.6.4/share/Geant4-9.6.4/geant4make
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
