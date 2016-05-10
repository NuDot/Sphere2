name := sphere2
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = /u/nobackup/lwinslow/apps/geant4/geant4.9.6.p02-install/share/Geant4-9.6.2/geant4make
#  G4INSTALL=/code/geant4.9.5_rhel6/share/Geant4-9.5.1/geant4make
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

ifdef ROOTSYS
  CPPFLAGS += -I$(ROOTSYS)/include/root
  LDFLAGS += -L$(ROOTSYS)/lib/root
  LOADLIBS += $(shell root-config --libs) 
endif                                      

G4BINDIR = .
