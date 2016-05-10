#!/bin/bash

export ROOTSYS=/code/root_v5.34.01_vnc
source $ROOTSYS/bin/thisroot.sh
#export LD_LIBRARY_PATH=$ROOTSYS/lib:$ROOTSYS/include:$LD_LIBRARY_PATH

source /code/geant4.9.6.4/bin/geant4.sh
source /code/geant4.9.6.4/share/Geant4-9.6.4/geant4make/geant4make.sh

export G4INSTALL=/code/geant4.9.6.4/share/Geant4-9.6.4/geant4make
#export G4LEDATA=/local/data1/elagin/0vbb/geant4/data/G4EMLOW6.32
export G4WORKDIR=/local/data2/elagin/g4w

