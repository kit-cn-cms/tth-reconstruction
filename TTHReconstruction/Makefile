MEDIR=../MEs/
LHAPDFDIR=/home/mildner/HEP/lib/LHAPDF/

CC=g++
LD=g++

ROOTLIBS := $(shell root-config --libs) 
ROOTC := $(shell root-config --cflags) 
ROOTLD := $(shell root-config --ldflags)

SOURCES:=$(wildcard src/*.cpp)
OBJECTS:=$(patsubst src/%.cpp,tmp/%.o,$(SOURCES))

MEINCL=-I${MEDIR}/src -I${MEDIR}/SubProcesses/P0_Sigma_sm_gg_ttxh/ -I${MEDIR}/SubProcesses/P0_Sigma_sm_gg_ttxbbx/ -I${MEDIR}/SubProcesses/P0_Sigma_sm_gg_ttxh_ttxbbx/ -I${MEDIR}/SubProcesses/P1_Sigma_HC_UFO_gg_ttxx0/ -I${MEDIR}/SubProcesses/P1_Sigma_HC_UFO_uux_ttxx0
MEOBJECTS=${MEDIR}/SubProcesses/P0_Sigma_sm_gg_ttxh/tthProcess.o ${MEDIR}/SubProcesses/P0_Sigma_sm_gg_ttxbbx/ttbbProcess.o ${MEDIR}/SubProcesses/P0_Sigma_sm_gg_ttxh_ttxbbx/tthbbProcess.o ${MEDIR}/SubProcesses/P1_Sigma_HC_UFO_gg_ttxx0/ggttx0Process.o ${MEDIR}/SubProcesses/P1_Sigma_HC_UFO_uux_ttxx0/qqttx0Process.o
MELIBS=-L${MEDIR}/lib -lmodel_sm
LHAPDFINCL=-I${LHAPDFDIR}/include/
LHAPDFLIB=-L${LHAPDFDIR}/lib -l LHAPDF

all : bin/test

tmp/%.o : src/%.cpp
	$(CC) $(ROOTC) -fPIC -I./interface $(MEINCL) $(LHAPDFINCL) -c -o $@ $<

bin/test : $(OBJECTS)
	$(LD) $(ROOTLD) -o bin/test tmp/*o $(MEOBJECTS) $(ROOTLIBS) $(MELIBS) $(LHAPDFLIB)

clean :
	rm -f tmp/*; rm -f bin/*
