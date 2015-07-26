MEDIR=../MEs/

CC=g++
LD=g++

ROOTLIBS := $(shell root-config --libs) 
ROOTC := $(shell root-config --cflags) 
ROOTLD := $(shell root-config --ldflags)

SOURCES:=$(wildcard src/*.cpp)
OBJECTS:=$(patsubst src/%.cpp,tmp/%.o,$(SOURCES))

MEINCL=-I${MEDIR}/src -I${MEDIR}/SubProcesses/P0_Sigma_sm_gg_ttxh/ -I${MEDIR}/SubProcesses/P0_Sigma_sm_gg_ttxbbx/ -I${MEDIR}/SubProcesses/P0_Sigma_sm_gg_ttxh_ttxbbx/
MEOBJECTS=${MEDIR}/SubProcesses/P0_Sigma_sm_gg_ttxh/tthProcess.o ${MEDIR}/SubProcesses/P0_Sigma_sm_gg_ttxbbx/ttbbProcess.o ${MEDIR}/SubProcesses/P0_Sigma_sm_gg_ttxh_ttxbbx/tthbbProcess.o
MELIBS=-L${MEDIR}/lib -lmodel_sm

all : bin/test

tmp/%.o : src/%.cpp
	$(CC) $(ROOTC) -fPIC -I./interface $(MEINCL) -c -o $@ $<

bin/test : $(OBJECTS)
	$(LD) $(ROOTLD) -o bin/test tmp/*o $(MEOBJECTS) $(ROOTLIBS) $(MELIBS)

clean :
	rm -f tmp/*; rm -f bin/*
