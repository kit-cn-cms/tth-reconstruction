CC=g++
LD=g++

ROOTLIBS := $(shell root-config --libs) 
ROOTC := $(shell root-config --cflags) 
ROOTLD := $(shell root-config --ldflags)

SOURCES:=$(wildcard src/*.cpp)
OBJECTS:=$(patsubst src/%.cpp,tmp/%.o,$(SOURCES))

all : bin/test

tmp/%.o : src/%.cpp
	$(CC) $(ROOTC) -fPIC -I./interface -c -o $@ $<

bin/test : $(OBJECTS)
	$(LD) $(ROOTLD) -o bin/test tmp/*o $(ROOTLIBS)

clean :
	rm -f tmp/*; rm -f bin/*
