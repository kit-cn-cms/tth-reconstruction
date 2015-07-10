CC= g++
LD= g++

ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LDFLAGS := $(shell root-config --ldflags)
ROOT_LIBS := $(shell root-config --libs) 
ROOT_INCDIR := $(shell root-config --incdir) 

INCLUDES = -I.. -I. -I./interface -I$(ROOT_INCDIR)
CFLAGS = -c -fPIC $(ROOT_CFLAGS) -g
LDFLAGS = $(ROOT_LDFLAGS) $(ROOT_LIBS)

SOURCES:=$(wildcard src/*.cpp)
OBJECTS:=$(patsubst src/%.cpp,tmp/%.o,$(SOURCES))

all : bin/test

tmp/%.o : src/%.cpp
	@echo creating $@ from $<; $(CC) $(CFLAGS) $(LDFLAGS) $(INCLUDES) -c -o $@ $<

bin/test : $(OBJECTS)
	@echo creating $@ from $+; $(LD) $(LDFLAGS) $(INCLUDES) -o $@ $+

clean :
	rm -f tmp/*; rm -f bin/*
