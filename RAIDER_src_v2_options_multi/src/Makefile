# Nathan Figueroa

CC = g++ -std=c++11

CXXFLAGS+=-Wall -W -O3 -DNDEBUG
CFLAGS = -Wall -c -DNDEBUG

CPPINCL+=-I../../seqan/include
CPPINCL+=-I../../seqan/extras/include

###########################################################
# RAIDER
###########################################################
all: main

main: Family.o LmerVector.o main.o
	@echo compiling RAIDER;
	$(CC) $(CPPINCL) $(CXXFLAGS) $^ -o ../../raiderv2_options

main.o: main.cpp SeedChain.h
	$(CC) $(CPPINCL) $(CFLAGS) main.cpp

Family.o: Family.cpp Family.h LmerVector.h
	$(CC) $(CFLAGS) Family.cpp

LmerVector.o: LmerVector.cpp LmerVector.h
	$(CC) $(CFLAGS) LmerVector.cpp
