#!/bin/bash

g++ -std=c++11 -c main.cpp libraryZ3.cpp math.cpp test.cpp isomorphism.cpp
g++ -o main main.o libraryZ3.o math.o test.o isomorphism.o -lz3 
./main