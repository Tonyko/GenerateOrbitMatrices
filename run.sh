#!/bin/bash

g++ -std=c++11 -c main.cpp libraryZ3.cpp
g++ -o main main.o libraryZ3.o -lz3
./main