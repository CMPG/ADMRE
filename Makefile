# makefile for main
admre: main.cpp
	g++ -Wall -o admre main.cpp world.cpp deme.cpp individual.cpp -std=c++11
