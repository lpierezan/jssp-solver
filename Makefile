CXX		  := g++
CXX_FLAGS := -O3 -Wall -Wno-sign-compare -Wextra -std=c++17 -ggdb

BIN		:= bin
SRC		:= src
INCLUDE	:= include
LIB		:= -L lib
INC		:= -I $(INCLUDE)

SRC_FILES := $(wildcard src/*.cpp)
SOURCES_NO_MAIN := $(filter-out src/main.cpp, $(SRC_FILES))

LIBRARIES	:=
EXECUTABLE	:= main

all: $(BIN)/$(EXECUTABLE)

run: clean all
	clear
	./$(BIN)/$(EXECUTABLE)

$(BIN)/$(EXECUTABLE): $(SRC)/*.cpp
	$(CXX) $(CXX_FLAGS) $(INC) $(LIB) $^ -o $@ $(LIBRARIES)

tester: $(SOURCES_NO_MAIN)
	$(CXX) $(CXX_FLAGS) $(INC) $(LIB) $^ test/tester.cpp -o bin/tester

clean:
	-rm $(BIN)/*
