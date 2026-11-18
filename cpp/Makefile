# Сборка: make            (g++ >= 9 / clang >= 12, Linux)
CXX      ?= g++
CXXFLAGS += -std=c++20 -O3 -march=native -Wall -Wextra -Iinclude
LDFLAGS  += -pthread

SRC := src/rst.cpp src/reader.cpp src/zlite.cpp src/main.cpp
OBJ := $(SRC:.cpp=.o)
BIN := rst_reader

all: $(BIN)

$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp include/rst.h include/reader.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(BIN)

.PHONY: all clean
