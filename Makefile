
CXX := g++

CXXFLAGS := -Wall -Wextra -std=c++17 -O3
SRC := generator.cpp expression.cpp matrix.cpp algorithm.cpp genetic.cpp
OBJ := $(SRC:.cpp=.o)
TARGET := main

dependency: $(OBJ)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


all: dependency
	$(CXX) $(CXXFLAGS) main.cpp $(OBJ) -o $(TARGET)

clean:
	rm -rf $(OBJ) $(TARGET)