
CXX := g++

CXXFLAGS := -Wall -Wextra -std=c++17 -O3
SRC := generator.cpp expression.cpp matrix.cpp algorithm.cpp genetic.cpp csv.cpp infra.cpp info.cpp
OBJ := $(SRC:.cpp=.o)
TARGET := main

dependency: $(OBJ)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


all: dependency
	$(CXX) $(CXXFLAGS) main.cpp $(OBJ) -o $(TARGET)

matrix: dependency
	$(CXX) $(CXXFLAGS) matrix_calculation.cpp $(OBJ) -o matrix

clean:
	rm -rf $(OBJ) $(TARGET) matrix