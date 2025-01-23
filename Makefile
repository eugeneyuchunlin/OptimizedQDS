
CXX := g++

CXXFLAGS := -Wall -Wextra -std=c++17 -g -O0
SRC := generator.cpp expression.cpp matrix.cpp
OBJ := $(SRC:.cpp=.o)
TARGET := main

dependency: $(OBJ)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


all: dependency
	$(CXX) $(CXXFLAGS) main.cpp $(OBJ) -o $(TARGET)

clean:
	rm -rf $(OBJ) $(TARGET)