CXX=g++
CXXFLAGS=-Wall -O3 -pedantic
LIBS= -lblas -fopenmp
TARGET=main
OBJS=main.o ShallowWater.o
HDRS=ShallowWater.h

default: $(TARGET)

%.o:%.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIBS)

$(TARGET):$(OBJS)
	$(CXX) -o $@ $^ $(LIBS)

export OMP_NUM_THREADS := 10

test1: $(TARGET)
	./$(TARGET) --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 1

test2: $(TARGET)
	./$(TARGET) --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 2

test3: $(TARGET)
	./$(TARGET) --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 3

test4: $(TARGET)
	./$(TARGET) --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 4

.PHONY: clean plot

clean:
	rm -rf $(TARGET) *.o

plot:
	gnuplot
