CXX = g++-8
HOME= /usr/local/include
LIB_HOME = flowstar-template
LIBS = -lflowstar -lmpfr -lgmp -lgsl -lgslcblas -lm -lglpk
CFLAGS = -I . -I $(HOME) -g -O3 -std=c++11
LINK_FLAGS = -g -L$(LIB_HOME) -L/usr/local/lib
OBJS = NeuralNetwork.o Activation.o Neuron.o NNTaylor.o domain_computation.o

all: runtime test

runtime: $(OBJS)
	g++-8 -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)

test: test_new_implementation.o
	g++-8 -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)

%.o: %.cc
	$(CXX) -O3 -c $(CFLAGS) -o $@ $<
%.o: %.cpp
	$(CXX) -O3 -c $(CFLAGS) -o $@ $<
%.o: %.c
	$(CXX) -O3 -c $(CFLAGS) -o $@ $<


clean:
	rm -f *.o

