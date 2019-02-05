CC = g++
LFLAGS = -lgsl -lgslcblas -lm
LIBS = -lgsl

all: sam_LG_MC.cpp constants.h
	$(CC) -g sam_LG_MC.cpp -o do_mc.o $(LFLAGS) $(LIBS)

clean:
	rm *.o 
