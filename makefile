CC = g++
MyO =

all: HK_1d.o 
	$(CC)  -lgslcblas -lgsl -o HK_1d  HK_1d.o -lm -lfftw3

HK_1d.o :HK_1d.c ~/Documents/UC3M/Dottorato/Equations/equations.h
	$(CC) $(MyO) -lm -lfftw3 -lgslcblas -lgsl -c HK_1d.c ~/Documents/UC3M/Dottorato/Equations/equations.h
	
conv_2d: test_conv_2d.o
	$(CC)  -lgslcblas -lgsl -o test_conv_2d  test_conv_2d.o -lm -lfftw3
	
test_conv_2d.o: test_conv_2d.c ~/Documents/UC3M/Dottorato/Equations/equations.h
	$(CC) $(MyO) -lm -lfftw3 -lgslcblas -lgsl -c test_conv_2d.c ~/Documents/UC3M/Dottorato/Equations/equations.h

.PHONY : clean

clean:
	\rm *.o
