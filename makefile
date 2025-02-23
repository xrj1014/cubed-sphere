.PHONY: clean 

default: cubed-sphere

cubed-sphere: cubed-sphere.cpp cmdline.h
	g++ -std=c++11 -o cubed-sphere cubed-sphere.cpp -g -Wall -O0
	
clean: 
	rm -f cubed-sphere *.o