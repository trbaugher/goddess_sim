
ROOTCFLAGS := $(shell root-config --cflags) 
ROOTLIBS := $(shell root-config --libs)  
ROOTGLIBS := $(shell root-config --glibs) -lSpectrum -lThread -lpthread
ROOTINC := -I$(shell root-config --incdir)

CXX = g++ -Wall -g -O3 -ggdb 
CPPFLAGS = -fPIC $(ROOTCFLAGS)

.PHONY : all clean

all : godsort

#GodHit.o : GodHit.h
#$(CXX) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

#godsort : godsort.cpp GodHit.o
#	@echo "BUILDING MAIN APPLICATION"
#	@echo "dependencies: $^"
#	$(CXX) $(CPPFLAGS) $(ROOTLIBS) $(ROOTGLIBS) GodHit.o $<  -o $@
#	@echo "done"
godsort : godsort.cpp
	@echo "BUILDING MAIN APPLICATION"
	@echo "dependencies: $^"
	$(CXX) $(CPPFLAGS) $(ROOTLIBS) $(ROOTGLIBS) $<  -o ../$@
	@echo "done"

clean :
	rm -f *.o

