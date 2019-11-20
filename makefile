# Written by Vineet Yadav
# Date : July 26 2015 

#choice for compiler
Compiler=g++

# Parent Source Directory (Change as per use case)
sourcedir=/home/vineet/sparsecheck

# Include path for header files specific to this sourcecode
# classdefinitions.h functiondefinions.h
IDIR=$(sourcedir)/include

# Path for storing executable (Specific to melvin2.jpl.nasa.gov)
Execpath=$(sourcedir)

# sourcecode
srcfiles=$(sourcedir)/src/memfunctions.cpp \
        $(sourcedir)/src/print.cpp \
        $(sourcedir)/src/sparse-sparse=dense.cpp \
        $(sourcedir)/src/sparse-sparse=sparse.cpp \
        $(sourcedir)/src/sparsework.cpp \
        $(sourcedir)/src/workdivision.cpp \
        $(sourcedir)/src/filereading_main.cpp 
# This filereading_main.cpp main function should be used if you want to read matrices from file
# comment out main.cpp
# $(sourcedir)/src/main.cpp \
        
        
# Compiler flags
CXXFLAGS+= -I$(IDIR) -Wall -ggdb -fopenmp -lm -O3

# Create Executable
all: 
	$(Compiler) $(srcfiles) $(CXXFLAGS) -o $(Execpath)/sparse.exe

# Rule to clean existing compiled file
.PHONY: clean

clean:
	rm -f $(Execpath)/sparse.exe


