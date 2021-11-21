# Configure compiler settings
CXX=clang++
CXXFLAGS=-g -std=c++11 -Werror
# The header files for the program
# HFILES = $(shell ls -1 *.h)
# UnitTest++ keeps its object files in this directory.
UNITTEST_LIB = -lUnitTest++

MAINPROG = step1
SIEVEPROG = step2


all: $(MAINPROG) $(SIEVEPROG)

# This target builds your main program.
$(MAINPROG): $(MAINPROG).o
	$(CXX) $(CXXFLAGS) -o $@ $(MAINPROG).o -lgmpxx -lgmp
# This target describes how to compile a .o file from a .cpp file.
%.o: %.cpp $(HFILES)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# This target builds your main program.
$(SIEVEPROG): $(SIEVEPROG).o
	$(CXX) $(CXXFLAGS) -o $@ $(SIEVEPROG).o -lgmpxx -lgmp
# This target describes how to compile a .o file from a .cpp file.
%.o: %.cpp $(HFILES)
	$(CXX) $(CXXFLAGS) -c -o $@ $<



# This target deletes the temporary files we have built.
.PHONY: clean all

clean:
	rm -f *.o
	rm -f $(MAINPROG)
