# Configure compiler settings
CC = mpicc
CFLAGS = -g -Wall -Wvla -Werror -Wno-error=unused-variable -Wno-error=unused-but-set-variable

LIBDIR = -L./
INCDIR = -I./


MAINPROG = step1
SIEVEPROG = step2
GAUSSPROG = step3


all: $(MAINPROG) $(SIEVEPROG) $(GAUSSPROG)

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

# This target builds your main program.
$(GAUSSPROG): $(GAUSSPROG).o
	$(CXX) $(CXXFLAGS) -o $@ $(GAUSSPROG).o -lgmpxx -lgmp
# This target describes how to compile a .o file from a .cpp file.
%.o: %.cpp $(HFILES)
	$(CXX) $(CXXFLAGS) -c -o $@ $<



# This target deletes the temporary files we have built.
.PHONY: clean all

clean:
	rm -f *.o
	rm -f *.txt
	rm -f $(MAINPROG)
	rm -f $(SIEVEPROG)
