LIBS = -lcuba -lm

all: single_unresolved.exe double_unresolved.exe

%.exe : %.o ../src/PhaseSpace.o ../src/Utilities.o ../src/Tree.o ../VEGAS_interface.o
		g++  -std=c++17 -O3 -o $@ $^ $(LIBS)

%.o : %.cpp ../src/PhaseSpace.cpp ../src/Utilities.cpp ../src/Tree.cpp ../VEGAS_interface.cpp
		g++ -std=c++17 -O3 -o $@ -c $< $(LIBS)

clean:
		rm -f *.o *.exe

.PHONY: all clean
.PRECIOUS: ../src/PhaseSpace.o ../src/Utilities.o ../src/Tree.o ../VEGAS_interface.o