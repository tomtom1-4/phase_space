LIBS = -lcuba -lm

main.exe : main.o src/PhaseSpace.o src/Utilities.o src/Tree.o
		g++  -std=c++17 -O3 -o $@ $^ $(LIBS)

main.o : main.cpp src/PhaseSpace.cpp src/Utilities.cpp src/Tree.cpp
		g++ -std=c++17 -O3 -o $@ -c $< $(LIBS)

clean:
		rm -f *.o *.exe

.PHONY: all clean
.PRECIOUS: src/PhaseSpace.o src/Utilities.o src/Tree.o