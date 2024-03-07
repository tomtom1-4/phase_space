main.exe : main.o src/PhaseSpace.o src/Utilities.o
		g++  -std=c++17 -O3 -o $@ $^ $(LIBS)

main.o : main.cpp src/PhaseSpace.cpp src/Utilities.cpp
		g++ -std=c++17 -O3 -o $@ -c $< $(LIBS)

clean:
		rm -f *.o *.exe

.PHONY: all clean
.PRECIOUS: src/PhaseSpace.o src/Utilities.o