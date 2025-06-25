
# optimization level, default is 3
LEVEL = 3

USE_CUBA = true # true/false interface to CUBA MC integration using the VEGAS algorithm

LIBS = -L/home/tom/Uni/phd/NnLO/phase_space_parameterization -lphase_space -lm
CXXFLAGS = -std=c++17 -O$(LEVEL)

ifeq ($(strip $(USE_CUBA)),true)
	LIBS += -lcuba
	CXXFLAGS += -DUSE_CUBA
endif

export USE_CUBA CXXFLAGS LIBS

all: phase_space.exe


# Build main.exe; only link in CUBA-related objects if USE_CUBA=true
%.exe: %.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

# Generic rule to build .o files from .cpp files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean:
		rm -f *.o *.exe

.PHONY: all clean