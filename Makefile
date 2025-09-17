# optimization level, default is 3
LEVEL = 3

USE_CUBA = true # true/false interface to CUBA MC integration using the VEGAS algorithm

LIBS = -lm
CXXFLAGS = -std=c++17 -O$(LEVEL) -fPIC

ifeq ($(strip $(USE_CUBA)),true)
	LIBS += -lcuba
	CXXFLAGS += -DUSE_CUBA -I$(CUBA_PATH)
	CUBA_SOURCES = VEGAS_interface.cpp
	CUBA_OBJECTS = VEGAS_interface.o
else
	CUBA_SOURCES =
	CUBA_OBJECTS =
endif

export USE_CUBA CXXFLAGS LIBS

SOURCES = src/PhaseSpace.cpp src/Utilities.cpp src/Tree.cpp
OBJECTS = $(SOURCES:.cpp=.o)

# Add CUBA objects if USE_CUBA=true
OBJECTS += $(CUBA_OBJECTS)

all: libphase_space.so

# Build a shared library libmylib.so from the object files in src/
%.so: $(OBJECTS)
	$(CXX) -shared -o $@ $^

# Build main.exe; only link in CUBA-related objects if USE_CUBA=true
%.exe: %.o $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

# Generic rule to build .o files from .cpp files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Build object files for sources in src/
src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o *.exe src/*.o

tests:
	$(MAKE) -C tests

.PRECIOUS: src/PhaseSpace.o src/Utilities.o src/Tree.o VEGAS_interface.o
.PHONY: all clean tests
