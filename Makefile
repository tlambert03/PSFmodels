# @(#)makefile
# after conda install gcc libtiff

ifeq ($(shell uname),Darwin)
  PLATFORM = darwin64
endif
ifeq ($(shell uname),Linux)
  PLATFORM = linux64
endif

CC = gcc

CFLAGS      = -I$(CONDA_PREFIX)/include -Wall 
LDFLAGS = ""
ifeq ($(shell uname),Darwin)
    LDFLAGS += -L$(CONDA_PREFIX)/lib/
    LDFLAGS += "-Wl,-rpath,$(CONDA_PREFIX)/lib/"
	CXXFLAGS += -I/opt/local/include
endif


LIBS = -lboost_program_options -ltiff -lX11 -lpthread

all: psf

psf: main.cpp
	$(CXX) -o psf $^ $(LDFLAGS) $(CFLAGS) $(LIBS)

vectorialpsf: vectorialpsf.cpp
  c++ -O3 -Wall -shared -std=c++11 -fPIC -Wl,-undefined,dynamic_lookup `python3 -m pybind11 --includes` $^ -o vectorialpsf`python3-config --extension-suffix`

clean:
	$(RM) *.o *.exe *~ psf *.so
