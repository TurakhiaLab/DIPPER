# ==== Compiler ====
CXX := g++

CXXFLAGS = -O3 -std=c++17 -march=native -Wall -Wextra \
    -I$(CONDA_PREFIX)/include
LDFLAGS  = -L$(CONDA_PREFIX)/lib -llapack -llapack -lblas -lm

INCDIR  := $(CONDA_PREFIX)/include
LIBDIR  := $(CONDA_PREFIX)/lib

CXXFLAGS += -I$(INCDIR)
LDFLAGS  += -L$(LIBDIR) -Wl,-rpath,$(LIBDIR) -llapack -llapack -lblas -lm


# ==== Libraries ====
# 若你用 OpenBLAS + LAPACKE：
LIBS :=  -llapack -lm
# 若你的系統需要 pthread / gfortran：
# LIBS := -lopenblas -llapacke -lpthread -lgfortran -lm

# ==== Targets ====
TARGET := test_pmatrix
SRC := pmat.cpp test_pmatrix.cpp 
OBJ := $(SRC:.cpp=.o)

# ==== Rules ====
all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS)  -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)

.PHONY: all clean
