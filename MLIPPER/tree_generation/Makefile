# ==========================================
# CUDA
# ==========================================
CUDA_HOME ?= /usr/local/cuda-12
NVCC      := $(CUDA_HOME)/bin/nvcc
CUDA_INC  := -I$(CUDA_HOME)/include
CUDA_LIB  := -L$(CUDA_HOME)/lib64
NVCC_ARCH ?= -gencode arch=compute_86,code=sm_86

# ==========================================
# Compiler
# ==========================================
CXX := g++

# ==========================================
# libpll
# ==========================================
PLL_INC_DIR ?= /usr/local/include
PLL_LIB_DIR ?= /usr/local/lib

# ==========================================
# System LAPACK / BLAS
# ==========================================
LAPACK_LIBS := -L/usr/lib/x86_64-linux-gnu \
               -llapack -lblas -lm

# ==========================================
# Include paths（這裡是重點）
# ==========================================
INCLUDES := -I. \
            -I../pmatrix \
            -I../partial_CUDA \
            -I../core_CUDA \
            -I$(PLL_INC_DIR) \
            -I$(CUDA_HOME)/include \
            -I$(CONDA_PREFIX)/include

# ==========================================
# Flags
# ==========================================
CXXFLAGS := -O3 -std=c++17 -fno-omit-frame-pointer $(INCLUDES)

NVCCFLAGS := -O3 -std=c++17 -ccbin $(CXX) \
             -Xcompiler -fno-omit-frame-pointer \
             $(INCLUDES) $(CUDA_INC) $(NVCC_ARCH)

# ==========================================
# Sources
# ==========================================
CPU_SRCS  := ../pmatrix/pmat.cpp \
             tree_generation.cpp \
             tree.cpp \
             parse_file.cpp

CUDA_SRCS := tree_generation_device.cu \
             ../partial_CUDA/partial_likelihood.cu \
             root_likelihood.cu

CPU_OBJS  := $(CPU_SRCS:.cpp=.o)
CUDA_OBJS := $(CUDA_SRCS:.cu=.cu.o)
OBJS      := $(CPU_OBJS) $(CUDA_OBJS)

TARGET := test_tree_gen

# ==========================================
# Rules
# ==========================================
# ---- C++ -> .o ----
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ---- CUDA -> .cu.o ----
%.cu.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

# ---- Link ----
$(TARGET): $(OBJS)
	$(NVCC) -ccbin $(CXX) -o $@ $^ \
		-L$(PLL_LIB_DIR) -lpll \
		$(CUDA_LIB) -lcudart -lcuda \
		$(LAPACK_LIBS)

clean:
	rm -f $(CPU_OBJS) $(CUDA_OBJS) $(TARGET)
