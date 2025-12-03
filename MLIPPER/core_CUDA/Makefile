# ==== Compilers ====
NVCC := nvcc
CXX  := g++

# ==== CUDA paths ====
CUDA_HOME ?= /usr/local/cuda
CUDA_INC  := -I$(CUDA_HOME)/include
CUDA_LIB  := -L$(CUDA_HOME)/lib64
CUDA_RPATH := -Wl,-rpath,$(CUDA_HOME)/lib64

# ==== GPU arch ====
ARCH := \
   	-gencode arch=compute_75,code=sm_75 \
   	-gencode arch=compute_80,code=sm_80 \
  	-gencode arch=compute_86,code=sm_86

# ==== libpll paths ====
PLL_INC_DIR ?= /usr/local/include
PLL_LIB_DIR ?= /usr/local/lib
PLL_RPATH   := -Wl,-rpath,$(PLL_LIB_DIR)

# ==== Flags (Release) ====）
CXXFLAGS   := -O3 -std=c++17 -march=native -ffast-math -fno-omit-frame-pointer -MMD -MP
NVCCFLAGS_COMPILE := -O3 -std=c++17 $(ARCH) -rdc=true \
    -Xcompiler -fno-omit-frame-pointer \
    -Xptxas=-O3,-dlcm=ca,-v \
    -MMD -MP

NVCCFLAGS_LINK := -O3 -std=c++17 $(ARCH) -rdc=true  \
    -Xcompiler -fno-omit-frame-pointer \
    -Xptxas=-O3,-dlcm=ca,-v 

CXXWARN    := -Wall -Wextra -Wno-unused-parameter
NVCCWARN   := -Xcompiler -Wall -Xcompiler -Wextra 

PROJECT_INC := -I.
INC       := $(PROJECT_INC) -I$(PLL_INC_DIR) $(CUDA_INC)

LDFLAGS := $(CUDA_LIB) -L$(PLL_LIB_DIR) -lpll -lcudart -lm \
  -Xlinker -rpath -Xlinker $(CUDA_HOME)/lib64 \
  -Xlinker -rpath -Xlinker $(PLL_LIB_DIR)

# ==== Sources ====
SRCS_GPU := \
  core_likelihood.cu

SRCS_TEST := \
	test_tip_inner_lnl.cpp
# 	test_root_lnl.cpp
# 	test_inner_inner_lnl.cpp
# 	test_tip_inner_lnl.cpp

#   test_inner_inner_lnl.cpp
  

OBJS_GPU  := $(SRCS_GPU:.cu=.o)
OBJS_TEST := $(SRCS_TEST:.cpp=.o)
OBJS      := $(OBJS_GPU) $(OBJS_TEST)
DEPS      := $(OBJS:.o=.d)


# ==== Targets ====
TARGET_GPU := test_inner_inner_lnl

.PHONY: all clean debug sanitize run run_cpu

all: $(TARGET_GPU)

$(TARGET_GPU): $(OBJS) cuda_link.o
	$(NVCC) $(NVCCFLAGS_LINK) $(NVCCWARN) -o $@ $(OBJS) cuda_link.o $(LDFLAGS)

cuda_link.o: $(OBJS_GPU)
	$(NVCC) $(NVCCFLAGS_LINK) -dlink $(OBJS_GPU) -o $@ $(LDFLAGS)
%.o: %.cu
	$(NVCC) $(NVCCFLAGS_COMPILE) $(NVCCWARN) $(INC) -c $< -o $@

%.o: %.cpp
	$(CXX)  $(CXXFLAGS)  $(CXXWARN) $(INC) -c $< -o $@
-include $(DEPS)
# ---- Debug build ----
debug: clean
	$(MAKE) \
	  NVCCFLAGS_COMPILE="-O0 -g -G -lineinfo $(ARCH) -rdc=true -Xcompiler -fno-omit-frame-pointer -MMD -MP" \
	  NVCCFLAGS_LINK="-O0 -g -lineinfo $(ARCH) -rdc=true -Xcompiler -fno-omit-frame-pointer" \
	  CXXFLAGS="-O0 -g -std=c++17 -fno-omit-frame-pointer -MMD -MP" \
	  all

# ---- Compute Sanitizer ----
sanitize: debug
	@echo "== Running with compute-sanitizer (memcheck) =="
	CUDA_LAUNCH_BLOCKING=1 compute-sanitizer --tool memcheck ./$(TARGET_GPU)

run: all
	./$(TARGET_GPU)

# $(TARGET_CPU): $(SRCS_TEST:.cpp=.o)
# 	$(CXX) $(CXXFLAGS) $(CXXWARN) -o $@ $^ -L$(PLL_LIB_DIR) -lpll -lm $(PLL_RPATH)

# run_cpu: $(TARGET_CPU)
# 	./$(TARGET_CPU)

clean:
		rm -f $(OBJS) $(DEPS) cuda_link.o $(TARGET_GPU)
# 	rm -f $(OBJS) $(TARGET_GPU) $(TARGET_CPU)