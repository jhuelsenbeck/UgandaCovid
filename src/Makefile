# =============================================================================
# Makefile
# =============================================================================
# Phylogeographic CTMC inference for large-scale viral evolution
# 
# Features:
#   - Automatic platform detection (macOS, Linux, Windows/MinGW)
#   - Multiple BLAS/LAPACK backends (Accelerate, OpenBLAS, MKL, system)
#   - Debug and Release build configurations
#   - Proper dependency tracking
#   - Parallel build support
#
# Usage:
#   make                  # Build with auto-detected settings
#   make release          # Optimized release build
#   make debug            # Debug build with symbols
#   make BLAS=openblas    # Force OpenBLAS backend
#   make BLAS=mkl         # Use Intel MKL
#   make clean            # Remove build artifacts
#   make info             # Show detected configuration
#
# =============================================================================

# -----------------------------------------------------------------------------
# Project Configuration
# -----------------------------------------------------------------------------

PROGRAM      = cov
VERSION      = 2.0.0

# Source files
SOURCES      = main.cpp \
               CondLikeJob.cpp \
               CondLikeJobMngr.cpp \
               Container.cpp \
               FileManager.cpp \
               GPUMatrixExponentialBatch.cpp \
               History.cpp \
               HistorySummary.cpp \
               MathCache.cpp \
               MathCacheAccelerated.cpp \
               Mcmc.cpp \
               McmcInfo.cpp \
               McmcPhase.cpp \
               MetaData.cpp \
               Model.cpp \
               Msg.cpp \
               Node.cpp \
               Probability.cpp \
               RandomVariable.cpp \
               RateMatrix.cpp \
               Samples.cpp \
               SteppingStones.cpp \
               Stone.cpp \
               Threads.cpp \
               TransitionProbabilities.cpp \
               TransitionProbabilitiesMngr.cpp \
               Tree.cpp \
               UserSettings.cpp

OBJECTS      = $(SOURCES:.cpp=.o)
DEPENDS      = $(SOURCES:.cpp=.d)

# -----------------------------------------------------------------------------
# Platform Detection
# -----------------------------------------------------------------------------

UNAME_S := $(shell uname -s 2>/dev/null || echo Windows)
UNAME_M := $(shell uname -m 2>/dev/null || echo x86_64)

# Detect if running on Apple Silicon and which chip variant
ifeq ($(UNAME_S),Darwin)
    ifeq ($(UNAME_M),arm64)
        APPLE_SILICON = 1
        # Extract chip variant (M1, M2, M3, M4, M5, etc.) from CPU brand string
        APPLE_CHIP := $(shell sysctl -n machdep.cpu.brand_string 2>/dev/null | grep -oE 'M[0-9]+' | head -1)
    endif
endif

# -----------------------------------------------------------------------------
# Compiler Selection
# -----------------------------------------------------------------------------

# Default compilers by platform
ifeq ($(UNAME_S),Darwin)
    CXX ?= clang++
else ifeq ($(UNAME_S),Linux)
    CXX ?= g++
else
    # Windows/MinGW
    CXX ?= g++
endif

# Check for available compilers if default not found
ifeq ($(shell which $(CXX) 2>/dev/null),)
    ifeq ($(shell which clang++ 2>/dev/null),)
        CXX = g++
    else
        CXX = clang++
    endif
endif

# -----------------------------------------------------------------------------
# BLAS/LAPACK Backend Selection
# -----------------------------------------------------------------------------

# Auto-detect BLAS backend if not specified
ifndef BLAS
    ifeq ($(UNAME_S),Darwin)
        # macOS: Always use Accelerate
        BLAS = accelerate
    else ifeq ($(UNAME_S),Linux)
        # Linux: Try to detect available BLAS
        ifneq ($(shell pkg-config --exists openblas 2>/dev/null && echo yes),)
            BLAS = openblas
        else ifneq ($(shell ldconfig -p 2>/dev/null | grep -c libopenblas),0)
            BLAS = openblas
        else ifneq ($(shell pkg-config --exists mkl 2>/dev/null && echo yes),)
            BLAS = mkl
        else ifneq ($(shell ldconfig -p 2>/dev/null | grep -c libblas),0)
            BLAS = system
        else
            BLAS = none
        endif
    else
        BLAS = none
    endif
endif

# -----------------------------------------------------------------------------
# Build Configuration
# -----------------------------------------------------------------------------

# Default to release build
BUILD ?= release

# Base flags for all builds
CXXFLAGS_BASE = -std=c++20 -Wall -Wextra -Wpedantic

# Release flags
CXXFLAGS_RELEASE = -O3 -DNDEBUG

# Debug flags  
CXXFLAGS_DEBUG = -g -O0 -DDEBUG -fsanitize=address -fsanitize=undefined

# Profile flags (release + debug symbols)
CXXFLAGS_PROFILE = -O3 -g -DNDEBUG

# Architecture-specific optimizations
ifeq ($(UNAME_S),Darwin)
    ifeq ($(APPLE_SILICON),1)
        # Apple Silicon - select appropriate -mcpu flag based on chip
        ifeq ($(APPLE_CHIP),M1)
            CXXFLAGS_ARCH = -mcpu=apple-m1
        else ifeq ($(APPLE_CHIP),M2)
            CXXFLAGS_ARCH = -mcpu=apple-m2
        else ifeq ($(APPLE_CHIP),M3)
            CXXFLAGS_ARCH = -mcpu=apple-m3
        else ifeq ($(APPLE_CHIP),M4)
            CXXFLAGS_ARCH = -mcpu=apple-m4
        else ifeq ($(APPLE_CHIP),M5)
            CXXFLAGS_ARCH = -mcpu=apple-m5
        else
            # Future chip or detection failed - use native (safe fallback)
            CXXFLAGS_ARCH = -mcpu=native
        endif
    else
        # Intel Mac
        CXXFLAGS_ARCH = -march=native
    endif
else
    # Linux/Windows - use native arch
    CXXFLAGS_ARCH = -march=native
endif

# Link-time optimization (release only)
ifeq ($(BUILD),release)
    CXXFLAGS_LTO = -flto
    LDFLAGS_LTO = -flto
else ifeq ($(BUILD),profile)
    CXXFLAGS_LTO = -flto
    LDFLAGS_LTO = -flto
endif

# Fast math (can be disabled with FASTMATH=0)
FASTMATH ?= 1
ifeq ($(FASTMATH),1)
    CXXFLAGS_MATH = -ffast-math
endif

# -----------------------------------------------------------------------------
# BLAS/LAPACK Configuration
# -----------------------------------------------------------------------------

ifeq ($(BLAS),accelerate)
    # Apple Accelerate Framework
    CXXFLAGS_BLAS = -DUSE_ACCELERATE
    LDFLAGS_BLAS = -framework Accelerate
    BLAS_NAME = Apple Accelerate

else ifeq ($(BLAS),openblas)
    # OpenBLAS
    OPENBLAS_INC ?= $(shell pkg-config --cflags openblas 2>/dev/null || echo -I/usr/include/openblas)
    OPENBLAS_LIB ?= $(shell pkg-config --libs openblas 2>/dev/null || echo -lopenblas)
    CXXFLAGS_BLAS = -DUSE_OPENBLAS -DUSE_CBLAS $(OPENBLAS_INC)
    LDFLAGS_BLAS = $(OPENBLAS_LIB)
    BLAS_NAME = OpenBLAS

else ifeq ($(BLAS),mkl)
    # Intel MKL
    MKL_ROOT ?= /opt/intel/oneapi/mkl/latest
    CXXFLAGS_BLAS = -DUSE_MKL -DUSE_CBLAS -I$(MKL_ROOT)/include
    # Link with MKL (sequential for thread pool compatibility)
    LDFLAGS_BLAS = -L$(MKL_ROOT)/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
    BLAS_NAME = Intel MKL

else ifeq ($(BLAS),system)
    # System BLAS/LAPACK
    CXXFLAGS_BLAS = -DUSE_CBLAS
    LDFLAGS_BLAS = -lblas -llapack
    BLAS_NAME = System BLAS

else
    # No BLAS - portable fallback
    CXXFLAGS_BLAS =
    LDFLAGS_BLAS =
    BLAS_NAME = Portable (no BLAS)
endif

# -----------------------------------------------------------------------------
# Threading Configuration
# -----------------------------------------------------------------------------

ifeq ($(UNAME_S),Darwin)
    # macOS uses Grand Central Dispatch internally, but we still need pthreads
    LDFLAGS_THREAD = -lpthread
else ifeq ($(UNAME_S),Linux)
    LDFLAGS_THREAD = -lpthread
else
    # Windows/MinGW
    LDFLAGS_THREAD = -lpthread
endif

# -----------------------------------------------------------------------------
# Assemble Final Flags
# -----------------------------------------------------------------------------

ifeq ($(BUILD),debug)
    CXXFLAGS = $(CXXFLAGS_BASE) $(CXXFLAGS_DEBUG) $(CXXFLAGS_BLAS)
    LDFLAGS = $(LDFLAGS_BLAS) $(LDFLAGS_THREAD) -fsanitize=address -fsanitize=undefined
else ifeq ($(BUILD),profile)
    CXXFLAGS = $(CXXFLAGS_BASE) $(CXXFLAGS_PROFILE) $(CXXFLAGS_ARCH) $(CXXFLAGS_LTO) $(CXXFLAGS_MATH) $(CXXFLAGS_BLAS)
    LDFLAGS = $(LDFLAGS_LTO) $(LDFLAGS_BLAS) $(LDFLAGS_THREAD)
else
    # Release build
    CXXFLAGS = $(CXXFLAGS_BASE) $(CXXFLAGS_RELEASE) $(CXXFLAGS_ARCH) $(CXXFLAGS_LTO) $(CXXFLAGS_MATH) $(CXXFLAGS_BLAS)
    LDFLAGS = $(LDFLAGS_LTO) $(LDFLAGS_BLAS) $(LDFLAGS_THREAD)
endif

# Add math library
LDFLAGS += -lm

# -----------------------------------------------------------------------------
# Build Targets
# -----------------------------------------------------------------------------

.PHONY: all release debug profile clean distclean info help install

# Default target
all: $(PROGRAM)

# Release build
release:
	@$(MAKE) BUILD=release $(PROGRAM)

# Debug build
debug:
	@$(MAKE) BUILD=debug $(PROGRAM)

# Profile build (optimized with debug symbols)
profile:
	@$(MAKE) BUILD=profile $(PROGRAM)

# Link the program
$(PROGRAM): $(OBJECTS)
	@echo "Linking $(PROGRAM)..."
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $(PROGRAM)
	@echo "Build complete: $(PROGRAM)"

# Compile source files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

# Include dependency files
-include $(DEPENDS)

# -----------------------------------------------------------------------------
# Cleaning
# -----------------------------------------------------------------------------

clean:
	@echo "Cleaning build artifacts..."
	rm -f $(OBJECTS) $(DEPENDS) $(PROGRAM)
	@echo "Clean complete."

distclean: clean
	rm -f *.orig *.bak *~ core
	rm -rf *.dSYM

# -----------------------------------------------------------------------------
# Installation
# -----------------------------------------------------------------------------

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin

install: $(PROGRAM)
	@echo "Installing $(PROGRAM) to $(BINDIR)..."
	install -d $(BINDIR)
	install -m 755 $(PROGRAM) $(BINDIR)/$(PROGRAM)
	@echo "Installation complete."

uninstall:
	@echo "Removing $(PROGRAM) from $(BINDIR)..."
	rm -f $(BINDIR)/$(PROGRAM)

# -----------------------------------------------------------------------------
# Information
# -----------------------------------------------------------------------------

info:
	@echo "=============================================="
	@echo "Uganda Covid $(VERSION) Build Configuration"
	@echo "=============================================="
	@echo ""
	@echo "Platform:"
	@echo "  OS:           $(UNAME_S)"
	@echo "  Architecture: $(UNAME_M)"
ifdef APPLE_SILICON
	@echo "  Apple Silicon: Yes ($(APPLE_CHIP))"
endif
	@echo ""
	@echo "Compiler:"
	@echo "  CXX:          $(CXX)"
	@echo "  Version:      $(shell $(CXX) --version | head -1)"
	@echo ""
	@echo "BLAS Backend:"
	@echo "  Selected:     $(BLAS_NAME)"
	@echo "  CXXFLAGS:     $(CXXFLAGS_BLAS)"
	@echo "  LDFLAGS:      $(LDFLAGS_BLAS)"
	@echo ""
	@echo "Build Type:     $(BUILD)"
	@echo ""
	@echo "Flags:"
	@echo "  CXXFLAGS:     $(CXXFLAGS)"
	@echo "  LDFLAGS:      $(LDFLAGS)"
	@echo ""
	@echo "Source Files:   $(words $(SOURCES))"
	@echo "=============================================="

help:
	@echo "Uganda Covid $(VERSION) - Build System"
	@echo ""
	@echo "Usage: make [target] [options]"
	@echo ""
	@echo "Targets:"
	@echo "  all (default)  Build the program"
	@echo "  release        Build optimized release version"
	@echo "  debug          Build with debug symbols and sanitizers"
	@echo "  profile        Build optimized with debug symbols"
	@echo "  clean          Remove build artifacts"
	@echo "  distclean      Remove all generated files"
	@echo "  install        Install to PREFIX (default: /usr/local)"
	@echo "  uninstall      Remove installed files"
	@echo "  info           Show build configuration"
	@echo "  help           Show this help message"
	@echo ""
	@echo "Options:"
	@echo "  BUILD=release|debug|profile   Build configuration"
	@echo "  BLAS=accelerate|openblas|mkl|system|none"
	@echo "                                BLAS backend selection"
	@echo "  CXX=compiler                  C++ compiler to use"
	@echo "  FASTMATH=0|1                  Enable fast-math (default: 1)"
	@echo "  PREFIX=/path                  Installation prefix"
	@echo ""
	@echo "Examples:"
	@echo "  make                          # Auto-detect and build"
	@echo "  make -j8                      # Parallel build with 8 jobs"
	@echo "  make debug                    # Debug build"
	@echo "  make BLAS=openblas            # Force OpenBLAS"
	@echo "  make CXX=g++-13 release       # Use specific compiler"
	@echo "  make PREFIX=~/local install   # Install to home directory"

# -----------------------------------------------------------------------------
# Header Dependencies (explicit for critical headers)
# -----------------------------------------------------------------------------

# Core dependencies
Container.o: Container.cpp Container.hpp
Node.o: Node.cpp Node.hpp
Tree.o: Tree.cpp Tree.hpp Node.hpp
Threads.o: Threads.cpp Threads.hpp

# Math/BLAS dependencies
MathCache.o: MathCache.cpp MathCache.hpp Container.hpp
MathCacheAccelerated.o: MathCacheAccelerated.cpp MathCacheAccelerated.hpp MathCache.hpp Container.hpp
GPUMatrixExponentialBatch.o: GPUMatrixExponentialBatch.cpp GPUMatrixExponentialBatch.hpp MathCacheAccelerated.hpp Container.hpp Threads.hpp

# Transition probability dependencies
TransitionProbabilities.o: TransitionProbabilities.cpp TransitionProbabilities.hpp MathCacheAccelerated.hpp Container.hpp
TransitionProbabilitiesMngr.o: TransitionProbabilitiesMngr.cpp TransitionProbabilitiesMngr.hpp GPUMatrixExponentialBatch.hpp TransitionProbabilities.hpp

# Conditional likelihood dependencies
CondLikeJob.o: CondLikeJob.cpp CondLikeJob.hpp Node.hpp TransitionProbabilities.hpp Threads.hpp
CondLikeJobMngr.o: CondLikeJobMngr.cpp CondLikeJobMngr.hpp CondLikeJob.hpp Tree.hpp Node.hpp

# Model dependencies
Model.o: Model.cpp Model.hpp Tree.hpp Node.hpp RateMatrix.hpp TransitionProbabilitiesMngr.hpp CondLikeJobMngr.hpp
RateMatrix.o: RateMatrix.cpp RateMatrix.hpp Container.hpp MathCache.hpp

# MCMC dependencies
Mcmc.o: Mcmc.cpp Mcmc.hpp Model.hpp
McmcInfo.o: McmcInfo.cpp McmcInfo.hpp

# Metadata dependencies
MetaData.o: MetaData.cpp MetaData.hpp Node.hpp Tree.hpp

# Utility dependencies
Probability.o: Probability.cpp Probability.hpp
RandomVariable.o: RandomVariable.cpp RandomVariable.hpp
Msg.o: Msg.cpp Msg.hpp
UserSettings.o: UserSettings.cpp UserSettings.hpp

# History dependencies
History.o: History.cpp History.hpp
HistorySummary.o: HistorySummary.cpp HistorySummary.hpp
Samples.o: Samples.cpp Samples.hpp

# Main
main.o: main.cpp Model.hpp Mcmc.hpp UserSettings.hpp

# -----------------------------------------------------------------------------
# Special Targets
# -----------------------------------------------------------------------------

# Rebuild everything
rebuild: clean all

# Check code with static analyzer (clang only)
analyze:
	@echo "Running static analysis..."
	$(CXX) --analyze $(CXXFLAGS_BASE) $(CXXFLAGS_BLAS) $(SOURCES)

# Format code (requires clang-format)
format:
	@echo "Formatting source code..."
	clang-format -i $(SOURCES) *.hpp

# Count lines of code
loc:
	@echo "Lines of code:"
	@wc -l $(SOURCES) *.hpp | tail -1

# Generate documentation (requires doxygen)
docs:
	@echo "Generating documentation..."
	doxygen Doxyfile

# Run the program (for testing)
run: $(PROGRAM)
	./$(PROGRAM)

# Valgrind memory check (Linux only)
memcheck: debug
	valgrind --leak-check=full --show-leak-kinds=all ./$(PROGRAM)

# Perf profiling (Linux only)  
perf: profile
	perf record -g ./$(PROGRAM)
	perf report

# -----------------------------------------------------------------------------
# Platform-Specific Metal Support (macOS only)
# -----------------------------------------------------------------------------

ifeq ($(UNAME_S),Darwin)
# Metal compute shader compilation (future GPU support)
METAL_SOURCES = $(wildcard *.metal)
METAL_LIBS = $(METAL_SOURCES:.metal=.metallib)

%.metallib: %.metal
	xcrun -sdk macosx metal -c $< -o $(@:.metallib=.air)
	xcrun -sdk macosx metallib $(@:.metallib=.air) -o $@
	rm -f $(@:.metallib=.air)

metal: $(METAL_LIBS)
endif

# -----------------------------------------------------------------------------
# Continuous Integration Helpers
# -----------------------------------------------------------------------------

# CI build (strict warnings, all checks)
ci:
	@$(MAKE) BUILD=release CXXFLAGS_BASE="$(CXXFLAGS_BASE) -Werror" $(PROGRAM)
	@echo "CI build successful."

# Test build on all supported configurations
test-configs:
	@echo "Testing release build..."
	@$(MAKE) clean && $(MAKE) release
	@echo "Testing debug build..."
	@$(MAKE) clean && $(MAKE) debug
	@echo "All configurations build successfully."
