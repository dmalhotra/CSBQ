SCTL_INCLUDE_DIR ?= ./SCTL/include
SCTL_DATA_PATH ?= ./data

CXX=g++ # requires g++-8 or newer / icpc (with gcc compatibility 7.5 or newer) / clang++ with llvm-10 or newer
CXXFLAGS = -std=c++11 -fopenmp -Wall -Wfloat-conversion # need C++11 and OpenMP

CXXFLAGS += -DSCTL_DATA_PATH=$(SCTL_DATA_PATH) # path for quadrature tables
CXXFLAGS += -DSCTL_QUAD_T=__float128 # enable quadruple precision (for stored quadrature tables)

#Optional flags
DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CXXFLAGS += -O0 -fsanitize=address,leak,undefined,pointer-compare,pointer-subtract,float-divide-by-zero,float-cast-overflow -fno-sanitize-recover=all -fstack-protector # debug build
  CXXFLAGS += -DSCTL_MEMDEBUG # Enable memory checks
else
	CXXFLAGS += -O3 -march=native -DNDEBUG # release build
endif

ifeq ($(shell uname -s),Darwin)
	CXXFLAGS += -g -rdynamic -Wl,-no_pie # for stack trace (on Mac)
else
	CXXFLAGS += -g -rdynamic # for stack trace
endif

CXXFLAGS += -DSCTL_PROFILE=5 -DSCTL_VERBOSE # enable profiling

#CXXFLAGS += -DSCTL_HAVE_MPI # enable MPI

CXXFLAGS += -lblas -DSCTL_HAVE_BLAS # use BLAS
CXXFLAGS += -llapack -DSCTL_HAVE_LAPACK # use LAPACK
#CXXFLAGS += -lopenblas -DSCTL_HAVE_BLAS -DSCTL_HAVE_LAPACK # use BLAS and LAPACK (OpenBLAS)
#CXXFLAGS += -mkl -DSCTL_HAVE_BLAS -DSCTL_HAVE_LAPACK # use MKL BLAS and LAPACK (Intel compiler)
#CXXFLAGS += -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -DSCTL_HAVE_BLAS -DSCTL_HAVE_LAPACK # use MKL BLAS and LAPACK (non-Intel compiler)
#CXXFLAGS += -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm -ldl -DSCTL_HAVE_BLAS -DSCTL_HAVE_LAPACK
#CXXFLAGS += -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl -DSCTL_HAVE_BLAS -DSCTL_HAVE_LAPACK

CXXFLAGS += -lfftw3 -DSCTL_HAVE_FFTW
CXXFLAGS += -lfftw3f -DSCTL_HAVE_FFTWF
CXXFLAGS += -lfftw3l -DSCTL_HAVE_FFTWL

#CXXFLAGS += -DSCTL_HAVE_SVML

#PVFMM_INC_DIR = $(PVFMM_DIR)/include
#PVFMM_LIB_DIR = $(PVFMM_DIR)/lib/.libs
#CXXFLAGS += -DSCTL_HAVE_PVFMM -I$(PVFMM_INC_DIR)
#LDLIBS += $(PVFMM_LIB_DIR)/libpvfmm.a

RM = rm -f
MKDIRS = mkdir -p

BINDIR = ./bin
OBJDIR = ./obj
INCDIR = ./include
TESTDIR = ./test
TUTORDIR = ./tutorial

TARGET_BIN = \
       $(BINDIR)/demo1-geometry \
       $(BINDIR)/demo2-quadrature

TEST_BIN = \
       $(BINDIR)/bvp-solve \
       $(BINDIR)/centerline-force \
       $(BINDIR)/greens-identity \
       $(BINDIR)/mobility \
       $(BINDIR)/spectrum \
       $(BINDIR)/stokes-drag \
       $(BINDIR)/stokes-mutual-drag \
       $(BINDIR)/tangle-adap-geom

all : $(TARGET_BIN)

test : $(TEST_BIN)

$(BINDIR)/%: $(OBJDIR)/%.o
	-@$(MKDIRS) $(dir $@)
	$(CXX) $^ $(LDLIBS) -o $@ $(CXXFLAGS)

$(OBJDIR)/%.o: $(TUTORDIR)/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -I$(SCTL_INCLUDE_DIR) -c $^ -o $@

$(OBJDIR)/%.o: $(TESTDIR)/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -I$(SCTL_INCLUDE_DIR) -c $^ -o $@

clean:
	$(RM) -r $(BINDIR)/* $(OBJDIR)/*
	$(RM) *~ */*~ */*/*~
