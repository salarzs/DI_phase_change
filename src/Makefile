# Makefile

ARCH ?= generic
APP ?= basic
USE_FAST_KERNELS ?= 0

include targets/target.$(ARCH)
include apps/$(APP)/app.$(APP)

POST_PATH = apps/$(APP)/postp.$(APP) -I./postprocessing

.PHONY: clean

.SUFFIXES :
.SUFFIXES : .o .c .f90

.f90.o:
	$(FC) $(FFLAGS) $(DFLAGS) $(CPP_FLAGS) -c -I$(POST_PATH) -o $@ $<

.c.o:
	$(CC) $(CFLAGS) -c -o $@ $<
#
# We first compile the modules that are alway needed
# di_acdi.o added here by salar
OBJS = profiler.o \
bound.o \
di_acdi.o \
chkdiv.o \
chkdt.o \
common_mpi.o \
correc.o \
debug.o \
fft.o \
fftw.o \
fillps.o \
mom.o \
funcs.o \
initflow.o \
initgrid.o \
initmpi.o \
initsolver.o \
load.o \
output.o \
apps/$(APP)/param.o \
rk.o \
source.o \
sanity.o \
2decomp/decomp_2d.o \
2decomp/io.o \
types.o
#
###########################################################
######## TO CUSTOMIZE BASED ON APPLICATION PHYSICS ########
###########################################################
#
ifeq ($(USE_FAST_KERNELS),1) 
	DFLAGS  += -D_FAST_KERNELS_1
	DFLAGS  += -D_FAST_KERNELS_2
	DFLAGS  += -D_FAST_KERNELS_3
	DFLAGS  += -D_FAST_KERNELS_4
endif

ifeq ($(TWOD),1)
	DFLAGS += -D_TWOD
endif

ifeq ($(NAVIER_STOKES),1)
	DFLAGS += -D_NAVIER_STOKES
endif

ifeq ($(TWO_PHASE),1)
	DFLAGS += -D_TWO_PHASE
ifeq ($(USE_VOF),1)
	DFLAGS += -D_USE_VOF
	OBJS += vof.o
endif
ifeq ($(USE_ACDI),1)
	DFLAGS += -D_USE_ACDI
endif
endif
ifeq ($(SCALAR_TRANS),1)    #added by salar
	DFLAGS += -D_SCALAR_TRANS
endif
ifeq ($(TURB_FORCING),1) 
	DFLAGS += -D_TURB_FORCING
endif

ifeq ($(DO_POSTPROC),1)
	OBJS += postprocessing/post.o
ifeq ($(TWO_PHASE),1)
	OBJS += postprocessing/tagging.o
endif
	DFLAGS += -D_DO_POSTPROC
endif

ifeq ($(INIT_MONTECARLO),1) 
	DFLAGS += -D_INIT_MONTECARLO
	OBJS += monte_sphere.o
endif

ifeq ($(CONSTANT_COEFFS_POISSON),1) 
	DFLAGS += -D_CONSTANT_COEFFS_POISSON 
	OBJS += solver_gpu.o solver_cpu.o
else 
	OBJS += solver_vc.o
	LDFLAGS += -lHYPRE
endif

ifeq ($(USE_NVTX),1)  
	DFLAGS += -D_USE_NVTX
	LDFLAGS += $(NVTX_LIB)
endif

ifeq ($(BENCHMARK_NO_IO),1)  
	DFLAGS += -D_BENCHMARK_NO_IO
endif
#

ifeq ($(CONTACT_LINE),1)       ##added by Armin
        DFLAGS += -D_CONTACT_LINE
endif


#
###########################################################
#
TARGET = flutas

all: $(TARGET).$(APP)
	ln -s $(TARGET).$(APP) $(TARGET)

main.o : apps/$(APP)/$(MAIN_SRC) $(OBJS)
	$(FC) $(FFLAGS) $(DFLAGS) $(CPP_FLAGS) -I$(POST_PATH) -c -o main.o $<

$(TARGET).$(APP): $(OBJS) main.o
	$(FC) $(FFLAGS) $(DFLAGS) $(OBJS) $(LDFLAGS) main.o -o $(TARGET).$(APP)

clean-obj:
	rm -rf *.o *.mod *dSYM && rm -rf 2decomp/*.{mod,d,o} && rm -rf postprocessing/*.{mod,d,o} && rm -rf apps/$(APP)/*.{mod,d,o}

clean: clean-obj
	rm -rf $(TARGET)*

include make.deps
