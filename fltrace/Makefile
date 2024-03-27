# Set compiler, flags and netcdf library according to machine we are using:

MPIF90 ?= /usr/lib64/openmpi/bin/mpif90
FFLAGS = -O3 -fcheck=all -Wuninitialized -march=native -fimplicit-none -Wall -Wextra -ffast-math -funroll-loops --param max-unroll-times=5
NETCDF = -I /usr/lib64/gfortran/modules
NETCDFLIB = -L/usr/lib64/libnetcdff.so.7 -lnetcdff

MODULEFLAG = -module $(OBJDIR)
TARGET = fltrace

# --------------------------------------------------
# Shouldn't need to touch below here
# --------------------------------------------------


all: main

SDF := SDF/FORTRAN
SDFMOD = $(SDF)/include/sdf.mod
SRCDIR = src
OBJDIR = obj
BINDIR = bin
FC = $(MPIF90)
DATE := $(shell date +%s)
MACHINE := $(shell uname -n)
PREPROFLAGS = $(DEFINES) $(D)_COMMIT='"$(COMMIT)"' $(D)_DATE=$(DATE) \
  $(D)_MACHINE='"$(MACHINE)"'



SRCFILES = shared_data.o gridpoint_tools.o grid.o import_export.o fltrace.o

OBJFILES := $(SRCFILES:.f90=.o)
OBJFILES := $(OBJFILES:.F90=.o)

FULLTARGET = $(BINDIR)/$(TARGET)

VPATH = $(SRCDIR):$(SRCDIR)/core:$(SDF)/src:$(OBJDIR)

-include $(SRCDIR)/COMMIT

ifeq ($(DONE_COMMIT),)
main: commit
else
main: $(FULLTARGET)
endif

# Rule to build the fortran files

%.o: %.f90
	@mkdir -p $(BINDIR) $(OBJDIR)
	$(FC) -c $(FFLAGS) $(NETCDF) -J $(OBJDIR) -o $(OBJDIR)/$@ $<

%.o: %.F90
	@mkdir -p $(BINDIR) $(OBJDIR)
	$(FC) -c $(FFLAGS) $(NETCDF) -J $(OBJDIR) -o $(OBJDIR)/$@ $(PREPROFLAGS) $<

$(FULLTARGET): $(OBJFILES)
	$(FC) $(FFLAGS) -J $(OBJDIR) -o $@ $(addprefix $(OBJDIR)/,$(OBJFILES)) $(NETCDFLIB)


clean:
	@rm -rf $(BINDIR) $(OBJDIR)


#$(OBJFILES): | $(OBJDIR)

$(OBJDIR):
	@mkdir -p $(OBJDIR)

commit: FORCE
	@sh $(SRCDIR)/gen_commit_string.sh && $(MAKE) $(MAKECMDGOALS) DONE_COMMIT=1

FORCE:

.PHONY: commit clean cleanall tidy datatidy visit visitclean main FORCE

shared_data.o: shared_data.f90
gridpoint_tools.o: gridpoint_tools.f90
grid.o: grid.f90 shared_data.o gridpoint_tools.o
import_export.o: import_export.f90 shared_data.o
fltrace.o: fltrace.f90 shared_data.o grid.o import_export.o


