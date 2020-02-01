# Make program
MAKE =  make

# GNU Fortran compiler
FC =  gfortran

# GNU C compiler
CC =  gcc

# Library directory
LIBDIR = ""

# ----------------------------------------------------------------------
#       Library directives
# ----------------------------------------------------------------------

LDLIBS = -lm

# ----------------------------------------------------------------------
#       Fortran compiler options
# ----------------------------------------------------------------------

#auto detect if mac or linux
UNAME := $(shell uname -s)

ifeq ($(UNAME), Linux)
	FFLAGS = -O -fPIC -Wall -Wuninitialized -fno-second-underscore --std=legacy
endif
ifeq ($(UNAME), Darwin)
	FFLAGS = -D_GNU_SOURCE -D_DARWIN -Wuninitialized -fno-second-underscore
endif

# ----------------------------------------------------------------------
#       Link editor options
# ----------------------------------------------------------------------

LDFLAGS = -Xlinker

# ----------------------------------------------------------------------
#       Fortran source files
# ----------------------------------------------------------------------

SOURCES = FLIP-CHEM.f KEMPRN.f PESIMP.f Rates.f RSPRIM.f MSIS00.f

# ----------------------------------------------------------------------
#       Object files
# ----------------------------------------------------------------------
OBJECTS = $(SOURCES:%.f=%.o)

# ----------------------------------------------------------------------
#       Libraries and Programs
# ----------------------------------------------------------------------

all: clean libflip.a libflip.so

libflip.a: $(OBJECTS)
	ar -cvq $@ $(OBJECTS)

libflip.so: $(OBJECTS)
	$(FC) -o $@ -shared $(OBJECTS)

clean:
	@echo "Cleaning ... "
	rm -f *~ *.o *.a *.so
