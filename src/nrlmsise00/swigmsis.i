%module c_msis

%{
#define SWIG_FILE_WITH_INIT
#include "nrlmsise-00.h"
%}

%include "nrlmsise-00.h"

%include "carrays.i"
%array_class(double, doubleArray);
%array_class(int, intArray);