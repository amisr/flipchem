# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.12
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info >= (2, 7, 0):
    def swig_import_helper():
        import importlib
        pkg = __name__.rpartition('.')[0]
        mname = '.'.join((pkg, '_c_msis')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_c_msis')
    _c_msis = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_c_msis', [dirname(__file__)])
        except ImportError:
            import _c_msis
            return _c_msis
        try:
            _mod = imp.load_module('_c_msis', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _c_msis = swig_import_helper()
    del swig_import_helper
else:
    import _c_msis
del _swig_python_version_info

try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_setattr_nondynamic(self, class_type, name, value, static=1):
    if (name == "thisown"):
        return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name, None)
    if method:
        return method(self, value)
    if (not static):
        if _newclass:
            object.__setattr__(self, name, value)
        else:
            self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)


def _swig_setattr(self, class_type, name, value):
    return _swig_setattr_nondynamic(self, class_type, name, value, 0)


def _swig_getattr(self, class_type, name):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    raise AttributeError("'%s' object has no attribute '%s'" % (class_type.__name__, name))


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except __builtin__.Exception:
    class _object:
        pass
    _newclass = 0

class doubleArray(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, doubleArray, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, doubleArray, name)
    __repr__ = _swig_repr

    def __init__(self, nelements):
        this = _c_msis.new_doubleArray(nelements)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _c_msis.delete_doubleArray
    __del__ = lambda self: None

    def __getitem__(self, index):
        return _c_msis.doubleArray___getitem__(self, index)

    def __setitem__(self, index, value):
        return _c_msis.doubleArray___setitem__(self, index, value)

    def cast(self):
        return _c_msis.doubleArray_cast(self)
    if _newclass:
        frompointer = staticmethod(_c_msis.doubleArray_frompointer)
    else:
        frompointer = _c_msis.doubleArray_frompointer
doubleArray_swigregister = _c_msis.doubleArray_swigregister
doubleArray_swigregister(doubleArray)

def doubleArray_frompointer(t):
    return _c_msis.doubleArray_frompointer(t)
doubleArray_frompointer = _c_msis.doubleArray_frompointer

class intArray(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, intArray, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, intArray, name)
    __repr__ = _swig_repr

    def __init__(self, nelements):
        this = _c_msis.new_intArray(nelements)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _c_msis.delete_intArray
    __del__ = lambda self: None

    def __getitem__(self, index):
        return _c_msis.intArray___getitem__(self, index)

    def __setitem__(self, index, value):
        return _c_msis.intArray___setitem__(self, index, value)

    def cast(self):
        return _c_msis.intArray_cast(self)
    if _newclass:
        frompointer = staticmethod(_c_msis.intArray_frompointer)
    else:
        frompointer = _c_msis.intArray_frompointer
intArray_swigregister = _c_msis.intArray_swigregister
intArray_swigregister(intArray)

def intArray_frompointer(t):
    return _c_msis.intArray_frompointer(t)
intArray_frompointer = _c_msis.intArray_frompointer


def PyFloat_AsDouble(arg1):
    return _c_msis.PyFloat_AsDouble(arg1)
PyFloat_AsDouble = _c_msis.PyFloat_AsDouble
class nrlmsise_flags(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, nrlmsise_flags, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, nrlmsise_flags, name)
    __repr__ = _swig_repr
    __swig_setmethods__["switches"] = _c_msis.nrlmsise_flags_switches_set
    __swig_getmethods__["switches"] = _c_msis.nrlmsise_flags_switches_get
    if _newclass:
        switches = _swig_property(_c_msis.nrlmsise_flags_switches_get, _c_msis.nrlmsise_flags_switches_set)
    __swig_setmethods__["sw"] = _c_msis.nrlmsise_flags_sw_set
    __swig_getmethods__["sw"] = _c_msis.nrlmsise_flags_sw_get
    if _newclass:
        sw = _swig_property(_c_msis.nrlmsise_flags_sw_get, _c_msis.nrlmsise_flags_sw_set)
    __swig_setmethods__["swc"] = _c_msis.nrlmsise_flags_swc_set
    __swig_getmethods__["swc"] = _c_msis.nrlmsise_flags_swc_get
    if _newclass:
        swc = _swig_property(_c_msis.nrlmsise_flags_swc_get, _c_msis.nrlmsise_flags_swc_set)

    def __init__(self):
        this = _c_msis.new_nrlmsise_flags()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _c_msis.delete_nrlmsise_flags
    __del__ = lambda self: None
nrlmsise_flags_swigregister = _c_msis.nrlmsise_flags_swigregister
nrlmsise_flags_swigregister(nrlmsise_flags)

class ap_array(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ap_array, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ap_array, name)
    __repr__ = _swig_repr
    __swig_setmethods__["a"] = _c_msis.ap_array_a_set
    __swig_getmethods__["a"] = _c_msis.ap_array_a_get
    if _newclass:
        a = _swig_property(_c_msis.ap_array_a_get, _c_msis.ap_array_a_set)

    def __init__(self):
        this = _c_msis.new_ap_array()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _c_msis.delete_ap_array
    __del__ = lambda self: None
ap_array_swigregister = _c_msis.ap_array_swigregister
ap_array_swigregister(ap_array)

class nrlmsise_input(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, nrlmsise_input, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, nrlmsise_input, name)
    __repr__ = _swig_repr
    __swig_setmethods__["year"] = _c_msis.nrlmsise_input_year_set
    __swig_getmethods__["year"] = _c_msis.nrlmsise_input_year_get
    if _newclass:
        year = _swig_property(_c_msis.nrlmsise_input_year_get, _c_msis.nrlmsise_input_year_set)
    __swig_setmethods__["doy"] = _c_msis.nrlmsise_input_doy_set
    __swig_getmethods__["doy"] = _c_msis.nrlmsise_input_doy_get
    if _newclass:
        doy = _swig_property(_c_msis.nrlmsise_input_doy_get, _c_msis.nrlmsise_input_doy_set)
    __swig_setmethods__["sec"] = _c_msis.nrlmsise_input_sec_set
    __swig_getmethods__["sec"] = _c_msis.nrlmsise_input_sec_get
    if _newclass:
        sec = _swig_property(_c_msis.nrlmsise_input_sec_get, _c_msis.nrlmsise_input_sec_set)
    __swig_setmethods__["alt"] = _c_msis.nrlmsise_input_alt_set
    __swig_getmethods__["alt"] = _c_msis.nrlmsise_input_alt_get
    if _newclass:
        alt = _swig_property(_c_msis.nrlmsise_input_alt_get, _c_msis.nrlmsise_input_alt_set)
    __swig_setmethods__["g_lat"] = _c_msis.nrlmsise_input_g_lat_set
    __swig_getmethods__["g_lat"] = _c_msis.nrlmsise_input_g_lat_get
    if _newclass:
        g_lat = _swig_property(_c_msis.nrlmsise_input_g_lat_get, _c_msis.nrlmsise_input_g_lat_set)
    __swig_setmethods__["g_long"] = _c_msis.nrlmsise_input_g_long_set
    __swig_getmethods__["g_long"] = _c_msis.nrlmsise_input_g_long_get
    if _newclass:
        g_long = _swig_property(_c_msis.nrlmsise_input_g_long_get, _c_msis.nrlmsise_input_g_long_set)
    __swig_setmethods__["lst"] = _c_msis.nrlmsise_input_lst_set
    __swig_getmethods__["lst"] = _c_msis.nrlmsise_input_lst_get
    if _newclass:
        lst = _swig_property(_c_msis.nrlmsise_input_lst_get, _c_msis.nrlmsise_input_lst_set)
    __swig_setmethods__["f107A"] = _c_msis.nrlmsise_input_f107A_set
    __swig_getmethods__["f107A"] = _c_msis.nrlmsise_input_f107A_get
    if _newclass:
        f107A = _swig_property(_c_msis.nrlmsise_input_f107A_get, _c_msis.nrlmsise_input_f107A_set)
    __swig_setmethods__["f107"] = _c_msis.nrlmsise_input_f107_set
    __swig_getmethods__["f107"] = _c_msis.nrlmsise_input_f107_get
    if _newclass:
        f107 = _swig_property(_c_msis.nrlmsise_input_f107_get, _c_msis.nrlmsise_input_f107_set)
    __swig_setmethods__["ap"] = _c_msis.nrlmsise_input_ap_set
    __swig_getmethods__["ap"] = _c_msis.nrlmsise_input_ap_get
    if _newclass:
        ap = _swig_property(_c_msis.nrlmsise_input_ap_get, _c_msis.nrlmsise_input_ap_set)
    __swig_setmethods__["ap_a"] = _c_msis.nrlmsise_input_ap_a_set
    __swig_getmethods__["ap_a"] = _c_msis.nrlmsise_input_ap_a_get
    if _newclass:
        ap_a = _swig_property(_c_msis.nrlmsise_input_ap_a_get, _c_msis.nrlmsise_input_ap_a_set)

    def __init__(self):
        this = _c_msis.new_nrlmsise_input()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _c_msis.delete_nrlmsise_input
    __del__ = lambda self: None
nrlmsise_input_swigregister = _c_msis.nrlmsise_input_swigregister
nrlmsise_input_swigregister(nrlmsise_input)

class nrlmsise_output(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, nrlmsise_output, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, nrlmsise_output, name)
    __repr__ = _swig_repr
    __swig_setmethods__["d"] = _c_msis.nrlmsise_output_d_set
    __swig_getmethods__["d"] = _c_msis.nrlmsise_output_d_get
    if _newclass:
        d = _swig_property(_c_msis.nrlmsise_output_d_get, _c_msis.nrlmsise_output_d_set)
    __swig_setmethods__["t"] = _c_msis.nrlmsise_output_t_set
    __swig_getmethods__["t"] = _c_msis.nrlmsise_output_t_get
    if _newclass:
        t = _swig_property(_c_msis.nrlmsise_output_t_get, _c_msis.nrlmsise_output_t_set)

    def __init__(self):
        this = _c_msis.new_nrlmsise_output()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _c_msis.delete_nrlmsise_output
    __del__ = lambda self: None
nrlmsise_output_swigregister = _c_msis.nrlmsise_output_swigregister
nrlmsise_output_swigregister(nrlmsise_output)


def gtd7(input, flags, output):
    return _c_msis.gtd7(input, flags, output)
gtd7 = _c_msis.gtd7

def gtd7d(input, flags, output):
    return _c_msis.gtd7d(input, flags, output)
gtd7d = _c_msis.gtd7d

def gts7(input, flags, output):
    return _c_msis.gts7(input, flags, output)
gts7 = _c_msis.gts7

def ghp7(input, flags, output, press):
    return _c_msis.ghp7(input, flags, output, press)
ghp7 = _c_msis.ghp7
# This file is compatible with both classic and new-style classes.


