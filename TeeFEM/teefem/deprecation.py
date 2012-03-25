# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 07:47:37 2012

DeprecationWarning decorator

"""


import warnings
import functools

def deprecated(func):
    """This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used."""
    
    @functools.wraps(func)
    def new_func(*args, **kwargs):
        warnings.warn_explicit(
            "Call to deprecated function %(funcname)s." % {
                'funcname': func.__name__,
            },
            category=DeprecationWarning,
            filename=func.func_code.co_filename,
            lineno=func.func_code.co_firstlineno + 1
        )
        return func(*args, **kwargs)
    return new_func


if __name__ == '__main__':
    
    ## Usage examples ##
    @deprecated
    def my_func():
        print("my_func()")
    
    my_func()