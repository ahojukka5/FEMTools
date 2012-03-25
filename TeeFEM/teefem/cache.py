# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 14:20:39 2012

@author: Jukka Aho

Cache

"""

import functools
import teefem

class Cache(object):
    """Decorator that caches a function's return value each time it is called.
    If called later with the same arguments, the cached value is returned, and
    not re-evaluated.
    """

    def __init__(self, func):
        self.func = func
        self.cache = {}
    
    def __call__(self, *args):
        try:
#            if kwds.get('recache', False): 
#                raise KeyError
            return self.cache[args]
        except KeyError:
            value = self.func(*args)
            self.cache[args] = value
            return value
        except TypeError:
            # uncachable -- for instance, passing a list as an argument.
            # Better to not cache than to blow up entirely.
            return self.func(*args)
    
    def __repr__(self):
        """Return the function's docstring."""
        return self.func.__doc__
        
    def __get__(self, obj, objtype):
        """Support instance methods."""
        return functools.partial(self.__call__, obj)

cache = Cache

@cache
def fibonacci(n):
   "Return the nth fibonacci number."
   if n in (0, 1):
      return n
   return fibonacci(n-1) + fibonacci(n-2)

if __name__ == '__main__':
    import time
    t = 0
    for i in range(5):
        t0 = time.clock()
        fib = fibonacci(30)
        t1 = time.clock()
        t += t1-t0
        print("{0}. That took {1} seconds".format(fib, t1-t0))
    print("Total {0} seconds".format(t))

teefem.log.info("Module {0} loaded.".format(__file__))