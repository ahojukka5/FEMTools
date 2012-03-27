# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 17:34:59 2012

@author: Jukka Aho

Työkalu matriisien manipulointiin

"""

from __future__ import division
import numpy as np

def removedofs(A, doflist):
    ''' 
    Poistaa vapausasteita listan mukaisesti matriisista / taulukosta / vektorista 
    
    >>> U = np.array(range(25))
    >>> K = U.reshape(5,5)
    >>> U
    array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
           17, 18, 19, 20, 21, 22, 23, 24])
    >>> K
    array([[ 0,  1,  2,  3,  4],
           [ 5,  6,  7,  8,  9],
           [10, 11, 12, 13, 14],
           [15, 16, 17, 18, 19],
           [20, 21, 22, 23, 24]])
    >>> removedofs(U, [4,1,8])
    array([ 0,  2,  3,  5,  6,  7,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
           20, 21, 22, 23, 24])
    >>> removedofs(K, [3,1])
    array([[ 0,  2,  4],
           [10, 12, 14],
           [20, 22, 24]])
           
    '''
    # np.compress
    dim = len(A.shape)
    if dim == 1:
        return np.delete(A,np.s_[doflist])
    else:
        return np.delete(np.delete(A,np.s_[doflist],axis=0), np.s_[doflist], axis=1)

def insertdofs(U, doflist):
    '''
    Lisää vapausasteita (nollia) listan mukaisesti vektoriin
    >>> U = np.array([ 0,  2,  3,  5,  6,  7,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24])
    >>> doflist = np.sort(np.array([4,1,8]))
    >>> insertdofs(U, doflist)
    array([ 0,  0,  2,  3,  0,  5,  6,  7,  0,  9, 10, 11, 12, 13, 14, 15, 16,
           17, 18, 19, 20, 21, 22, 23, 24])
    
    '''
    # np.compress
    dim = len(U.shape)
    ilist = list(np.zeros(len(doflist)))
    if dim == 1:
        return np.insert(U,np.sort(doflist)-range(len(doflist)),ilist)
    else:
        raise ValueError

def test1():
    doflist = np.array([1,3]) # 1 -> poistaa 2. ja 4. rivin
    doflist = np.sort(doflist)
    K = np.random.rand(5,5)
    U = np.random.rand(5)
    print K
    print U
    newK = removedofs(K, doflist)
    newU = removedofs(U, doflist)
    print newK
    print newU

#    doflist -= range(len(doflist)) # Uudet vapausasteet
    new2U = insertdofs(newU, doflist)
#    print new2U
    

if __name__ == '__main__':
    #test1()
    import doctest
    doctest.testmod()