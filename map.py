# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 15:36:47 2021

@author: Nikit
"""
import numpy as np
import time
from numba import njit

list1 = np.linspace(1,10, 5)
list2 = np.linspace(1,60, 5)
list3 = np.linspace(1,120, 5)
list4 = np.array(list(map(lambda x, y: x+y, list1, list2)))

@njit
def fun(x, y):
    return x + y

start = time.time()
for i in range(10000000):
    fun(10, 5)
print(time.time() - start)