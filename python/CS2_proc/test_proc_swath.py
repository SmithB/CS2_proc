#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  1 20:22:30 2024

@author: ben
"""

from CS2_L1 import CS2_L1
import pointCollection as pc
import numpy as np
import matplotlib.pyplot as plt
from proc_swath import proc_swath
L1D_file='../notebooks/CS_LTA__SIR_SIN_1B_20180710T101404_20180710T101529_D001.nc'
D_L1 = CS2_L1().from_nc(L1D_file)
D_sw = proc_swath(D_L1, [1])
