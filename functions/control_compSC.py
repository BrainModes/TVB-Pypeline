#!/bin/python
# -*- coding: utf-8 -*-

import multiprocessing
from computeSC_clusterDK import comp_sc_row

# Create worker pool
pool = multiprocessing.Pool()

# Params
path = '/Users/srothmei/Desktop/charite/toronto/AJ_20140516_1600/mrtrix_68/'
subID = 'AJ_20140516_1600'

# Fire it up
for roi in range(68):
    arguments = path, roi + 1, subID
    print('Iteration ' + str(roi) + ' / 68')
    pool.apply_async(comp_sc_row, arguments)

pool.close()
pool.join()




