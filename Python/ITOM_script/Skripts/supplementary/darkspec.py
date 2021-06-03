# -*- coding: utf-8 -*-
# =============================================================================
# Created By  : Maximilian Sender
# Copyright   : Copyright 2021, Institute of Chemical Engineering, Prof. Dr. Dirk Ziegenbalg, Ulm University'
# License     : GNU LGPL
# =============================================================================
"""The module Darkspec complements the module Radiometric_measurement and creates a dark spectrum for a measurement
with no switchable light source. Must be run in an ITOM environment.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
# =============================================================================
# Imports
# =============================================================================
import numpy as np
import queue
import time
import matplotlib.pyplot as plt

try:
    spec = dataIO("AvantesAvaSpec", 6546, 1639)
except:
    del spec
    spec = dataIO("AvantesAvaSpec", 6546, 1639)
    
int_time = 0.35
averages = 2


def dunkel():
    spec.acquire()
    time.sleep(int_time*averages)
    lambda_table = spec.getParam("lambda_table")
    data = dataObject([1, 2048], 'float32')
    spec.copyVal(data)
    dunkel_spec = (list(data))
    global dunkel_spec
    del dunkel_spec[1062:2048]
    del dunkel_spec[0:203]
    print("dunkel")
    # print ("dunkelend")

spec.startDevice()
spec.setParam("integration_time", int_time)
spec.setParam("dark_correction", 0)
spec.setParam("average", averages)
print("wait")
time.sleep(5)

dunkel()


all_data = dataObject([1,859],'float32', data=dunkel_spec)
#  all_data = dataObject([1,1062],'float32', data = dunkel_spec)
#  plot(all_data, className = 'itom1dqwtplot')
plt.plot(all_data)
savetext = str("TransferCurves/dark_"+str(int_time).replace('.', '_')+'_av'+str(averages)+'.csv')
np.savetxt(savetext, np.transpose(all_data), delimiter=';')

spec.stopDevice()
del spec
