# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 15:51:02 2017

@author: Caroline
"""

import glob
import csv
import re
import sys
import pandas as pd
sys.path.append("../ana")
sys.path.append("../filter/")
sys.path.append("../sim_work/")
import dsphelper
import filters as filter
import numpy as np
import analysis_library as al
import matplotlib.pyplot as pyp
import endpoints
import numpy.fft as fft
import scipy.signal

z = 319
head, ec = endpoints.getEndpointWave(z, dirc='/u/project/saltzber/swissel/slac/analysis/slac/t510/endpoint/Mag_2400A_Run2_Pol_01ns_2/')

times2 = []
eh2 = []
ev2 = []
ex2 = []
ey2 = []
ez2 = []
print ec
# simple downsample by a factor of 2
# just skip every other one
# not as good as using the scipy.signal.resample function
# so just using this to grab the right time values
for i in range(len(ec.wave)):
	if( (i) % 2 == 0):
		times2.append(ec.wave.time.values[i])
'''		eh2.append(ec.wave.Eh.values[i])
		ev2.append(ec.wave.Ev.values[i])
		ex2.append(ec.wave.Ex.values[i])
		ey2.append(ec.wave.Ey.values[i])
		ez2.append(ec.wave.Ez.values[i])
'''

eh2 = scipy.signal.resample(ec.wave.Eh.values, len(ec.wave.Eh.values)/2)
ev2 = scipy.signal.resample(ec.wave.Ev.values, len(ec.wave.Ev.values)/2)
#times2 = scipy.signal.resample(ec.wave.time.values, len(ec.wave.time.values)/2)

#freq = fft.fftfreq(n=9999, d=0.2e-9)
ffth = fft.fft(eh2, n=9999)
fftv = fft.fft(ev2, n=9999)
freq = fft.fftfreq(n=9999, d=0.1e-9*2)
#ffth = fft.rfft(eh2, n=10000)
#fftv = fft.rfft(ev2, n=10000)
#freq = fft.fftfreq(n=5001, d=0.1e-9)
print freq[0:10]
print ffth[0:10]
print fftv[0:10]
print len(freq), len(ffth), len(fftv)

pyp.figure(1)
dt = ec.wave.time.values[1]-ec.wave.time.values[0]
dt2 = times2[1]-times2[0]
pyp.semilogy(ec.fd.freq.values, abs(ec.fd.Ev)*dt, 'b', linewidth=2, alpha=0.5)
pyp.semilogy(ec.fd.freq.values, abs(ec.fd.Eh)*dt, 'r', linewidth=2, alpha=0.5)
pyp.semilogy(fft.fftshift(freq), abs(fft.fftshift(ffth))*dt2, 'r--', linewidth=4, alpha=0.75)
pyp.semilogy(fft.fftshift(freq), abs(fft.fftshift(fftv))*dt2, 'b--', linewidth=4, alpha=0.75)
pyp.xlabel("Frequency (Hz)")
pyp.ylabel("|FFT|")
#pyp.ylim(10, 2000)
#pyp.xlim(0, 2.5e9)

pyp.figure(2)
pyp.plot(ec.wave.time.values, ec.wave.Eh, 'b', linewidth=2, alpha=0.5)
pyp.plot(ec.wave.time.values, ec.wave.Ev, 'r', linewidth=2, alpha=0.5)
pyp.plot(times2, eh2, 'b--', linewidth=4, alpha=0.75)
pyp.plot(times2, ev2, 'r--', linewidth=2, alpha=0.75)
pyp.xlim(4.40e-8, 4.70e-8)
pyp.show()

df = pd.DataFrame({'freq':freq, 'fft_Eh_real':ffth.real, 'fft_Eh_imag':ffth.imag, 'fft_Ev_real':fftv.real, 'fft_Ev_imag':fftv.imag}, index=freq)
df.to_csv("Mag_2400A_Run2_Pol_01ns_2_%d_downsampled_5Gsps.csv"%z)

#### this reads the file as a pandas dataframe
df = pd.read_csv("Mag_2400A_Run2_Pol_01ns_2_%d_downsampled_5Gsps.csv"%z)
# this part makes the complex numbers
print df.head()
df['fft_Eh'] = df['fft_Eh_real'] + 1j * df['fft_Eh_imag']
df['fft_Ev'] = df['fft_Ev_real'] + 1j * df['fft_Ev_imag']
print df.head()
#### this gives you pandas series
print df.freq
print df.fft_Eh
print df.fft_Ev
### this is to get numpy arrays
print df.freq.values
print df.fft_Eh.values
print df.fft_Ev.values
