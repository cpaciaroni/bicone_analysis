# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 15:02:16 2017

@author: Caroline
"""
import numpy as np
import scipy.integrate 
import scipy.fftpack
import matplotlib.pyplot as plt
import pandas as pd

def zero_offset(volts):
    #subtract DC offset to keep pulse centered at zero
    avg = np.mean(volts)
    volts0 = volts - avg
    return volts0

def waveform_reader(filename): 
    """Read csv from oscilloscope """
    rdr = pd.read_csv(filename, header=0, names = [0,1,2,3,4],
                      delimiter=',',usecols = [3,4])
                     
    return rdr

def get_data(fname):
    data = waveform_reader(fname)
    
    time = data[3]
    volts = data[4]
    return time, volts

def waveform_plotter(start,stop,step):
	pastel_colors = ['#FF6666','#FFCC66','#CCFF66','#66FF66','#66FFCC','#66FFFF','#66CCFF']
	rev_pastel_colors = pastel_colors[::-1]
	#'#6666FF','#CC66FF','#FF66FF','#FF6FCF']
	an_range = range(start,stop+step,step)
	rev_an_range = an_range[::-1]
	for i in rev_an_range:
		if i == 0:
			fname = '/Users/Caroline/Desktop/Antenna/Bicone_data_0622/CoPol_500averages_00' + str(i) + 'deg_Ch1.csv'
		elif i <= 90:
			fname = '/Users/Caroline/Desktop/Antenna/Bicone_data_0622/CoPol_500averages_0' + str(i) + 'deg_Ch1.csv'
		elif i == 360:
			continue
		else:
			fname = '/Users/Caroline/Desktop/Antenna/Bicone_data_0622/CoPol_500averages_' + str(i) +'deg_Ch1.csv'							
		time,volts_raw = get_data(fname)
		volts = zero_offset(volts_raw)
    
		t = time + np.abs(np.min(time));
		plt.figure(1, figsize=(8,6))
		ax = plt.subplot()
		plt.hold(True)
		plt.plot(t*10**9,volts,c=rev_pastel_colors[(i-90)/15],linewidth=2,label = (str(i) + '$^\circ$'))
		plt.legend(loc="upper right",fontsize=12)
		plt.xlabel('Time (ns)',fontsize=20)
		plt.ylabel('Voltage (V)',fontsize=18)
		plt.title('Voltage on $R_x$ bicone from 90-180$^\circ$',fontsize=22)
		plt.xlim(890,1100)
		plt.ylim(-0.20,0.20)
        
		ax.tick_params(axis='both', which='major', labelsize=14,width=2)
		ax.tick_params(axis='both', which='minor', labelsize=12,width=2)
		ax.spines['top'].set_linewidth(3)
		ax.spines['right'].set_linewidth(3)
		ax.spines['bottom'].set_linewidth(3)
		ax.spines['left'].set_linewidth(3)
        
		plt.savefig('Bicone_IR_pretty_90t180.png',format='png',dpi=500)
	return