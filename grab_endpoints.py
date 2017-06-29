# -*- coding: utf-8 -*-
"""
Created on Fri May 12 13:25:32 2017

Reads endpoint files from anne & returns horizontal and vertical components 
of field

@author: Caroline
"""

import pandas as pd
import matplotlib.pyplot as plt
import endpoints
import numpy.fft as fft
import scipy.signal


# folder for 0 B field: Mag_0A_Caroline/
#folder for Bfield on: VHF_Sims/

def GetEndpoints(z,files=True, plots=False):
    head, ec = endpoints.getEndpointWave(z, dirc='/Users/Caroline/Desktop/bicones/New_Sims/Mag_0A_Caroline/')

    times2 = []
    eh2 = []
    ev2 = []
    ex2 = []
    ey2 = []
    ez2 = []
    

    # simple downsample by a factor of 2
    # just skip every other one
    # not as good as using the scipy.signal.resample function
    # so just using this to grab the right time values
    
    #
    #for i in range(len(Endpoints.index)):
    #	if( (i) % 2 == 0):
    #		times2.append(Endpoints.time[i])
    #		
    #'''		eh2.append(ec.wave.Eh.values[i])
    #		ev2.append(ec.wave.Ev.values[i])
    #		ex2.append(ec.wave.Ex.values[i])
    #		ey2.append(ec.wave.Ey.values[i])
    #		ez2.append(ec.wave.Ez.values[i])
    #'''
    
    eh2 = scipy.signal.resample(ec.wave.Eh.values, len(ec.wave.Eh.values)/2)
    ev2 = scipy.signal.resample(ec.wave.Ev.values, len(ec.wave.Ev.values)/2)
    times2 = scipy.signal.resample(ec.wave.time.values, len(ec.wave.time.values)/2)
    
    #freq = fft.fftfreq(n=9999, d=0.2e-9)
    ffth = fft.fft(eh2, n=9999)
    fftv = fft.fft(ev2, n=9999)
    freq = fft.fftfreq(n=9999, d=0.1e-9*2)
    #ffth = fft.rfft(eh2, n=10000)
    #fftv = fft.rfft(ev2, n=10000)
    #freq = fft.fftfreq(n=5001, d=0.1e-9)
    
    
    #print freq[0:10]
    #print ffth[0:10]
    #print fftv[0:10]
    #print len(freq), len(ffth), len(fftv)

    if plots:
        plt.figure(1)
        dt = ec.wave.time.values[1]-ec.wave.time.values[0]
        dt2 = times2[1]-times2[0]
        plt.semilogy(ec.fd.freq.values, abs(ec.fd.Ex)*dt, 'r', linewidth=2, alpha=0.5)
        plt.semilogy(ec.fd.freq.values, abs(ec.fd.Ey)*dt, 'b', linewidth=2, alpha=0.5)
        plt.semilogy(ec.fd.freq.values, abs(ec.fd.Ez)*dt, 'g', linewidth=2, alpha=0.5)
        #plt.semilogy(fft.fftshift(freq), abs(fft.fftshift(ffth))*dt2, 'b--', linewidth=4, alpha=0.75)
        #plt.semilogy(fft.fftshift(freq), abs(fft.fftshift(fftv))*dt2, 'r--', linewidth=4, alpha=0.75)
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("|FFT|")
        #plt.ylim(10, 2000)
        #plt.xlim(0, 2.5e9) 
        
        plt.figure(2)
        ax = plt.subplot()
        plt.plot(ec.wave.time.values, ec.wave.Ex, 'r', linewidth=2, alpha=0.5,label='$E_x$')
        plt.plot(ec.wave.time.values, ec.wave.Ey, 'b', linewidth=2, alpha=0.5,label='$E_y$')
        plt.plot(ec.wave.time.values, ec.wave.Ez, 'g', linewidth=2, alpha=0.5,label='$E_z$')
        #plt.plot(times2, eh2, 'b--', linewidth=4, alpha=0.75) #downsampled
        #plt.plot(times2, ev2, 'r--', linewidth=2, alpha=0.75)
        plt.title('E-field components with B-field on (I = 2400A)',fontsize=14)
        plt.xlabel("Time (s)",fontsize = 14)
        plt.ylabel("Electric Field (V/m)",fontsize = 14)
        plt.legend(loc="lower left",fontsize=14)
        ax.tick_params(axis='both', which='major', labelsize=14,width=2)
        ax.tick_params(axis='both', which='minor', labelsize=12,width=2)
        plt.xlim(4.8e-8, 5.25e-8)
        ax.spines['top'].set_linewidth(1.5)
        ax.spines['right'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_linewidth(1.5)
        
        #plt.savefig('/Users/Caroline/Desktop/Bicones/sim_plots/E_comps_B_off.png',format='png',dpi=500)
                
        
        plt.figure(3)
        ax = plt.subplot()
        plt.plot(ec.wave.time.values, ec.wave.Eh, 'r', linewidth=2, alpha=0.5, label = '$E_H$')
        plt.plot(ec.wave.time.values, ec.wave.Ev, 'b', linewidth=2, alpha=0.5, label = '$E_V$')
        #plt.plot(times2, eh2, 'b--', linewidth=4, alpha=0.75)
        #plt.plot(times2, ev2, 'r--', linewidth=2, alpha=0.75)
        plt.title('Simulated E-fields with B-field off (I = 0A)',fontsize=14)
        plt.xlabel("Time (s)",fontsize = 14)
        plt.ylabel("Electric Field (V/m)",fontsize = 14)
        plt.legend(loc="upper left",fontsize=14)
        ax.tick_params(axis='both', which='major', labelsize=14,width=2)
        ax.tick_params(axis='both', which='minor', labelsize=12,width=2)
        plt.xlim(4.8e-8, 5.25e-8)
        #plt.ylim(0,35)
        ax.spines['top'].set_linewidth(1.5)
        ax.spines['right'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_linewidth(1.5)
        
        #plt.savefig('/Users/Caroline/Desktop/Bicones/E_Sims_B_off.png',format='png',dpi=500)
        
    if files:
        
        df = pd.DataFrame({'freq':freq, 'fft_Eh_real':ffth.real, 'fft_Eh_imag':ffth.imag, 'fft_Ev_real':fftv.real, 'fft_Ev_imag':fftv.imag}, index=freq)
        
        filename = "/Users/Caroline/Desktop/bicones/Downsampled_VHFsims/Mag_2400A_%d_downsampled_5Gsps.csv"%z
        df.to_csv(filename)
        
        '''#### this reads the file as a pandas dataframe
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
        print df.fft_Ev.values '''
        
        print filename
        
        return filename
    

