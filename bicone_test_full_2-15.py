# -*- coding: utf-8 -*-
"""
2/9/17 22:46:05 2017

@author: Caroline
"""
import numpy as np
import scipy.integrate 
import scipy.fftpack
import scipy.constants
import matplotlib.pyplot as plt
import pandas as pd
import seaborn.apionly as sns
import cmath
import sys
sys.path.append('/Users/Caroline/A3ImpulseResponse')
import tfUtils as tf


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
    
def get_power(volts,times):
    V = [v**2 for v in volts]
    integral = scipy.integrate.simps(V, times,even='last')
    power = integral/50
    return power
        
 
def waveform_plotter(start,stop,step):
    for i in range(start,stop+step,step):
        if i == 0:
		fname = '/Users/Caroline/Desktop/Antenna/Bicone_data_0622/CoPol_500averages_00' + str(i) + 'deg_Ch1.csv'
		col = '#29A2E8'							
        elif i <= 90:
		fname = '/Users/Caroline/Desktop/Antenna/Bicone_data_0622/CoPol_500averages_0' + str(i) + 'deg_Ch1.csv'
		col = '#9129E8'								
        elif i == 360:
		continue
        else:
		fname = '/Users/Caroline/Desktop/Antenna/Bicone_data_0622/CoPol_500averages_' + str(i) +'deg_Ch1.csv'
		col = '#E83229'								
        time,volts_raw = get_data(fname)
        volts = zero_offset(volts_raw)
    
        t = time + np.abs(np.min(time));
        
        
        plt.figure(1, figsize=(8,6))
        ax = plt.subplot()
        plt.hold(True)
        plt.plot(t*10**9,volts,c=col,linewidth=0.5,label = (str(i) + '$^\circ$'))
        plt.legend(loc="upper right",fontsize=14)
        plt.xlabel('Time (ns)',fontsize=20)
        plt.ylabel('Voltage (V)',fontsize=18)
        plt.title('Voltage on $R_x$ bicone from 0-360$^\circ$',fontsize=22)
        plt.xlim(890,1100)
        plt.ylim(-0.20,0.20)
        
        ax.tick_params(axis='both', which='major', labelsize=14,width=2)
        ax.tick_params(axis='both', which='minor', labelsize=12,width=2)
        ax.spines['top'].set_linewidth(3)
        ax.spines['right'].set_linewidth(3)
        ax.spines['bottom'].set_linewidth(3)
        ax.spines['left'].set_linewidth(3)
        
        plt.savefig('Bicone_impulse_response1.png',format='png',dpi=500)
        
    return 
    
def FFT(time,volts):
    #make the time array start at 0 - it's indexing off of the trigger and starts negative
    t = []
    for i in time:
        delta = abs(min(time))
        ti = i + delta
        t.append(ti)
    

    N = len(volts)
    
    fft = np.fft.fft(volts)
    #sum squared amplitude normalization
    p = (1/N)*abs(fft**2.)
    #sample rate - 
    s_space = 1./(5*10**9) ##changed from 2.5
    fr = np.fft.fftfreq(len(volts),s_space)
				
    return p,fr,fft 
    
def Gauss(A,s,b,x):
    # return a Gaussian with given range of x values
	# s = width, A = amplitude, b = center (shift from zero)
    sig2 = 2*s*s
    G = A*np.exp(-(x-b)**2/sig2)
    return G
    
def zero_offset(volts):
    #subtract DC offset to keep pulse centered at zero
    avg = np.mean(volts)
    volts0 = volts - avg
    return volts0
    
    
def computeTF(inF,inFFT,outF,outFFT):
    """
    Generating the transfer function is always the same, so it gets its own class
    Comparable to "deconvwnr()" in matlab (nearly identical)
    
    It has a bunch of other things you can uncomment to use, but most end up being dumb
    """
        #then generate the transfer function!
        #In(f)*H(f) = Out(f) -> H(f) = Out(f)/In(f)
    tfFFT = np.divide(outFFT,inFFT)
    tfF = inF
    
    tfY = tf.fftw.irfft(tfFFT)        

#    print tfY

    tfX = np.arange(0,len(tfY))*(1./(2.*tfF[-1]))


    tfF,tfFFT = tf.genFFT(tfX,tfY)
    
    return tfX,tfY,tfF,tfFFT
				
			

def pulse_plot(time,volts,G,winV):
	""" Plot pulse, window  """
	plt.figure(1,figsize=(8,6))
	ax = plt.subplot()
	plt.plot(time*10**9, volts,linewidth=2,label='Unfiltered impulse')
	plt.hold(True)
	plt.plot(time*10**9, G,'k:',linewidth=2,label='Gaussian window')
	plt.xlim(80,1100)
	plt.ylim(min(volts)-0.05,max(volts)+0.05)
		
		
	plt.hold(True)
	plt.plot(time*10**9,winV,linewidth=2,label='Windowed impulse')
	#plt.xlim(850,1000)
	#plt.ylim(-max(winV)-0.01,max(winV)+0.01)
	
	plt.title('Windowing of recieved voltage at '+ str(angle) + ' Degrees',fontsize=20)
	plt.xlabel('Time (ns)',fontsize=18)
	plt.ylabel('Voltage (V)',fontsize=18)
	plt.legend(loc='upper right',fontsize=16)
		
	ax.tick_params(axis='both', which='major', labelsize=14, width=2)
	ax.tick_params(axis='both', which='minor', labelsize=12, width=2)
	ax.spines['top'].set_linewidth(3)
	ax.spines['right'].set_linewidth(3)
	ax.spines['bottom'].set_linewidth(3)
	ax.spines['left'].set_linewidth(3)
		
	saveopt = input("Do you want to save the plot? 'Y/N'")
	if saveopt == "Y":
		fname = ('bicone_pulse_' + str(angle) + 'deg.png')
		plt.savefig(fname,format='png',dpi=500)
		print ('Saved as : ' + str(fname)) 
		


def PSD_plot(F, P, Fwin, Pwin):
	plt.figure(2,figsize=(10,8))
	ax = plt.subplot()
	plt.plot(F*10**-9,10*np.log10(P),'r',linewidth = 1,label='PSD before windowing')

			
	plt.hold(True)
	plt.plot(Fwin*10**-9,10*np.log10(Pwin),'b',linewidth = 1,label='PSD after windowing')
	plt.xlim(0,.4) #cut to range
	plt.legend(loc='lower left',fontsize=20)
	title = ('Power Spectral Density at' + str(angle) + ' Degrees')
	plt.title(title,fontsize=30)
	plt.xlabel('Frequency (GHz)',fontsize=28)
	plt.ylabel('Power/Frequency (dB/Hz)',fontsize=28)
			
			
	ax.tick_params(axis='both', which='major', labelsize=18, width=2)
	ax.tick_params(axis='both', which='minor', labelsize=18, width=2)
	ax.spines['top'].set_linewidth(3)
	ax.spines['right'].set_linewidth(3)
	ax.spines['bottom'].set_linewidth(3)
	ax.spines['left'].set_linewidth(3)
			
	saveopt = input("Do you want to save the plot? 'Y/N'")
	if saveopt == "Y":
		fname = ('bicone_PSD_' + str(angle) + 'deg.png')
		plt.savefig(fname,format='png',dpi=500)
		print ('Saved as : ' + str(fname)) 
		

def 	effh_plot(freq,h,angle):
	lab = (str(an)+ ' Degrees')
	plt.figure(10,figsize=(8,6))
	plt.semilogy(freq*10**-6,h,label = lab)
	ax = plt.subplot()
	plt.xlim([0,1*10**3])
	#plt.ylim([0,0.06])
	plt.xlabel('Frequency (MHz)',fontsize = 20)
	plt.ylabel('Effective Height (m)',fontsize = 20)
	#plt.title('Bicone effective height at '+ str(angle) + ' degrees',fontsize = 22)
	major_ticks = np.arange(0, 1000, 100)                                              
	minor_ticks = np.arange(0, 1000, 50)                                               
		
	ax.set_xticks(major_ticks)                                                       
	ax.set_xticks(minor_ticks, minor=True) 
	
	ax.tick_params(axis='both', which='major', labelsize=18, width=2)
	ax.tick_params(axis='both', which='minor', labelsize=16, width=2)
	ax.spines['top'].set_linewidth(3)
	ax.spines['right'].set_linewidth(3)
	ax.spines['bottom'].set_linewidth(3)
	ax.spines['left'].set_linewidth(3)
	
#	saveopt = input("Do you want to save the plot? 'Y/N'")
#	if saveopt == "Y":
#		filename = ('eff_height'+ str(angle) + '.png')
#		plt.savefig(filename,format='png',dpi=500)
#		print ('Saved as : ' + str(filename)) 


	

def get_pulse(pulseplot=False,PSDplot=False):
	""" Import transmitted pulse data, apply Gaussian window, take FFTs """

	fname3 = '/Users/Caroline/Desktop/Antenna/Bicone_data_0622/FID_F0512_63dBatten_500averages_5GSPS_10kRL000_Ch1.csv'
	time,unatt_volts = get_data(fname3)
	Tvolts = unatt_volts
	
	#Gaussian window
	A = 1.
	s = 0.5
	x = time*10**9
	GP = Gauss(A,s,0,x)
	
	#windowd pulse
	pulse_win = GP*Tvolts

	#Pulse FFT
	PT_raw,FT_raw,fftT_raw = FFT(time,Tvolts)

	# pulse fft after windowing
	PTwin,FTwin,FFTTwin = FFT(time,pulse_win)
	
	if pulseplot: 
		figure6 = pulse_plot(time,Tvolts,GP,pulse_win)
		return figure6
	if PSDplot:
		figure7 = PSD_plot(FT_raw, PT_raw, FTwin, PTwin)
		return figure7
		
	return PTwin, FTwin, FFTTwin
	

def bicone_response(angle,pulseplot=False,PSDplot=False):
	""" Import recieved pulse data, apply Gaussian window, take FFTs """
	
	# import data from filename
	if angle == 0:
		fname = '/Users/Caroline/Desktop/Antenna/Bicone_data_0622/CoPol_500averages_00' + str(angle) + 'deg_Ch1.csv'
								
	elif angle <= 90:
		fname = '/Users/Caroline/Desktop/Antenna/Bicone_data_0622/CoPol_500averages_0' + str(angle) + 'deg_Ch1.csv'
												
	else:
		fname = '/Users/Caroline/Desktop/Antenna/Bicone_data_0622/CoPol_500averages_' + str(angle) +'deg_Ch1.csv'

	t_raw,volts_raw = get_data(fname)
	Rtime = t_raw + abs(np.min(t_raw))

	centered_volts = zero_offset(volts_raw)
	
	Rvolts = centered_volts
	
	peakind = np.argmax(abs(Rvolts))
	
	centertime = Rtime[peakind]
	## autocorrelation to find peak
	# find abs. max
	
	
	# apply the Gaussian window
	A = 1.
	s = 13.5
	x = Rtime*10**9
	b = (centertime*10**9)
	G = Gauss(A,s,b,x)        
	
	winV = G*Rvolts

	#take the FFTs
	Pr, Fr, FFTr = FFT(Rtime,Rvolts)
 
	# after windowing
	PRwin, FRwin, FFTRwin = FFT(Rtime,winV)
	
	if pulseplot:
		figure8 = pulse_plot(Rtime,Rvolts,G,winV)
		return figure8
		
	if PSDplot:
		figure9 = PSD_plot(Fr, Pr, FRwin, PRwin)
		return figure9
			
	return PRwin, FRwin,FFTRwin
	
def eff_height(FTwin, FFTTwin, FRwin, FFTRwin, angle, plot=False):
	# takes frequency & FFT of transmit pulse (T) and recieved pulse (R)

	inF = FTwin
	inFFT =  FFTTwin
	outF =  FRwin
	outFFT = FFTRwin

	tfX,tfY,tfF,tfFFT = computeTF(inF,inFFT,outF,outFFT)
	
	
	# attenuation correct
	transfn_raw = np.fft.fftshift(tfFFT)
	freq = np.fft.fftshift(inF)
	Freq_mhz = freq*10**-6
	
	#correcting for attenuation
	Rx_cable_attp100 = (0.242080 * np.sqrt(Freq_mhz)) + (0.000330 * Freq_mhz) #from cable spec sheet 
		#dB per 100 ft - the 0.25 -> 25ft.
	Rx_cable_att = 3. + (0.25 * Rx_cable_attp100)
	
	pulse_att = 43.
	
	dBloss = 10**((Rx_cable_att - pulse_att)/20.)
	transfn = dBloss * transfn_raw
	
	#introduce our constants	
	ld = scipy.constants.c/freq 
	d = 54.31 #distance b/w antennas
	R = 50 #Ohms
	Z0 = 377 #Ohms - impedence of space
	
	Pratio = np.absolute(np.square(transfn)) 
	A = (d*ld)*np.sqrt(Pratio)

	h = np.sqrt((2*R*A)/Z0)
	

	
	if plot:
		fig10 = effh_plot(freq, h, angle)
		return fig10
	return h
	
#%% for single angle

#PTwin, FTwin, FFTTwin = get_pulse(pulseplot=False,PSDplot=False)

#angle = input("What angle?")
angle = 0
PRwin, FRwin, FFTRwin = bicone_response(angle,pulseplot=False,PSDplot=False)
#
h_e = eff_height(FTwin, FFTTwin, FRwin, FFTRwin, angle, plot=True)
print h_e

#%% eff height plots for all angles

angle = np.arange(0,225,45)
print angle

# 
#h = [];
#leg = [];
#j = 0;
#for i in range(0,len(angle)):
#	an = angle[i]
#		
#	PRwin, FRwin, FFTRwin = bicone_response(an,pulseplot=False,PSDplot=False)
#
#	hj = eff_height(FTwin, FFTTwin, FRwin, FFTRwin, an, plot=True)
#	plt.hold(True)
#	plt.title('')
#	h.append(hj)
#	j = j + 1
#plt.legend(loc = 'best')
##plt.savefig('EffHeight_angles.png',format='png',dpi=500)
#

	




	
