# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 12:08:14 2017

@author: Caroline
"""

import matplotlib.pyplot as plt
import scipy.integrate 
import numpy as np
import sys
sys.path.append("/Users/Caroline/Desktop/code/")
import T510_look_at_run
import pandas as pd
import scipy.fftpack
import scipy.constants
import seaborn.apionly as sns
import cmath
import sys
sys.path.append('/Users/Caroline/A3ImpulseResponse')
import tfUtils as tf

#SLAC Data Analysis:
def get_height():
    H0 = [0.76,2.79,4.805,6.84,8.85,10.82,12.864]
    H = [h - 2.07 for h in H0]
    theta = math.radians(26.)
    delta_y = (1.06)*math.cos(theta)
    height = []
    for i in H:
        y = i + delta_y
        height.append(y)
    return height     

def get_times( scopes, event_ind, scope_ind, ch_ind):
	event = scopes.events[event_ind]
	times = event.waves[ch_ind].times
	return times				
				
def get_volts( scopes, event_ind, sc_ind, ch_ind, fs=5.0e9, fpass=None, window=False, dwindow=5e-9, sband_scale=1.):
    event = scopes.events[event_ind]
    # zero mean
    volts = event.waves[ch_ind].volts -  np.mean(event.waves[ch_ind].volts)
    #window if desired
    #might be able to make this make sure that times is odd for simp integra
    if( window ):
        times = get_times(scopes, event_ind, sc_ind, ch_ind)
        times, volts = filters.rectWindow2(times, volts, dwindow, dwindow)
    return volts
				

def get_power(volts,times):
    V = [v**2 for v in volts]
    integral = scipy.integrate.simps(V, times, even='last')
    power = integral/50
    return power
				
def analyze_run(run1, run2):
    # run 1 is negative, 2 is positive B field
	scope1 = T510_look_at_run.Scope(run1, '694')
	scope2 = T510_look_at_run.Scope(run2, '694')
    
	num_events1 = len(scope1.events)
	num_events2 = len(scope2.events)
	if num_events1 <= num_events2:
		num_events = num_events1
	else:
		num_events = num_events2
        
    #print num_events1, num_events2
    
	vpol_power = []
	hpol_power = []
            #get_times(scope, event ind, scind, chind)
     
	for i in range(num_events):
		times_vpol1 = get_times(scope1,i,0,0)
		#times_vpol2 = get_times(scope2,i,0,0)
		volts_vpol1 = get_volts(scope1,i,0,0)
		volts_vpol2 = get_volts(scope2,i,0,0)
           # hpol  
		times_hpol1 = get_times(scope1,i,0,1)
		#times_hpol2 = get_times(scope2,i,0,1)
		volts_hpol1 = get_volts(scope1,i,0,1)
		volts_hpol2 = get_volts(scope2,i,0,1)
      #  if times_vpol1 != times_vpol2:
       #     print "error: times"
		times_vpol = times_vpol1
		volts_vpol = (volts_vpol2 + volts_vpol1)/2 #average volts for Vpol
        #can't correct cross-pol in vpol
       #volts_vpol = volts_vpol1
    
		times_hpol = times_hpol1
		volts_hpol = (volts_hpol2 - volts_hpol1)/2 #remove cross-pol
		#cross_hpol = (volts_hpol2 + volts_hpol1)/2 #this is the cross pol
									
		vpol_P = get_power(volts_vpol, times_vpol)
		hpol_P = get_power(volts_hpol, times_hpol)
		
		vpol_power.append(vpol_P)
		hpol_power.append(hpol_P)
		if i == 50:
			returnV = volts_hpol
			returnT = times_hpol
			return returnV, returnT
		# returns voltage & time 
	
#Impolse Response Tests: 
	
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
				
def complex_square_and_phase(c):
	mod = abs(c)**2
	phase = np.angle(c)
	#phase=0
	return mod*np.cos(phase)+1j*mod*np.sin(phase)

def complex_square_root(c):
    r = abs(c)
    theta = np.angle(c)
    
    #unwrap the phases
    #theta = np.unwrap(theta)
    #theta=0
    return np.sqrt(r)*( np.cos( theta/2 ) + 1j*np.sin(theta/2) )
    
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
    p = (1/N)*complex_square_and_phase(fft)
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
		

def 	effh_plot(freq,h,an):
	lab = (str(an)+ ' Degrees')
	plt.figure(10,figsize=(8,6))
	plt.plot(freq*10**-6,h,label = lab)
	ax = plt.subplot()
	#plt.xlim([0,1*10**3])
#	plt.ylim([0,0.06])
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

	A = inFFT
	B = outFFT
	C = np.divide(B,A)
	transfn_raw = np.fft.fftshift(C)
	
	# attenuation correct

	freq = np.fft.fftshift(inF)
	Freq_mhz = freq*10**-6
	
	#correcting for attenuation
	Rx_cable_attp100 = (0.242080 * np.sqrt(abs(Freq_mhz)) + (0.000330 * abs(Freq_mhz))) #from cable spec sheet 
		#dB per 100 ft - the 0.25 -> 25ft.
	Rx_cable_att = 3. + (0.25 * Rx_cable_attp100)
	
	pulse_att = 43.
	
	dBloss = 10**((Rx_cable_att - pulse_att)/20.)
	transfn = dBloss * transfn_raw
	plt.plot(freq, transfn)
	
	#introduce our constants	
	ld = scipy.constants.c/freq 
	d = 54.31 #distance b/w antennas
	R = 50 #Ohms
	Z0 = 377 #Ohms - impedence of space
	
	Pratio = complex_square_and_phase(transfn)
	A = (d*ld)*complex_square_root(Pratio)

	h = complex_square_root((2*R*A)/Z0)
	print h[-10:-1]
	print h[0:10]
	
	#convert nan to 0s 
	
	
	if plot:
		fig10 = effh_plot(freq, abs(h), angle)
		
	h = np.nan_to_num(h)
	h = np.fft.ifftshift(h)
	return h
				
				
#%%
#HPOL - 				
H0 = [0.76,2.79,4.805,6.84,8.85,10.82,12.864]
H_Height = [h-2.07 for h in H0]

run_319 = '/Users/Caroline/Desktop/Antenna/VHF_data/run_319_scope_694_date_2_15_2014_time_3_43_PM.lvm'
run_320 = '/Users/Caroline/Desktop/Antenna/VHF_data/run_320_scope_694_date_2_15_2014_time_3_52_PM.lvm'

#Hvolts, Htimes = analyze_run(run_319, run_320)

#convert from list to array
#meas_v = np.asarray(Hvolts)
#meas_t = np.asanyarray(Htimes)
#print type(meas_v)
#print meas_t[0:10]

fig = plt.figure(1)
plt.plot(meas_t*10**9, meas_v)


#%%

#%% Get the electric field sims: 

#### this reads the file as a pandas dataframe

df = pd.read_csv("~/Desktop/bicones/Mag_2400A_Run2_Pol_01ns_2_413_downsampled_5Gsps.csv")
# this part makes the complex numbers
#print df.head()
df['fft_Eh'] = df['fft_Eh_real'] + 1j * df['fft_Eh_imag']
df['fft_Ev'] = df['fft_Ev_real'] + 1j * df['fft_Ev_imag']
#print df.head()
#### this gives you pandas series
#print df.freq
#print df.fft_Eh
#print df.fft_Ev
### this is to get numpy arrays

plt.plot(df.freq.values,abs(df.fft_Eh.values ))
print df.freq.values[0:10]
#return df.freq.values, df.fft_Eh.values, df.fft_Ev.values


#%% for 90

PTwin, FTwin, FFTTwin = get_pulse(pulseplot=False,PSDplot=False)

#angle = input("What angle?")
angle = 0
color = '#006EBF'
PRwin, FRwin, FFTRwin = bicone_response(angle,pulseplot=False,PSDplot=False)

h_e = eff_height(FTwin, FFTTwin, FRwin, FFTRwin, angle, plot=True)


#figure(2)
#plt.plot(df.freq.values, h_e)


#%% convolution!
	#print type(h_e)
	#print type(df.fft_Eh.values)

FFT_sim_unshifted = np.multiply(h_e, df.fft_Eh.values)

	#print FFT_sim_unshifted[0:10]
FFT_sim = np.fft.fftshift(FFT_sim_unshifted)
freq = np.fft.fftshift(df.freq.values)
plt.plot(df.freq.values, abs(FFT_sim_unshifted))

# take the ifft to get back to time-domain signal
sim_volts = np.fft.ifft(FFT_sim_unshifted)
delta_t = 1./(5*10**9)
N = len(sim_volts)
t = np.arange(0,N*delta_t,delta_t)

#dt = .2ns
fig2 = plt.figure(2,figsize=(8,6))
plt.plot(meas_t*10**9-75, meas_v,label='Measured Voltage')
plt.plot((t*10**9)-1,sim_volts,'r', linewidth=2,label='Simulated Voltage')
plt.xlim([0,100])
#the subtractions are to get the graphs to line up

plt.xlabel('Time (ns)',fontsize = 22)
plt.ylabel('Voltage (V)',fontsize = 22)
plt.legend(loc='best',fontsize=16)
#plt.title('Bicone effective height at '+ str(angle) + ' degrees',fontsize = 22)
major_ticks = np.arange(0, 100, 20)                                              
minor_ticks = np.arange(0, 100, 10)
                                               
ax = plt.subplot()
ax.set_xticks(major_ticks)                                                       
ax.set_xticks(minor_ticks, minor=True) 

ax.tick_params(axis='both', which='major', labelsize=18, width=2)
ax.tick_params(axis='both', which='minor', labelsize=16, width=2)
ax.spines['top'].set_linewidth(3)
ax.spines['right'].set_linewidth(3)
ax.spines['bottom'].set_linewidth(3)
ax.spines['left'].set_linewidth(3)

filename = 'HPOL_volts_simcomparison.png'
plt.savefig(filename,format='png',bbox_inches='tight',dpi=500)



## expect to be within factor of 2 based off of reflections 
# E sim smaller by 2 that recieved V






