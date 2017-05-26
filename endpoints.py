import glob
import csv
import re
import sys
import pandas as pd
import dsphelper
import filters as filter
import numpy as np
import analysis_library as al

def getEndpointWave(anty, dirc='/Users/Caroline/Desktop/bicones/VHF_Sims/', NFFT=None, window=None):
	
	headers, endpoints = getAllEndpoints(dirc, NFFT=NFFT, window=window)
	#print "back at getEndpointWave"
	
	# choose the cherenkov cone waveform
	chidx = ''
	for idx in headers.index:
	    head = headers.ix[idx]
	    if( head['anty'] == anty ):
		chidx = idx
		#print "Using Endpoints Simulation: ", head
	if chidx == '':
		return None, None
	ec = endpoints.ix[chidx]
	head = headers.ix[chidx]
	return head, ec

def getInfoFromName(fname, refract=True):
        # remove path name:
        #print 'at get info from name'
        name = re.split('\/', fname)[-1]
        m = re.search("Endpoint_([-+]?\d+)G_(\d+\.\d*)GeV_(\d+)e_([-+]?\d+)-([-+]?\d+)-([-+]?\d+)_(\d+)_Antenna(\d+).dat", name)
	
	try:        
		bfield = int(m.group(1))
		#bfield = -0.38165104 current +  0.25708565]
		current = (float(bfield) - 0.25708565) / -0.38165104 
		energy = float(m.group(2))
		electrons = int(m.group(3))
		antx = int(m.group(4))
		anty = int(m.group(5))
		antz = int(m.group(6))
		something = int(m.group(7))
		antenna = int(m.group(8))
        except:
		m = re.search("Endpoint_([-+]?\d+)A_(\d+\.\d*)GeV_(\d+)e_([-+]?\d+)-([-+]?\d+)-([-+]?\d+)_(\d+)_Antenna(\d+).dat", name)
		current = int(m.group(1))
		bfield = -0.38165104 * float(current) + 0.25708565
                energy = float(m.group(2))
                electrons = int(m.group(3))
                antx = int(m.group(4))
                anty = int(m.group(5))
                antz = int(m.group(6))
                something = int(m.group(7))
                antenna = int(m.group(8))

        return {'name':name,'bfield':bfield,'current':current, 'energy':energy, 'electrons':electrons, 'antx':antx, 'anty':anty, 'antz':antz, 'something':something, 'antenna':antenna, refract:refract}
        
def readEndpointFiles(fname, norm=131.e-12/1.602e-19/5000., NFFT=None, refract=True, window=None, vector_potential=False ):
        #print "Reading endpoint file ", fname 
        time = []
        Ex = [];Ey=[];Ez=[];
        with open(fname, 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='|')
            for row in reader:
                
                time.append(float(row[0])*1e-9)
		# 
		#first: arrival time at antenna in ns
		#second: EField component horizontal to beam axis in V/m for 500e (has to be scaled to 131pC), X
		#third: EField component vertical to beam axis in V/m   for 500e (has to be scaled to 131pC), Z
		#fourth: EField component along the beam axis  in V/m  for 500e (has to be scaled to 131pC), Y
		#
		# following conventions in SLAC standards document:
		# z vertical to beam axis
		# y along beam
		# x horizontal to beam axis
																
		#changed following Anne's simulations: 
		# z along beam
		# y vertical to beam axis
		# x horizontal to beam axis
		# sims: t, Ex, Ey, Ez
                Ex.append(norm*float(row[1]))
                Ey.append(norm*float(row[2]))
                Ez.append(norm*float(row[3]))

	Ex = np.array(Ex)
	Ey = np.array(Ey)
	Ez = np.array(Ez)
	
	header = getInfoFromName(fname, refract=refract)

        #if( vector_potential ):
	#	vp = al.VectorPotential(fnm=fname,num_particles=header['electrons'], charge_norm=131.e-12)
	#	time = vp.t*1e-9
	#	Ex = vp.Ex
	#	Ey = vp.Ey
	#	Ez = vp.Ez
	Ex = dsphelper.nanscan(Ex, before=False)
	Ey = dsphelper.nanscan(Ey, before=False)
	Ez = dsphelper.nanscan(Ez, before=False)

        if( NFFT==None ):
		NFFT = dsphelper.nextpower2(len(Ex))

	if( window != None):
		tmaxx = dsphelper.get_tmax(time, Ex)
		tmaxy = dsphelper.get_tmax(time, Ey)
		tmaxz = dsphelper.get_tmax(time, Ez)
		time, Ex = filter.rectWindow(time, Ex, tmaxx-window[0], tmaxx+window[1])
		time, Ey = filter.rectWindow(time, Ey, tmaxy-window[0], tmaxy+window[1])
		time, Ez = filter.rectWindow(time, Ez, tmaxz-window[0], tmaxz+window[1])

	# get the vpol & hpol components
	Eh = Ex
	Ev = -Ey * (np.sin(19.6*np.pi/180.)) + Ez * np.cos(19.6*np.pi/180.)

	Ex = dsphelper.nanscan(Ex, before=False)
	Ey = dsphelper.nanscan(Ey, before=False)
	Ez = dsphelper.nanscan(Ez, before=False)
	Ev = dsphelper.nanscan(Ev, before=False)
	Eh = dsphelper.nanscan(Eh, before=False)
	t,ex,fr,FEx = dsphelper.compute_fourier_pairs(time,Ex,len(Ex),NFFT)
	t,ey,fr,FEy = dsphelper.compute_fourier_pairs(time,Ey,len(Ey),NFFT)
	t,ez,fr,FEz = dsphelper.compute_fourier_pairs(time,Ez,len(Ez),NFFT)
	t,ev,fr,FEv = dsphelper.compute_fourier_pairs(time,Ev,len(Ev),NFFT)
	t,eh,fr,FEh = dsphelper.compute_fourier_pairs(time,Eh,len(Eh),NFFT)  
       
        return {'header': header, 'wave': pd.DataFrame({'time':time, 'Ex':Ex, 'Ey':Ey, 'Ez':Ez, 'Ev':Ev, 'Eh':Eh}, index=time), 'fd': pd.DataFrame({'freq':fr, 'Ex':FEx, 'Ey':FEy, 'Ez':FEz, 'Ev':FEv, 'Eh':FEh}, index=fr)}

def getAllEndpoints(dirc ='/Users/Caroline/Desktop/bicones/VHF_Sims/', NFFT=None, window=None, norm=131.e-12/1.602e-19/5000.,  vector_potential=False ):
    files = glob.glob(dirc+"/Endpoint*.dat")

    epts = []
    names = []
    headers = []
    for fname in files:   
            ep = readEndpointFiles(fname, NFFT=NFFT, window=window, norm=norm, vector_potential=vector_potential)
            epts.append(ep)
            names.append(ep['header']['name'])
            headers.append(ep['header'])
    #print 'get all Endpoints:'
    headers = pd.DataFrame(headers, index=names)
    endpoints = pd.DataFrame(epts, index=names)
    return headers, endpoints

