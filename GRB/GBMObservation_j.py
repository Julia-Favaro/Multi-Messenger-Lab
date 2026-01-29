import os 
import sys

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.table import Table
import pandas as pd

import scipy.signal as signal
import inspect

class GRBAnalysis:
    '''This is class runs time domain analysis of a GRB event detected by GBM of the Fermi Mission.
    
    Attributes:
    -self.GRB=string, name of the event
    -self.trig=float, trig time in MET 
    -self.TTE=path string, Time Tagged Event directory
    -self.trigdet=array of floats, NaI detectors that have activated the trigger
    -self.detected=array of filenames strings, NaI and BGO detectors that have seen the event
    -self.PAT=array of floats, counts for cumulative light curve
    -self.times= array of floats, bin edges of the time bins around the event 
    -self.phot= array of floats, bin counts around the event time with background removal
    
    Methods summary:
    -plots hist of all light curves, in function of time and energy for the detectors who have triggered the observation
    -finds which of the NaI and BGO detected the event
    -plots the cumulative background and restricts to the event time interval
    -calculates the T90  
    -identifies and plots peaks in some detector's light curves
    -saves images in the correct folders'''
    
    def __init__(self,name,trigger,TTE_dir):
        '''Initialization of fundamental information of the GRB.
        Target name, trigger in MET time
        and Time Tagged Event folder for all GBM detectors.'''
        self.GRB=name
        self.trig=trigger
        self.TTE=TTE_dir
        
    def plot_hist(self,bins,hist,acceptance=5):
        """Plots histogram for cumulative light curve. 
        bins: array, edges of bins calculated by another function
        hist: array, counts calculated by another function
        acceptance=float, default to 5, level of confidence about background for it to be called an event"""
        plt.figure(figsize=(12,6))
        plt.bar(bins,hist)
        plt.title('GRB light curve', fontsize=15)
        plt.xlabel('Time since T0 [s]', fontsize=15)
        plt.ylabel('Counts/Bin ', fontsize=15)
        plt.text(55.5, np.max(hist)-3000,'acceptance='+str(acceptance)+'$\sigma$', fontsize=12,bbox=dict(alpha=0.5))
        
    def save_image(self, foldername,title):
        """Saves the graph as a png
        foldername: string, subfolder name
        title: string, name for the file"""
        main_dir = os.getcwd()
        results_dir = os.path.join(main_dir,"results")
        subfolder_dir=os.path.join(results_dir,foldername)
        if not os.path.exists(subfolder_dir):
            os.mkdir(subfolder_dir)
        output_filename = os.path.join(subfolder_dir,title)
        plt.savefig(output_filename+".jpg")
    
    #---------------------------------------------------------------------------------#
    def light_curve_analysis(self,window_beginning,window_ending,trigNaI, acceptance=5):
        '''Creates a plot of the light curves of all the NaI and BGO detectors of the GBM. 
        Identifies the detectors who have seen the event.
        Finds the NaI files of the detectors that have sent the trigger.
        window_beginning,window_ending: floats, time start and stop for the evaluation of background
        acceptance: float, default to 5, level of confidence about background for it to be called an event
        Raises: save_image'''
        self.detected=[] #NaI and BGO detectors that have seen the GRB
        files = sorted(os.listdir(self.TTE))
        NAI_files=[]        
        BG0_files=[]
        for file in files:
            if "glg_tte_n" in file: #selecting NAI detectors
                NAI_files=np.append(NAI_files, file) 
            if "glg_tte_b" in file: #selecting BGO detectors
                BG0_files=np.append(BG0_files, file) 
                
        fig=plt.figure(figsize=(15,15))
        plt.title('NaI detectors light curves')
        ii=0
        for file in NAI_files:
            input_filename = os.path.join(self.TTE,file)
            data=fits.getdata(input_filename, ext=2)
            n = Table(data)
            dfn= n.to_pandas()
            t=dfn['TIME']
            dt= t-np.min(t) #time from initial observation

            #light curve
            ax=plt.subplot(6,2,ii+1)
            bin_heights, bin_counts,_=plt.hist(dt,bins=1500)
            plt.ylim([0,np.max(bin_heights)+200])
            ax.set_title('n%d' %(ii), fontsize=12)
            ax.set_ylabel('Counts/Bin ', fontsize=12)
            
            #adjusting subplots
            if ii==11:
                ax.set_xlabel('Time since T0 [s]', fontsize=12)
            elif ii==10: 
                ax.set_xlabel('Time since T0 [s]', fontsize=12)
                plt.subplots_adjust(wspace=0.4)
            elif ii==0:
                 plt.subplots_adjust(wspace=0.4)
            elif (ii+1) % 2 == 0 and ii!=1:
                plt.subplots_adjust(hspace=0.4)
            else:
                plt.subplots_adjust(wspace=0.4)
                plt.subplots_adjust(hspace=0.4)

            # selecting window for computing background
            start = bin_counts>window_beginning 
            stop = bin_counts<window_ending
            mask_dt = np.logical_and(start, stop)

            bg=np.mean(bin_heights[mask_dt[:len(mask_dt)-1]]) #background
            err_bg=np.std(bin_heights[mask_dt[:len(mask_dt)-1]]) # background error
            print('n %d bg %f \pm %f' %(ii,bg,err_bg))

            mask_event=bin_heights>acceptance*err_bg+bg

            if len(bin_heights[mask_event])!=0:
                self.detected=np.append(self.detected,file)
            ii+=1
        self.save_image('light_curve','NAI_light_curve_'+str(acceptance)+'sigma')
        
        fig=plt.figure(figsize=(12,6))
        plt.title('BGO detectors light curves')
        ii=0
        for file in BG0_files:
            input_filename = os.path.join(self.TTE,file)
            data=fits.getdata(input_filename, ext=2)
            n = Table(data)
            dfn= n.to_pandas()
            t=dfn['TIME']
            dt= t-np.min(t) #time from start of observation

            #light curve
            ax=plt.subplot(1,2,ii+1)
            bin_heights, bin_counts,_=plt.hist(dt,bins=1500)
            plt.ylim([0,np.max(bin_heights)+200])
            ax.set_title('b%d' %(ii), fontsize=12)
            ax.set_ylabel('Counts/Bin ', fontsize=12)
            
            #adjusting subplots
            ax.set_xlabel('Time since T0 [s]', fontsize=12)
            if (ii+1) % 2 != 0:
                plt.subplots_adjust(wspace=0.4)
                
            # selecting window for computing background
            start = bin_counts>window_beginning 
            stop = bin_counts<window_ending
            mask_dt = np.logical_and(start, stop)

            bg=np.mean(bin_heights[mask_dt[:len(mask_dt)-1]]) #background
            err_bg=np.std(bin_heights[mask_dt[:len(mask_dt)-1]]) # background error
            print('b %d bg %f \pm %f' %(ii,bg,err_bg))

            mask_event=bin_heights>acceptance*err_bg+bg

            if len(bin_heights[mask_event])!=0:
                self.detected=np.append(self.detected,file)
            ii+=1
        self.save_image('light_curve','BG0_light_curve_'+str(acceptance)+'sigma')
        
        print('Detectors that have seen the event', len(self.detected))
        print('Detectors NaI',trigNaI,'have sent the trigger')
        self.trigdet=[]
        for ii in trigNaI:
            self.trigdet=np.append(self.trigdet,NAI_files[ii]) #NaI detectors of the trigger
        
    def cumulative_light_curve(self,window_beginning,window_ending,acceptance=5):
        '''Plots the cumulative light curve, with all the photons from the 
        detectors with a relevant number of counts from the event.
        window_beginning,window_ending: floats, time start and stop for the evaluation of background
        acceptance: float, default to 5, level of confidence about background for it to be called an event
        raises: save_image, plot_hist'''
        self.pat=[] #photon arrival time 
        for file in self.detected:
            input_filename = os.path.join(self.TTE,file)
            data=fits.getdata(input_filename, ext=2)
            n = Table(data)
            dfn= n.to_pandas()
            t=dfn['TIME']
            dt= t-np.min(t) #time from start of observation
            self.pat=np.append(self.pat,dt)
            
        self.trig=self.trig-np.min(t) #rescaling the trig time just for easier plots

        #background dell'istogramma cumulativo
        hist, bins=np.histogram(self.pat,bins=1500)
        start = bins>window_beginning 
        stop = bins<window_ending
        mask_dt = np.logical_and(start, stop)

        bg=np.mean(hist[mask_dt[:len(mask_dt)-1]])
        err_bg=np.std(hist[mask_dt[:len(mask_dt)-1]])
        print('bg %f \pm %f' %(bg,err_bg))

        adj_hist=hist-bg #bg removal
        mask_event=hist>acceptance*err_bg+bg 
        #in order to match the mask on hist, we need to consider one less bin edge
        self.time= bins[:-1][mask_event]
        self.phot= adj_hist[mask_event]
        
        #plots the cumulative GRB
        self.plot_hist(self.time,self.phot,acceptance)
                                   
    def calculate_T90(self):
        '''Computing T90, as the difference in time bins that corresponde to the 5% and 95% 
        of the total photon count. Calculates uncertainty from the time resolution of histogram.
        Plots the cumulative light curve with the start and end point for the T90.
        Checks what is the total error if considering also the poissonian and maximum error on calculation.
        raises: plot_hist,save_image and error_on_T90'''
        photons_tot = np.sum(self.phot)
        start_count=photons_tot*5/100
        stop_count=photons_tot*95/100
        print('Total counts above background = %.2f' %photons_tot)
        print('Five percent of counts = %.2f' %start_count)
        print('Ninety five percent of counts = %.2f' %stop_count)
        
        dt=300/1500
        print('Time resolution=%.2f' %dt)
        
        #find the bin that corresponds to 5% of the photon count
        for ii in range(len(self.phot)):
            somma=np.sum(self.phot[:ii])
            cond1=somma<start_count+1000 #define a maximum range of variability
            cond2=somma>start_count-1000
            cond=np.logical_and(cond1,cond2)
            if cond==True:
                start_time = self.time[ii]
                print('Start time from trigger= %f' % (start_time-self.trig))
                break

        #find the bin that corresponds to 95% of the photon count
        for jj in range(ii,len(self.phot)):
            somma=np.sum(self.phot[:jj])
            cond1=somma<stop_count+500 #steeper slope that we have to take to account to
            cond2=somma>stop_count-500
            cond=np.logical_and(cond1,cond2)
            if cond==True:
                stop_time = self.time[jj]
                break
                
        t_90 = stop_time - start_time
        print('T90 is %.2f +/- %.2f s' % (t_90, np.sqrt(2)*dt))
        if t_90>2:
            print('Long GRB!')
        else:
            print('Short GRB!')
        
        #checks if there really is a difference between just the time resolution and other 
        #uncertainty computing methods
        self.error_on_T90(start_count,np.sum(self.phot[:ii]),self.phot[ii],stop_count,np.sum(self.phot[:jj]),self.phot[jj],dt)
        
        #plots the cumulative GRB
        self.plot_hist(self.time,self.phot)
        plt.axvline(x=start_time, linestyle='dashed', color='orange', alpha=0.8, label='5% of counts')
        plt.axvline(x=stop_time, linestyle='dashed', color='red', alpha=0.5, label='95% of counts')
        plt.axvline(x=self.trig, linestyle='dashed', color='green', alpha=0.5, label='trigger time')
        plt.legend(loc='best', fontsize=12)
        self.save_image('light_curve','cum_light_curve_'+str(5)+'sigma')
        
    def error_on_T90(self,exp_start,meas_start,first_bin_phot,exp_stop,meas_stop,last_bin_phot,dt):
        '''Various way to determinate the error on the T90 should take into account the poissonian error
        on the photon counts and any computational error.
        exp_start: float, 5% of the total photon count, what we expect the first bin to be
        meas_start: float, the bin closest to the expected first one
        first_bin_phot: float, height of the first bin
        exp_stop: float, 95% of the total photon count, what we expect the last bin to be
        meas_stop: float, height of the last bin
        last_bin_phot: float, the bin closest to the expected last one
        dt: float, time resolution'''
        print('---------------------------------------------')
        #maximum error (computation method on the binned Time)
        print('With maximum error taken into account')
        err_start= abs(exp_start-meas_start) /np.sum(self.phot) *max(self.time)
        err_stop= abs(exp_stop-meas_stop) /np.sum(self.phot) *max(self.time)

        err_tot_start = np.sqrt(err_start**2 + dt**2)
        err_tot_stop = np.sqrt(err_stop**2 + dt**2)
        print('error t90= %.2f' %(np.sqrt(err_tot_start**2+err_tot_stop**2)))
        
        #statistical error (poissonian error of the start and stop binned photon count)
        print('With poissonian error taken into account')
        err_start = np.sqrt(first_bin_phot) /np.sum(self.phot) *max(self.time)
        err_stop = np.sqrt(last_bin_phot) /np.sum(self.phot) *max(self.time)
              
        err_tot_start = np.sqrt(err_start**2 + dt**2)
        err_tot_stop = np.sqrt(err_stop**2 + dt**2)
        print('error t90 = %.2f' %(np.sqrt(err_tot_start**2+err_tot_stop**2)))
        
    def light_curve_energy_analysis(self,range1_i, range2_i,range3_i, range3_f, event):
        '''Selecting a triggering NaI and BGO detector.
        Calculates counts in different energy ranges.
        rangex_i: float, minimum energy of the x energy range
        range3_f: float, maximum energy
        event: list of integer, time range we consider for the plot
        Raises: plot_energy_lc and save_image'''
        #selecting a NaI triggering and BGO detector filename
        files = sorted(os.listdir(self.TTE))
        BGO_files=[]
        for file in files:
            if "glg_tte_b" in file: #selecting BGO detectors
                BGO_files=np.append(BGO_files, file) 
    
        Ndetfile=self.trigdet[0]
        Bdetfile=BGO_files[0]
        
        #NaI energy dataframe
        input_filename = os.path.join(self.TTE,Ndetfile)
        data_energy = fits.getdata(input_filename, ext=1) 
        energy_bounds = Table(data_energy)
        df = energy_bounds.to_pandas()
        channel = df['CHANNEL']
        e_min = df['E_MIN']
        e_max = df['E_MAX']
        
        #selecting channels in the wanted energy range
        mask= np.logical_and(e_min>range1_i,e_max<range2_i)
        wanted_channels1=channel[mask]
        mask= np.logical_and(e_min>range2_i,e_max<range3_i)
        wanted_channels2=channel[mask]
        
        #NaI count dataframe
        data_counts = fits.getdata(input_filename, ext=2) 
        counts = Table(data_counts)
        df = counts.to_pandas()
        pha = df['PHA']
        time = df['TIME']

        #selecting counts from those channels
        mask_pha1=pha.isin(wanted_channels1)
        selected_time = np.array(df[mask_pha1]['TIME'])
        dt1= selected_time- np.min(selected_time)
        mask_pha2=pha.isin(wanted_channels2)
        selected_time = np.array(df[mask_pha2]['TIME'])
        dt2= selected_time- np.min(selected_time)

        #BGO energy dataframe
        input_filename = os.path.join(self.TTE,Bdetfile)
        data_energy = fits.getdata(input_filename, ext=1) 
        energy_bounds = Table(data_energy)
        df = energy_bounds.to_pandas()
        channel = df['CHANNEL']
        e_min = df['E_MIN']
        e_max = df['E_MAX']
        
        #selecting channels in the wanted energy range
        mask= np.logical_and(e_min>range3_i,e_max<range3_f)
        wanted_channels3=channel[mask]
        
        #BGO count dataframe
        data_counts = fits.getdata(input_filename, ext=2) 
        counts = Table(data_counts)
        df = counts.to_pandas()
        pha = df['PHA']
        time = df['TIME']

        #selecting counts from those channels
        mask_pha3=pha.isin(wanted_channels3)
        selected_time = np.array(df[mask_pha3]['TIME'])
        dt3= selected_time- np.min(selected_time)
        
        self.plot_energy_lc(dt1,dt2,dt3,range1_i, range2_i,range3_i, range3_f, event)
        self.save_image('light_curve','energy_range')
        
    def plot_energy_lc(self,dt1,dt2,dt3,range1_i, range2_i,range3_i, range3_f, event):
        '''Creates a plot of the light curves of a triggering NaI and BGO detector in different energy ranges.
        dt_x=counts in the x energy range
        rangex_i: float, minimum energy of the x energy range
        range3_f: float, maximum energy
        event: list of integer, time range we consider for the plot'''
        common_xlim = event
        fig, (ax1, ax2,ax3) = plt.subplots(3,figsize=(12,8),constrained_layout=True )
        fig.suptitle("Light curve in various energy ranges",fontsize=15)

        # Graph ax1
        ax1.set_ylabel('Counts/Bin ', fontsize=15)
        bin_heights, bin_counts, _ = ax1.hist(dt1,bins=1500)
        plt.ylim([np.min(bin_heights),np.max(bin_heights)+10])
        ax1.text(44.7, 95, 'Light curve n3 '+str(range1_i)+'-'+ str(range2_i)+' keV', fontsize=13,bbox=dict(alpha=0.5))
        ax1.set_xlim(common_xlim)
        
        # Graph ax2
        ax2.set_ylabel('Counts/Bin ', fontsize=15)
        bin_heights, bin_counts,_=ax2.hist(dt2,bins=1500)
        plt.ylim([np.min(bin_heights),np.max(bin_heights)+150])
        ax2.text(44.45, 1450, 'Light curve n3 '+str(range2_i)+'-'+ str(range3_i)+' keV', fontsize=13,bbox=dict(alpha=0.5))
        ax2.set_xlim(common_xlim)
        
        # Graph ax3
        ax3.set_xlabel('Time since T0 [s]', fontsize=15)
        ax3.set_ylabel('Counts/Bin ', fontsize=15)
        bin_heights, bin_counts,_=ax3.hist(dt3,bins=1500)
        plt.ylim([np.min(bin_heights),np.max(bin_heights)+150])
        ax3.text(43.17, 1050, 'Light curve b0 '+str(range3_i)+' keV-'+ str(range3_f/1000)+' MeV', fontsize=13,bbox=dict(alpha=0.5))
        ax3.set_xlim(common_xlim)
        
    def lc_peaks(self,file,event):
        '''Selecting relevant peaks in the light curve of some detectors.
        file= string, filename
        event: list of integer, time range we consider for the plot'''
        #count dataframe
        input_filename = os.path.join(self.TTE,file)
        data_counts = fits.getdata(input_filename, ext=2) 
        counts = Table(data_counts)
        df_c = counts.to_pandas()
        dt = np.array(df_c['TIME']) - np.min(df_c['TIME'])
        
        fig=plt.figure(figsize=(12,6))
        bin_heights, bin_counts,_=plt.hist(dt,bins=1500)
        
        # get bin centers
        bin_width = (np.max(bin_counts)-np.min(bin_counts))/1500
        bin_centers = bin_counts + bin_width/2
        bin_centers = bin_centers[:-1] 
        peaks,_ = signal.find_peaks(bin_heights, prominence=200) #find peaks
        plt.errorbar(bin_centers[peaks], bin_heights[peaks], fmt='.', color='red')
        print('Times of selected peaks', bin_centers[peaks])
        
        words=file.split("_")
        plt.title('Light curve '+words[2], fontsize=15)
        plt.xlabel('Time since T0 [s]', fontsize=15)
        plt.ylabel('Counts/Bin ', fontsize=15)
        plt.ylim([0,np.max(bin_heights)+200])
        plt.xlim(event)