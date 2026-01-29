import os 

from astropy.io import fits
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy import optimize

import functions as f

class CCD_camera:
    """ This class will:
    -read images from FITS Format from a given folder and their median image.
    -calculate the master bias, master dark and master flat.
    -calculate their relative count histograms.
    -calculate mean and standard deviation of their distribution.
    -fit the histograms to a gaussian model.
    -calculate the dark current of a catalog of master darks.
    -calculate gain of a camera from flat and dark frames."""
    
    def __init__(self,name,bias_dir,dark_dir,flat_dir):
        """Constructor. Initializing name of the CCD camera and its data directories"""
        self.name=name
        self.bias_dir=bias_dir
        self.dark_dir=dark_dir
        self.flat_dir=flat_dir

    def median_image(self,data_dir):
        """Opens the files in the directory and stacks images.
        data_dir:path name
        median_image:result of median image filtering"""
        files = sorted(os.listdir(data_dir))
        single_ccd_exposure_filename=[]
        for file in files:
            if os.path.isfile(os.path.join(data_dir, file)): #eliminates checkpoints
                single_ccd_exposure_filename=np.append(single_ccd_exposure_filename, file) 
        
        single_image=[fits.getdata(os.path.join(data_dir,frame), ext=0) for frame in single_ccd_exposure_filename]
        median_image = np.asarray(np.median(single_image,axis=0),float)
        return median_image
    
    def get_average_value(self,image): 
        """Gets mean and standard deviation from a single image.
        Calculates them on the basis of a gaussian fit on the corresponding histogram."""
        bin_heights, bin_borders=np.histogram(image.flatten(),bins='auto')
        r_min=np.min(image)
        r_max=np.max(image)
        r_mean = np.mean(image)
        r_std = np.std(image)
        
        bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
        init=[np.max(bin_heights),r_mean,r_std]
        popt,pcov= optimize.curve_fit(f.model_gauss,bin_centers, bin_heights,p0=init,maxfev=10000)
        error=np.sqrt(pcov.diagonal())
        
        mean=popt[1]
        err_mean=np.sqrt(pcov.diagonal())[1]
        sigma=popt[2]
        err_sigma=np.sqrt(pcov.diagonal())[2]
        return mean,err_mean,sigma,err_sigma
    
    def plot_median(self,median_image,r_ax1,title,color):
        """Plots median image.
        Customizable color and title.
        r_ax1: subplot on which you want to plot it"""
        r_ax1.set_title(title,fontsize=15)
        r_ax1.set_xlabel("X",fontsize=15)
        r_ax1.set_ylabel("Y",fontsize=15)
        r_img = r_ax1.imshow(median_image, cmap=color, norm=LogNorm())
        plt.colorbar(r_img)
    
    def plot_hist(self,median_image,r_ax2, color, fit='False', median='False'):
        """Plots relative histogram in ADUs. 
        Plots in a text box the mean and standard deviation of this distribution.
        You can specify if you want to fit it with a gaussian model.
        Customizable color. 
        r_ax2: subplot on which you want to plot it """
        r_ax2.set_title("Histogram",fontsize=15)
        r_ax2.set_xlabel("Counts [ADU]",fontsize=15)
        r_ax2.set_ylabel("n° pixels",fontsize=15)
        plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        bin_heights, bin_borders, _=r_ax2.hist(median_image.flatten(),bins=150,alpha=0.5)
        
        #Get some statistics
        r_min=np.min(median_image)
        r_max=np.max(median_image)
        r_mean = np.mean(median_image)
        r_std = np.std(median_image)  
        r_median=np.median(median_image)
        r_ax2.set_xlim([0,2*r_mean])
        
        if median=='True': #distribution is gaussian if mean and median correspond!
            r_text = '\n'.join((
            r'mean=%.f' % (r_mean ),
            r'median=%.f' %(r_median),
            r'std=%.f' % (r_std )))
            
            props = dict(boxstyle='round', alpha=0.5)
            r_ax2.text(0.03, 0.65, r_text, transform=r_ax2.transAxes, fontsize=12,
                    verticalalignment='top', bbox=props)
        else:
            r_text = '\n'.join((
                r'mean=%.f' % (r_mean ),
                r'std=%.f' % (r_std )))

            props = dict(boxstyle='round', alpha=0.5)
            r_ax2.text(0.65, 0.95, r_text, transform=r_ax2.transAxes, fontsize=12,
                    verticalalignment='top', bbox=props)
        plt.xlim([r_min,r_max]);
        
        if fit=='True':
            bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
            init=[np.max(bin_heights),r_mean,r_std]
            popt,pcov= optimize.curve_fit(f.model_gauss,bin_centers, bin_heights,p0=init)
            error=np.sqrt(pcov.diagonal())

            x= np.linspace(r_min, r_max, 10000)
            plt.plot(x,f.model_gauss(x,*popt),color='darkorange')

            #Place in a text box the results of fit
            fit_text = '\n'.join((
                '$\mu_{ADU}$=%.f $\pm$ %.f' % (popt[1],error[1]),
                '$\sigma_{ADU}$=%.f $\pm$ %.f' % (popt[2],error[2])))

            props = dict(boxstyle='round', facecolor='darkorange', alpha=0.5)
            r_ax2.text(0.03, 0.75, fit_text, transform=r_ax2.transAxes, fontsize=12,
                    verticalalignment='top', bbox=props)
        
    def save_image(self, title):
        """Saves the graph as a png"""
        main_dir = os.getcwd()
        results_julia_dir = os.path.join(main_dir,"results_julia")
        if not os.path.exists(results_julia_dir):
            os.mkdir(results_julia_dir)
        output_filename = os.path.join(results_julia_dir,title)
        plt.savefig(output_filename+".png")
        
    #-----------------------------------------------------------------#
    def plot_master_bias(self,color='gray', fit='False'):
        """Opens the image and histogram of the master bias. 
        Then saves the image.
        Default monochromatic camera, unless otherwise specified
        You can specify if you want to fit it with a guassian model""" 
        data_dir=self.bias_dir
        image=self.median_image(data_dir)
        r_fig = plt.figure(figsize=(12,5))
        title=self.name + " Master Bias"
        r_fig.suptitle(title+'\n\n',fontsize=15)
        r_ax1 = plt.subplot(1, 2, 1)
        self.plot_median(image,r_ax1,'Image',color)
        r_ax2 = plt.subplot(1, 2, 2)
        plt.subplots_adjust(wspace=0.4)
        self.plot_hist(image,r_ax2,color, fit)
        self.save_image(title)

    def master_dark_catalog(self):
        """Creates a catalog of master dark at different exposures
        and prints their average value and standard deviations,
        after removing the mean of the master bias."""
        bias_dir=self.bias_dir
        master_bias=self.median_image(bias_dir)
        
        data_dir=self.dark_dir
        times=['0.1s', '1s', '5s', '10s', '30s', '50s', '70s', '85s', '100s']
        exposure_times=['Dark_'+time for time in times]
        
        mean=np.ones(len(times))
        err_mean=np.ones(len(times))
        std=np.ones(len(times))
        err_std=np.ones(len(times))
        
        for ii in range(len(exposure_times)):
            exposure_dir = os.path.join(data_dir, exposure_times[ii])
            master_dark_uncal=self.median_image(exposure_dir)
            master_dark_cal=master_dark_uncal-master_bias
            mean[ii],err_mean[ii],std[ii],err_std[ii]=self.get_average_value(master_dark_cal)
            print(times[ii]+'-> mean=%.2f $\pm$ %.2f [ADU], std=%.2f $\pm$ %.2f [ADU]' %(mean[ii],err_mean[ii],abs(std[ii]),err_std[ii]))
        return mean,err_mean,abs(std)
    
    def plot_master_dark_catalog(self):
        """Create a plot of all the different histograms of the
        master darks at different time exposure.""" 
        data_dir=self.dark_dir
        times=['0.1s', '1s', '5s', '10s', '30s', '50s', '70s', '85s', '100s']
        exposure_times=['Dark_'+time for time in times]
    
        cm = plt.get_cmap('gist_rainbow')  #preparing colors
        colors = cm(np.linspace(0, 1, len(exposure_times)))
        
        r_fig = plt.figure(figsize=(12,5))
        title=self.name + " Master Dark catalog"
        plt.title(title+'\n\n',fontsize=15)
        plt.xlabel("Counts [ADU]",fontsize=15)
        plt.ylabel("n° pixels",fontsize=15)
        
        bias_dir=self.bias_dir
        master_bias=np.mean(self.median_image(bias_dir))
        #plt.axvline(master_bias,linestyle='--', alpha=0.5, label='bias')
        
        for ii in range(len(exposure_times)):
            exposure_dir = os.path.join(data_dir, exposure_times[ii])
            master_dark_uncal=self.median_image(exposure_dir)
            master_dark_cal=master_dark_uncal-master_bias
            plt.hist(master_dark_cal.flatten(),bins=3000,color=colors[ii],alpha=0.5,label=times[ii])
        
        plt.xlim([-100,200])
        plt.legend(loc='best')
        self.save_image(title)
    
    def dark_current(self, resolution):
        """Plots a graph of mean values of the master bias of different exposure times.
        Fits the trend with a linear function.
        Plots the best fit line and its residuals.
        Calculates the dark current of the camera, depending on the pixel number."""
        t=np.array([0.1,1,5,10,30,50,70,85,100])
        mean, err_mean, std=self.master_dark_catalog()
        
        #plot raw data
        plt.figure(figsize=(12,6))
        plt.subplot(2,1,1)
        title="CCD Atik460 Dark Current"
        plt.title(title+"\n", fontsize=15)
        plt.ylabel("Dark count [ADU]", fontsize=15)
        plt.errorbar(t, mean, yerr=std, fmt='.', label='data',alpha=0.5)
        
        #plot fit
        popt, pcov = optimize.curve_fit(f.line,t, mean, sigma=std,absolute_sigma=True)
        xx=np.linspace(min(t), max(t), 1000)
        plt.plot(xx,f.line(xx,*popt), color='darkorange', label='fit')
        plt.legend(loc='best',fontsize=12)
        
        dof=len(mean)-len(popt)
        res=(mean-f.line(t,*popt))/std
        chisq=(res**2).sum()
        chisq_red = chisq/dof
        print('Chisq ridotto %.2f' %(chisq_red))
        
        #plot residuals
        plt.subplot(2,1,2)
        plt.subplots_adjust(hspace=0.4)
        plt.errorbar(t, res, fmt='.',alpha=0.5)
        plt.ylabel('Norm. res.',fontsize=15)
        plt.xlabel("Exposure time [s]", fontsize=15)
        plt.axhline(0.,linestyle='--', color='darkorange')  
        
        #calculates the dark current
        m=popt[1]
        dm=np.sqrt(pcov.diagonal())[1]
        self.pixels=resolution #this are attributes of the class camera
        self.dark_current=m/self.pixels
        self.ddark_current=dm/self.pixels
        
        print('dark current: %.2f $\pm$ %.2f 10^{-6} $e{^-}/s/pix$' %(self.dark_current*1e6,self.ddark_current*1e6))
        print('dark current: %.2f $\pm$ %.2f $e{^-}/s$' %(self.dark_current*self.pixels,self.ddark_current*self.pixels))
        self.save_image(title)
        
    def plot_master_flat(self,color='gray',fit='False'):
        """Opens the image and histogram of the master flat. 
        Then saves the image.
        Default monochromatic camera, unless otherwise specified""" 
        data_dir=self.flat_dir
        image=self.median_image(data_dir)
        r_fig = plt.figure(figsize=(12,5))
        title=self.name + " Master Flat"
        r_fig.suptitle(title+'\n\n',fontsize=15)
        r_ax1 = plt.subplot(1, 2, 1)
        self.plot_median(image,r_ax1,'Image',color)
        r_ax2 = plt.subplot(1, 2, 2)
        plt.subplots_adjust(wspace=0.4)
        self.plot_hist(image,r_ax2,color, median='True')
        self.save_image(title)
    
    def calculate_gain_readnoise(self,bias1_filename, bias2_filename,flat1_filename,flat2_filename):
        """Calculates gain from the mean of B1,B2,F1,F2 images 
        and the standard deviation of the difference images F1-F2 e B1-B2.
        Calculates the readout noise from the gain so determined"""
        bias_dir=self.bias_dir
        flat_dir=self.flat_dir

        #calculating B1 e B2, B1-B2
        B1=fits.getdata(os.path.join(bias_dir,bias1_filename), ext=0)
        B2=fits.getdata(os.path.join(bias_dir,bias2_filename), ext=0)
        B1=np.asarray(B1,float) #converting to float to prevent overflow
        B2=np.asarray(B2,float)
        B12=B1-B2
        
        mb1,err_mb1,sb1,err_sb1=self.get_average_value(B1)
        print('mb1 %.2f $\pm$ %.2f' %(mb1, err_mb1))
        print('sb1 %.2f $\pm$ %.2f' %(sb1, err_mb1))
        mb2,err_mb2,sb2,err_sb2=self.get_average_value(B2)
        print('mb2 %.2f $\pm$ %.2f' %(mb2, err_mb2))
        print('sb2 %.2f $\pm$ %.2f' %(sb2, err_mb2))
        mb12,err_mb12,sb12,err_sb12=self.get_average_value(B12)
        print('mb12 %.2f $\pm$ %.2f' %(mb12, err_mb12))
        print('sb12 %.2f $\pm$ %.2f' %(sb12, err_mb12))
        
        #calculating F1 e F2, F1-F2
        F1=fits.getdata(os.path.join(flat_dir,flat1_filename), ext=0)
        F2=fits.getdata(os.path.join(flat_dir,flat2_filename), ext=0)
        F1=np.asarray(F1,float) 
        F2=np.asarray(F2,float)
        F12=F2-F1
        
        mf1,err_mf1,sf1,err_sf1=self.get_average_value(F1)
        print('mf1 %.2f $\pm$ %.2f' %(mf1, err_mf1))
        print('sf1 %.2f $\pm$ %.2f' %(sf1, err_sf1))
        mf2,err_mf2,sf2,err_sf2=self.get_average_value(F2)
        print('mf2 %.2f $\pm$ %.2f' %(mf2, err_mf2))
        print('sf2 %.2f $\pm$ %.2f' %(sf2, err_sf2))
        mf12,err_mf12,sf12,err_sf12=self.get_average_value(F12)
        print('mf12 %.2f $\pm$ %.2f' %(mf12, err_mf12))
        print('sf12 %.2f $\pm$ %.2f' %(sf12, err_sf12))
        
        gain,err_gain=f.gain(mb1,err_mb1,mb2,err_mb2,sb12,err_sb12, mf1,err_mf1,mf2,err_mf2,sf12,err_sf12)
        r_noise,err_rnoise=f.readout_noise(gain,err_gain,sb12,err_sb12)
        print('gain %.3f $\pm$ %.3f $e^{-s}/ADU$' %(gain, err_gain))
        print('readout noise %.2f $\pm$ %.2f ADU' %(r_noise,err_rnoise))