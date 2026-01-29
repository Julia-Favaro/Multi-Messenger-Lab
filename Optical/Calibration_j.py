import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from astropy.io import fits
from astropy.stats import sigma_clipped_stats, SigmaClip
from astropy.modeling import models, fitting

from photutils.detection import DAOStarFinder
from photutils.background import Background2D,MedianBackground
from photutils.centroids import centroid_com
from photutils.aperture import CircularAperture,CircularAnnulus, ApertureStats,aperture_photometry
from photutils.profiles import RadialProfile

class OptCCDObs:
    """ This class will:
    -calculate the master bias, master dark and master flat.
    -calibrate light frames and find the median object frame.
    -eliminating the background level and estimating background noise of an image
    -find the centroid and FWHM of the image with an iterative process from a rough guess from input
    -perform aperture photometry in a circular aperture around the source
    -calculates the local background in annulus aperture around the source
    -performes also the PSF photometry with a 2D Gaussian function fit
    -calculates the apparent instrumental magnitude and the costant of calibration
    -calculates the signal to noise ratio"""
    
    def __init__(self,name_camera,resolution,name_source,bias_dir,dark_dir,flat_dir, light_dir):
        """Constructor. Initializing name of the CCD camera and of the source
        and its data directories."""
        self.camera=name_camera
        self.source=name_source
        self.bias_dir=bias_dir
        self.dark_dir=dark_dir
        self.flat_dir=flat_dir
        self.light_dir=light_dir
        self.resolution=resolution

    def median_image(self,data_dir):
        """Opens the files in the directory (given its path name)
        and stacks images.
        median_image:result of median image filtering"""
        files = sorted(os.listdir(data_dir))
        single_ccd_exposure_filename=[]
        for file in files:
            if os.path.isfile(os.path.join(data_dir, file)): #eliminates checkpoints
                single_ccd_exposure_filename=np.append(single_ccd_exposure_filename, file) 
        
        single_image=[fits.getdata(os.path.join(data_dir,frame), ext=0) for frame in single_ccd_exposure_filename]
        median_image = np.asarray(np.median(single_image,axis=0),float)
        return median_image
    
    def removing_glob_bkg(self,data,var_bkg=False):
        """If considering an uniform background (default), substracts the background
        level, calculated over sigma clipped image, that doesn't take into
        account all pixels above a threshold of 3 sigma from the median
        over 5 consecutive iterations.
        If background varies, performs a more elaborated 2D background estimation 
        still on the basis of 3 sigma threshold."""
        if var_bkg:
            sigma_clip = SigmaClip(sigma=3.)
            bkg_estimator = MedianBackground()
            bkg = Background2D(data,(25, 25), filter_size=(3, 3),sigma_clip=sigma_clip,bkg_estimator=bkg_estimator)
            bkg_subtracted_image = data - bkg.background
            bkg_mean, bkg_median, bkg_std = sigma_clipped_stats(bkg_subtracted_image)
            print('Background noise', bkg_std)
        else:
            bkg_mean, bkg_median, bkg_std = sigma_clipped_stats(data, sigma=3.0)
            print('Global background on overall image', bkg_median,'$\pm$ Noise', bkg_std)
            bkg_subtracted_image=data-bkg_median
        return bkg_std,bkg_subtracted_image
    
    def finding_centroid_fwhm(self,xcen,ycen,radius,data):
        """Performs initial aperture photometry statistics with data provided by user 
        by looking at the graph in order to find the centroid and FWHM. 
        Using a couple iterations ApertureStats to get a better estimation. 
        Determinating FWHM from radial profile of star."""
        #initial guess for centroid and radius is exagerated
        #on purpose usually
        aperture = CircularAperture((xcen, ycen), radius/2) 
        aperstats = ApertureStats(data, aperture)
        #iteration for better statistics
        aperture=CircularAperture((aperstats.xcentroid, aperstats.ycentroid), aperstats.fwhm.value)
        aperstats = ApertureStats(data, aperture)
        
        #determinating the FWHM from radial profile of star
        edge_radii = np.arange(26)
        rp = RadialProfile(data,(aperstats.xcentroid, aperstats.ycentroid), edge_radii)
        print('Centroid', aperstats.xcentroid, aperstats.ycentroid)  
        print('FWHM',rp.gaussian_fwhm)
        return aperstats.xcentroid,aperstats.ycentroid,rp.gaussian_fwhm
    
    def plot_circle(self,centerx, centery,radius):
        """Self explanatory."""
        theta = np.linspace(0, 2*np.pi, 100)
        x = centerx+ radius * np.cos(theta)
        y = centery + radius * np.sin(theta)
        return x,y
        
    def plot_image(self,median_image,r_ax1,title):
        """Plots median image in black and white.
        Customizable title.
        r_ax1: subplot on which you want to plot it"""
        r_ax1.set_title(title,fontsize=15)
        r_ax1.set_xlabel("X",fontsize=15)
        r_ax1.set_ylabel("Y",fontsize=15)
        r_img = r_ax1.imshow(median_image, cmap='Grays', norm=LogNorm())
        plt.colorbar(r_img)
        
    def save_image(self, title):
        """Saves the graph as a png in results_julia directory."""
        main_dir = os.getcwd()
        results_julia_dir = os.path.join(main_dir,"results_julia")
        if not os.path.exists(results_julia_dir):
            os.mkdir(results_julia_dir)
        output_filename = os.path.join(results_julia_dir,title)
        plt.savefig(output_filename+".png")
        
    #-----------------------------------------------------------------#
    def calibrated_image(self):
        """Calculates the master bias, master dark and master flat.
        Calibrates the light images with these and gets the median frame.
        Returns raw_images if you want to plot and confront with the calibrated image."""
        master_bias=self.median_image(self.bias_dir)
        master_dark=self.median_image(self.dark_dir)
        master_flat=self.median_image(self.flat_dir)
        
        #getting all the raw_images
        files = sorted(os.listdir(self.light_dir))
        single_ccd_exposure_filename=[]
        for file in files:
            if os.path.isfile(os.path.join(self.light_dir, file)):
                single_ccd_exposure_filename=np.append(single_ccd_exposure_filename, file) 
        
        raw_images=[fits.getdata(os.path.join(self.light_dir,frame), ext=0) for frame in single_ccd_exposure_filename]
        
        #calibrating images
        cal_images=[np.divide(frame-master_dark,master_flat-master_bias) for frame in raw_images]
        self.cal_image= np.asarray(np.median(cal_images,axis=0),float)  
        return raw_images
    
    def plot_calibrated_image(self):
        """Draws a graph of the calibrated and uncalibrated image one next to the other
        in black and white (logaritmic color map)."""
        raw_images=self.calibrated_image()
        
        r_fig = plt.figure(figsize=(12,5))
        r_fig.suptitle(self.source+ " image reduction\n\n",fontsize=15)
        r_ax1 = plt.subplot(1, 2, 1)
        self.plot_image(raw_images[0],r_ax1,'Raw image')
        r_ax2 = plt.subplot(1, 2, 2)
        plt.subplots_adjust(wspace=0.4)
        self.plot_image(self.cal_image,r_ax2,'Cal. image')
        
        title=self.source + " calibration"
        self.save_image(title)
        
    def aperture_photometry(self,xcen,ycen,radius):
        """Creates a circular aperture for the source and an annulus aperture
        for the local background estimation, with radii based on the FWHM and centroid of image. 
        Performes aperture photometry on these apertures 
        and estimates the local background in the source circular aperture.
        Saves the image."""
        #finding FWHM and centroid
        bkg_std,bkg_subtracted_image=self.removing_glob_bkg(self.cal_image)
        xcentroid,ycentroid,fwhm=self.finding_centroid_fwhm(xcen,ycen,radius,bkg_subtracted_image)
        
        #using the FWHM as estimator of the radius of aperture and annulus
        aperture_radius=3*fwhm
        aperture_innersky=5*fwhm
        aperture_outersky=7*fwhm
        
        apertures = CircularAperture((xcentroid,ycentroid), r=aperture_radius)
        annulus_apertures = CircularAnnulus((xcentroid,ycentroid), r_in=aperture_innersky, r_out=aperture_outersky)
        
        #plotting the image around the centroid
        plt.figure(figsize=(7,7))
        plt.title(self.source+ " aperture photometry\n\n",fontsize=15)
        plt.xlabel("X",fontsize=15)
        plt.ylabel("Y",fontsize=15)
        bb = 10*fwhm
        plt.imshow(self.cal_image, cmap='Grays', norm=LogNorm())
        plt.xlim([np.floor(xcentroid)-bb,np.floor(xcentroid)+bb])
        plt.ylim([np.floor(ycentroid)-bb,np.floor(ycentroid)+bb])
        plt.colorbar()
        
        #plotting aperture, annulus and centroid
        plt.plot(xcentroid,ycentroid,marker='+', markersize=10, color='purple', label='Centroid')
        x,y = self.plot_circle(xcentroid,ycentroid,aperture_radius)
        plt.plot(x, y, color='purple', alpha=0.5, label='Aperture')
        x,y=self.plot_circle(xcentroid,ycentroid,aperture_innersky)
        plt.plot(x, y, color='green', alpha=0.5, label='Annulus')
        x,y=self.plot_circle(xcentroid,ycentroid,aperture_outersky)
        plt.plot(x, y, color='green', alpha=0.5)
        plt.legend(loc='best')
    
        #performing aperture photometry on circular aperture (only source)
        #performing aperture photometry on annulus aperture (local background)
        apers = [apertures, annulus_apertures]
        phot_table_local_bkg = aperture_photometry(self.cal_image, apers)
        print(phot_table_local_bkg)
        
        #counts in circular aperture
        self.source_counts=float(phot_table_local_bkg["aperture_sum_0"])
        
        #local background estimation in the aperture
        bkg_mean=float(phot_table_local_bkg['aperture_sum_1'])/annulus_apertures.area_overlap(self.cal_image)
        bkg_sum = bkg_mean * apertures.area_overlap(self.cal_image)
        self.local_sky_bkg=bkg_sum
        print('Sky background', self.local_sky_bkg)
                
        #calculating number of pixel in the circular aperture
        self.n_pixels=apertures.area_overlap(self.cal_image)/self.resolution
        
        title=self.source+ " aperture photometry"
        self.save_image(title)
        
    def PSF_photometry(self,xcen,ycen,radius, model):
        """Removes uniform background and zooms image in a square box around the 
        centroid of side about 2 times the FWHM.
        Performes a 2D gaussian fit over the PSF and plots the residual image. 
        Saves the image, model and residuals as png."""
        #Median bkg level, assuming uniform noise
        bkg_std,bkg_subtracted_image =self.removing_glob_bkg(self.cal_image)
        
        #finding FWHM and centroid from rough guesses given
        xcentroid,ycentroid,fwhm=self.finding_centroid_fwhm(xcen,ycen,radius,bkg_subtracted_image)
        xc = int(xcentroid)
        yc = int(ycentroid)

        # Cut out smaller box around PSF
        # Generate grid of same size like box to put the fit on
        bb = 30
        box = bkg_subtracted_image[yc-bb:yc+bb,xc-bb:xc+bb]
        yp, xp = box.shape
        y, x, = np.mgrid[:yp, :xp]
        
        # Fit the model to your data (box)
        if model == 'gaussian':
            model_init = models.Gaussian2D(amplitude=1, x_mean=30, y_mean=30, x_stddev=5, y_stddev=5)
        else:
            print('Model function not supported yet!')
        fitter= fitting.LevMarLSQFitter()
        fitted_model = fitter(model_init, x, y, box)
        
        # Plot the data with the best-fit model
        plt.figure(figsize=(12,5))
        plt.suptitle(self.source+ " PSF photometry",fontsize=15)
        plt.subplot(1, 3, 1)
        plt.imshow(box)
        plt.title("Data",fontsize=15)
        plt.subplot(1, 3, 2)
        plt.imshow(fitted_model(x, y))
        plt.title("Model",fontsize=15)
        plt.subplot(1, 3, 3)
        plt.imshow(box - fitted_model(x, y))
        plt.title("Residual",fontsize=15)
        
        title=self.source+ " PSF photometry"
        self.save_image(title)
        
    def find_stars(self,xcen,ycen,radius,var_bkg=True):
        """Estimates the FWHM of a single star to use for the entire cluster.
        Uses DAOStarFinder method from Photutils to detect sources above a 5
        times the noise of image with background subtracted (default, non uniform bkg).
        Plots the Circular Apertures around these sources.
        Returns the FWHM to the multiple aperture function."""
        #background subtraction, non uniform background by default
        bkg_std,bkg_subtracted_image=self.removing_glob_bkg(self.cal_image,var_bkg)
        
        #selecting a star from the calibrated image and findes the fwhm
        xcentroid,ycentroid,fwhm=self.finding_centroid_fwhm(xcen,ycen,radius,bkg_subtracted_image)
        
        #find sources that are reasonably above the sky background
        #uses the fwhm of the star selected as an example for all stars
        daofind = DAOStarFinder(threshold=5*bkg_std,fwhm=3*fwhm)  
        sources = daofind(bkg_subtracted_image)
        
        for col in sources.colnames:  
            if col not in ('id', 'npix'):# for consistent table output in printing
                sources[col].info.format = '%.2f'
        self.dict_aper = sources
        return fwhm
    
    def multapert_photometry(self,xcen,ycen,radius):
        """Finds the sources in the calibrated image after calculating the FWHM to use
        and background noise to take account for.
        Performes multiple aperture photometry around found centroids.
        Calculates total counts in the source aperture, the local background and
        the number of pixels in the area of the source."""
        #finding the centroids of the sources and a typical value for the FWHM
        fwhm=self.find_stars(xcen,ycen,radius)

        #performing aperture photometry on circular aperture (only source)
        #performing aperture photometry on annulus aperture (local background)
        positions = np.transpose((self.dict_aper['xcentroid'], self.dict_aper['ycentroid']))
        apertures = CircularAperture(positions, r=3*fwhm)
        annulus_apertures = CircularAnnulus(positions, r_in=5*fwhm, r_out=7*fwhm)
        
        apers = [apertures, annulus_apertures]
        phot_table_local_bkg = aperture_photometry(self.cal_image, apers)
        for col in phot_table_local_bkg.colnames:  
            phot_table_local_bkg[col].info.format = '%.2f'
        
        #counts in circular apertures
        self.source_counts=np.array(phot_table_local_bkg["aperture_sum_0"],dtype=float)
        
        #local background estimation in the aperture
        aperstats = ApertureStats(self.cal_image, annulus_apertures)
        bkg_mean = aperstats.median
        bkg_sum = bkg_mean * apertures.area_overlap(self.cal_image)
        self.local_sky_bkg=bkg_sum
        
        #calculating number of pixel in the circular aperture
        self.n_pixels=apertures.area_overlap(self.cal_image)/self.resolution

        plt.figure(figsize=(7,7))
        plt.title(self.source+ " aperture photometry\n\n",fontsize=15)
        plt.xlabel("X",fontsize=15)
        plt.ylabel("Y",fontsize=15)
        plt.imshow(self.cal_image, cmap='Grays', norm=LogNorm())
        plt.colorbar()
        
        #marking location of correct sources on image
        mask=self.source_counts-self.local_sky_bkg>0
        apertures = CircularAperture(positions[mask],3*fwhm)
        apertures.plot(color='red', lw=1.5,alpha=0.5)
        print("Number of sources found in the image:", len(self.dict_aper[mask]))
        
        title=self.source+ " multiple aperture photometry"
        self.save_image(title)
        
    def magnitude(self,gain,dark_noise=0,read_noise=0,magn_exp=0,C=0):
        """Calculates instrumental apparent magnitude based on source_counts.
        If given the expected magnitude, calculates the calibration costant.
        If calibration constant is already given, you receive the cal. magnitude.
        
        For the signal to Noise ratio, it takes into account the poissonian error
        from the source counts and the local sky background. If specified it
        takes into account readout noise and dark noise from the camera.
        NOTE:The values for read_noise and dark_noise are assumed to be used
        in units of photons or electrons, not in units of counts or ADUs. 
        We also need to consider the total noise provided by all pixels in the 
        source area. To convert to photons, use the gain of the camera."""
        #instrumental calibration
        mv=-2.5*np.log10(self.source_counts-self.local_sky_bkg) 
        if magn_exp!=0:
            C=magn_exp-mv
            print('Costant calibration', C)
        else:
            mv=mv+C
        
        #converting the counts from ADUs to e-
        tot_counts=self.source_counts*gain
        sky_noise=self.local_sky_bkg*gain
        #calculating the total readout and dark noise
        readout_noise=read_noise**2*self.n_pixels
        darkcurr_noise=dark_noise*self.n_pixels
        
        noise=sky_noise+tot_counts+readout_noise**2+darkcurr_noise
        SN=tot_counts/np.sqrt(noise)
        print('Signal to Noise ratio',SN)
        dmv=2.5*np.log10(np.exp(1))*1/SN
        print('M_V',mv[~np.isnan(mv)],'\ndmv', dmv[~np.isnan(mv)] )
        return mv[~np.isnan(mv)]