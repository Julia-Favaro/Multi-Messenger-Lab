import os 
import sys

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from scipy import optimize
from scipy.signal import find_peaks

import inspect
import functions_j as f

class Radio_fitter:
    """This is class runs the analysis of HI cloud in the Milky Way and map their distribution and velocity in our galaxy.

    Attributes:
    -self.x, self.y, self.dy= array of floats, x and y axis and uncertainty over y axis
    -self.mask= boolean array, mask that identifies tails from the event
    -self.fit_model= function, fit model function 
    -self.popt, self.pcov = parameters of best fit and covariance matrix
    -self.res= array floats, normalized residuals
    -self.chisq_red= float, reduced chisquare
    -self.xlabel, self.ylabel= strings, x-y axis labels for plots
    -self.name= string, name for results text and plots
    -self.l= float, galactic longitude
    -self.dl= float, uncertainty galactic longitude
    
    Methods summary:
    -Read data from TXT Format
    -Preprocesse the data, fit to a gaussian model given and print the results
    -Plot data and model
    -Compute residuals and plot them  
    -Compute the chisquare and print the result
    -Create output file where to save the image and one for the parameters of best fit
    -Find radius and velocity of the gas with maximum velocity 
    -Find radius and velocity of all the gases"""
        
    def __init__(self,x_label, y_label):
        """Constructor. Initializing names for plots"""
        self.xlabel=x_label
        self.ylabel=y_label
    
    def read_data(self,input_name):
        """This function creates the input filename to open and redirects to
        a function to create the data frame.
        input_name: string
        raises read_txt"""
        #define work, data directories and path of the input file
        main_dir = os.getcwd() 
        data_dir = os.path.join(main_dir,"data")
        input_filename = os.path.join(data_dir,input_name)
        
        #save the name for the output file name
        self.name=input_name
        self.read_txt(input_filename)
        
    def read_txt(self,input_filename):
        """This function should open the file and create the data frame from a txt file.
        input_filename: string"""
        header_length=8
        
        #getting the Glon of this acquisition
        input_file=open(input_filename, "r")
        catalog_lines = input_file.readlines() 
        coord=catalog_lines[4]
        words=coord.split("=")
        self.l=float(words[1])
        self.dl=0.2 #calibration error 
        
        df= pd.read_csv(input_filename, sep=' ', skiprows=header_length, names=[self.xlabel, self.ylabel])
        dy_value=0.25 #calibration error percentage
        errors=df[self.ylabel]*dy_value
        df.insert(2, "dy", errors , True)
        
        self.x=df[self.xlabel]
        self.y=df[self.ylabel]
        self.dy=df["dy"]  
        
    def plot_data(self):
        """This function should create a scatterplot of the raw data."""
        plt.figure(figsize=(12,6))
        plt.title("Raw data\nGlon= %.1f° $\pm$ %.1f° spectrum" %(self.l,self.dl),fontsize=20)
        plt.errorbar(self.x, self.y, fmt='.', label='data')
        plt.xlabel(self.xlabel,fontsize=20)
        plt.ylabel("uncal. antenna temp. [K]",fontsize=20)
        plt.minorticks_on()
            
    def baseline(self,beginning, end):
        """This function defines the tails of the function. 
        beginning: float, value where the first tail ends at
        end: float, value where to second tail starts at"""
        mask_tail=(self.x<beginning) | (self.x>end)
        self.mask=np.invert(mask_tail) #mask that cuts the tails
        
    def initial_guess(self,x,y,plot=False):
        """Finds peaks in smooth data and from that creates 
        the initial values for the next fit 
        and chooses what fit function to use"""
        index,_=find_peaks(y,prominence=5)
            
        if plot==True:
            plt.errorbar(x[index],y[index],fmt='o', label='peaks')
            
        guess=[0] #constant
        for ii in index:
            guess.append(y[ii]) #amplitude
            guess.append(x[ii]) #mean
            guess.append(1) #sigma
        
        if len(index)==1:
            self.fit_model=f.model_gauss
        if len(index)==2:
            self.fit_model=f.multi_gauss2
        if len(index)==3:
            self.fit_model=f.multi_gauss3
        if len(index)==4:
            self.fit_model=f.multi_gauss4
        return guess
        
    def preprocessing(self,plot=False):
        """It computes the popt of a blind fit over 
        smooth data with moving average (windowing of 2).
        the fit focuses only on the central part of the plots and 
        takes the initial values from another method.
        Returns the popt as the gaussian initial parameters 
        for the final fit."""
        x=f.moving_average(self.x[self.mask])
        y=f.moving_average(self.y[self.mask])
        dy=f.moving_average(self.dy[self.mask])
        
        if plot==True:
            plt.figure(figsize=(12,6))
            plt.title("Preprocessing\nGlon= %.1f° $\pm$ %.1f° spectrum" %(self.l,self.dl),fontsize=20)
            plt.errorbar(x, y, yerr=dy, fmt='.', label='data')
            plt.xlabel(self.xlabel,fontsize=20)
            plt.ylabel(self.ylabel,fontsize=20)
            plt.minorticks_on()
        
        if plot==True:
            init=self.initial_guess(x,y,plot=True)
        else:
            init=self.initial_guess(x,y)
            
        popt, pcov = optimize.curve_fit(self.fit_model, x,y, sigma=dy, 
                                        p0=init, absolute_sigma=True,maxfev = 10000)
        if plot==True:
            plt.plot(x,self.fit_model(x,*popt),label='fit')
        plt.legend(loc='best',fontsize=15)
        
        return popt
    
    def fit_data(self):
        """This function creates the fit depending on the 
        fitting function given after the preprocess routine.
        raises preprocessing"""
        plt.figure(figsize=(12,6))
        plt.subplot(2,1,1)
        plt.title("Glon= %.1f° $\pm$ %.1f° spectrum" %(self.l,self.dl),fontsize=15)
        plt.errorbar(self.x, self.y, yerr=self.dy, fmt='.', label='data')
        plt.xlabel(self.xlabel,fontsize=15)
        plt.ylabel(self.ylabel,fontsize=15)
        plt.minorticks_on()
        
        init=self.preprocessing()
        self.popt, self.pcov = optimize.curve_fit(self.fit_model, self.x[self.mask], self.y[self.mask], 
                                            sigma=self.dy[self.mask], p0=init, absolute_sigma=True, maxfev = 10000)
        
        #plotting the fit
        xx=np.linspace(min(self.x[self.mask]), max(self.x[self.mask]), 1000)
        plt.plot(xx,self.fit_model(xx,*self.popt), color='red', label='fit')
        plt.legend(loc='best',fontsize=15)
        
        #computing residuals and chisquare
        dof=len(self.y[self.mask])-len(self.popt)
        self.res=(self.y[self.mask]-self.fit_model(self.x[self.mask],*self.popt))/self.dy[self.mask]
        chisq=(self.res**2).sum()
        self.chisq_red = chisq/dof
        
    def plot_res(self):
        """This funtion plots the residuals on a subplot."""
        plt.subplot(2,1,2)
        plt.subplots_adjust(hspace=0.4)
        
        plt.errorbar(self.x[self.mask], self.res, fmt='.')
        plt.ylabel('Norm. res.',fontsize=15)
        plt.xlabel(self.xlabel,fontsize=15)
        plt.xlim([min(self.x),max(self.x)]) #so that it matches the plot above
        plt.axhline(0.,linestyle='--', color='red')
        plt.minorticks_on()
                  
    def print_fit_rotation_curve(self):
        """This function prints the parameter of best fit and their uncertainties. 
        It prints also the velocity and the distance from 
        GC for the cloud with highest velocity."""
        #get the name of the parameters of the function, 
        #but remember to skip the first one which is always x
        parameters_name=inspect.getargspec(self.fit_model).args
        matrix=self.pcov
        error=np.sqrt(matrix.diagonal())
        
        print("--------------------------\nGlon: %.2f $\pm$ %.2f" %(self.l, self.dl))
        print('Parameters of best fit:')
        for ii in range(len(self.popt)):
            print(parameters_name[ii+1], "%.2f $\pm$ %.2f" %(self.popt[ii], error[ii]))
        print('Chisq_red %.2f' %self.chisq_red)
        print('V_max %.2f $\pm$ %.2f' %(self.Vabs, self.dVabs))
        print('R %.4f $\pm$ %.4f' %(self.R, self.dR))
    
    def print_fit(self):
        """This function prints the parameter of best fit 
        and their uncertainties. """
        #get the name of the parameters of the function, 
        #but remember to skip the first one which is always x
        parameters_name=inspect.getargspec(self.fit_model).args
        matrix=self.pcov
        error=np.sqrt(matrix.diagonal())
        
        print("--------------------------\nGlon: %.2f $\pm$ %.2f" %(self.l, self.dl))
        print('Parameters of best fit:')
        for ii in range(len(self.popt)):
            print(parameters_name[ii+1], "%.2f $\pm$ %.2f" %(self.popt[ii], error[ii]))
        print('Chisq_red %.2f' %self.chisq_red)
                            
    def get_velocity(self):
        """This functions prints out all the velocities 
        corresponding to the mean of the gaussians."""
        matrix=self.pcov
        error=np.sqrt(matrix.diagonal())
        
        v=[]
        dv=[]
        #first 2 parameters are the constant and first amplitude
        for ii in range(2,len(self.popt),3):
            v=np.append(v,self.popt[ii])
            dv=np.append(dv,error[ii])
        return v, dv
    
    def abs_velocity(self):
        """Calculates the absolute velocity of the HI cloud with highest velocity 
        and its uncertainty in the galactic plane
        raises get_velocity"""
        v,dv=self.get_velocity()
        l=np.radians(self.l)
        dl=np.radians(self.dl)
        V0=220
        
        #the HI cloud with highest velocity has tangential velocity
        Vr_max=np.max(v) 
        dVr_max=dv[np.argmax(v)]
        print(dVr_max)
        print(dl)
        
        V=Vr_max+V0*np.sin(l) #absolute velocity 
        dV=np.sqrt(dVr_max**2+(V0*np.cos(l)*dl)**2)
        self.Vabs=V
        self.dVabs=dV
        return V,dV
    
    def distance_GC_tangent(self):
        """Calculates the distance from GC and its uncertainty
        for the HI cloud with highest velocity"""
        R0=8.5 
        l=np.radians(self.l)
        dl=np.radians(self.dl)
        
        R=R0*np.sin(l)
        dR=R0*np.cos(l)*dl
        self.R=R
        self.dR=dR
        return R,dR
    
    def orbit_radius(self):
        """Calculates the orbit radius and its uncertainty 
        for each HI cloud
        raises get_velocity"""
        R0=8.5 
        V0=220 
        l=np.radians(self.l)
        dl=np.radians(self.dl)
        Vr,dVr=self.get_velocity()
        
        num=R0*V0*np.sin(l)
        den=V0*np.sin(l)+Vr
        R=num/den
        
        derivata_l=(R0*V0*np.cos(l))/den-(num*V0*np.cos(l))/den**2
        derivata_Vr=-num/den**2
        dR=np.sqrt((derivata_l*dl)**2+(derivata_Vr*dVr)**2)
        
        return R,dR
    
    def distance_GC(self):
        """Calculates the distance from GC and its relative uncertainty
        for each HI cloud in cartesian coordinates
        raises orbit_radius and convert_distance_GC"""
        R0=8.5 
        l=np.radians(self.l)
        dl=np.radians(self.dl)
        R,dR=self.orbit_radius()
                     
        dpiu=+np.sqrt(R**2-R0**2*np.sin(l)**2)+R0*np.cos(l)
        dmeno=-np.sqrt(R**2-R0**2*np.sin(l)**2)+R0*np.cos(l)
        
        derivata_Rpiu=R/(np.sqrt(R**2-R0**2*np.sin(l)**2))
        derivata_lpiu=-(R0**2*np.sin(l)*np.cos(l))/(np.sqrt(R**2-R0**2*np.sin(l)**2))-R0*np.sin(l)
        derivata_Rmeno=-R/(np.sqrt(R**2-R0**2*np.sin(l)**2))
        derivata_lmeno=(R0**2*np.sin(l)*np.cos(l))/(np.sqrt(R**2-R0**2*np.sin(l)**2))-R0*np.sin(l)
        
        ddpiu=np.sqrt((derivata_Rpiu*dR)**2+(derivata_lpiu*dl)**2)
        ddmeno=np.sqrt((derivata_Rmeno*dR)**2+(derivata_lmeno*dl)**2)
        
        x=[]
        dx=[]
        y=[]
        dy=[]
        xdeg=[]
        dxdeg=[]
        ydeg=[]
        dydeg=[]

        for ii in range(len(R)):        
            if dmeno[ii]<0:
                xpoint,ypoint,err_xpoint,err_ypoint=self.convert_distance_GC(dpiu[ii],ddpiu[ii])
                x=np.append(x,xpoint)
                y=np.append(y,ypoint)
                dx=np.append(dx,err_xpoint) #relative error
                dy=np.append(dy,err_ypoint)
            elif dpiu[ii]<0:
                xpoint,ypoint,err_xpoint,err_ypoint=self.convert_distance_GC(dmeno[ii],ddmeno[ii])
                x=np.append(x,xpoint)
                y=np.append(y,ypoint)
                dx=np.append(dx,err_xpoint)
                dy=np.append(dy,err_ypoint)
            else:            
                if ((33<self.l) and (self.l<48)): #nubi lontane
                    xpointdeg,ypointdeg,drelx, drely=self.convert_distance_GC(dpiu[ii],ddpiu[ii])
                    xdeg=np.append(xdeg,xpointdeg)
                    ydeg=np.append(ydeg,ypointdeg)
                    dxdeg=np.append(dxdeg,drelx)
                    dydeg=np.append(dydeg,drely)
                else: #nubi vicine
                    xpointdeg,ypointdeg,drelx, drely=self.convert_distance_GC(dmeno[ii],ddmeno[ii])
                    xdeg=np.append(xdeg,xpointdeg)
                    ydeg=np.append(ydeg,ypointdeg)
                    dxdeg=np.append(dxdeg,drelx)
                    dydeg=np.append(dydeg,drely)  
                
        return x,y,xdeg,ydeg,dx,dy,dxdeg,dydeg
    
    def convert_distance_GC(self,d,dd):
        """Converts the distance in GC in cartesian coordinates.
        Calculates the relative error
        d: float, distance
        dd: float, uncertainty over distance"""
        R0=8.5
        l=np.radians(self.l)
        dl=np.radians(self.dl)
        
        x=d*np.cos(l-np.pi/2) #in kpc
        y=d*np.sin(l-np.pi/2)+R0
        
        dx=np.sqrt((np.cos(l-np.pi/2)*dd)**2+(-d*np.sin(l-np.pi/2)*dl)**2)
        dy=np.sqrt((np.sin(l-np.pi/2)*dd)**2+(d*np.cos(l-np.pi/2)*dl)**2)
        
        return x,y,dx,dy
                            
    def save_fit_rotation_curve(self):
        """This function creates the ouput txt filename 
        to create and save the graph as a png.
        It also create a txt file with the values of fit
        and Vmax and R
        raises write_txt"""
        # returns name without extension
        name=os.path.splitext(self.name)[0]  
        #creates the name of the output file depending on the input one
        output_name="res_"+name 
        
        main_dir = os.getcwd() 
        results_julia_dir = os.path.join(main_dir,"results_julia")
        if not os.path.exists(results_julia_dir):
            os.mkdir(results_julia_dir)
        output_filename = os.path.join(results_julia_dir,output_name)
        
        plt.savefig(output_filename+".png")
        self.write_txt(output_filename)
    
    def save_fit(self):
        """This function creates the ouput txt filename 
        to create and save the graph as a png"""
        name=os.path.splitext(self.name)[0] 
        output_name="res_"+name 
        
        main_dir = os.getcwd()
        results_julia_dir = os.path.join(main_dir,"results_julia")
        if not os.path.exists(results_julia_dir):
            os.mkdir(results_julia_dir)
        output_filename = os.path.join(results_julia_dir,output_name)
        
        plt.savefig(output_filename+".png")
        
    def write_txt(self,output_filename):
        """This function should write the relative velocity 
        and some other important values for the rotation curve
        to the output file.
        output_filename: string"""
        text=open(output_filename+".txt", "w")
        parameters_name=inspect.getargspec(self.fit_model).args
        matrix=self.pcov
        error=np.sqrt(matrix.diagonal())
        
        print("Fit model:", self.fit_model, '\n', file=text)
        print("Glon: %.2f $\pm$ %.2f" %(self.l,self.dl), file=text)
        for ii in range(len(self.popt)):
            print(parameters_name[ii+1], "%.2f $\pm$ %.2f" %(self.popt[ii], error[ii]), file=text)
        v,dv=self.get_velocity()
        for ii in range(len(v)):
            print("vr %d %.2f $\pm$ %.2f" %(ii, v[ii],dv[ii]), file=text)
        print('\nV_max %.0f $\pm$ %.0f' %(self.Vabs, self.dVabs), file=text)
        print('R %.4f $\pm$ %.4f' %(self.R, self.dR), file=text)
        
        text.close()