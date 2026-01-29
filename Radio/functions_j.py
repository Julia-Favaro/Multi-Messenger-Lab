import numpy as np

#preprocessing function
def moving_average(a):
    """A moving average is a way to smooth out data by calculating the average of a set of values over a specific window."""
    moving_average=np.array([]) 
    window_size=2
    for i in range(0,len(a)-1,window_size):
        this_window=a[i:i+window_size]
        window_average=np.mean(this_window)
        moving_average=np.append(moving_average,window_average)
    return moving_average

#postprocessing uncertainties (method that has been scarted)
def drel(a, mean):
    """Calculating standard deviation of points in a window around the mean. 
    Window size is defined by how many points you need to pass from the mean value to half of it """
    idx = (np.abs(a - mean)).argmin() #index of closest element to mean
    idx_half=(np.abs(a - mean/2)).argmin() #index of closest element to a third of the mean
    window_size=np.abs(idx-idx_half)
    window=a[idx-window_size: idx+window_size]
    error=np.std(window)
    return error

#--------------------------------------------------------------------------------------------------------------------------
#fitting functions
def model_gauss(x,c,A,x0,sigma):
    return c+A*np.exp(-np.power((x-x0),2.)/(2.*sigma**2))

def multi_gauss2(x,c0,c1,m1,s1,c2,m2,s2):
    gauss_1=c1*np.exp(-np.power((x-m1),2.)/(2.*s1**2))
    gauss_2=c2*np.exp(-np.power((x-m2),2.)/(2.*s2**2))
    return c0+gauss_1+gauss_2

def multi_gauss3(x,c0,c1,m1,s1,c2,m2,s2,c3,m3,s3):
    gauss_1=c1*np.exp(-np.power((x-m1),2.)/(2.*s1**2))
    gauss_2=c2*np.exp(-np.power((x-m2),2.)/(2.*s2**2))
    gauss_3=c3*np.exp(-np.power((x-m3),2.)/(2.*s3**2))
    return c0+gauss_1+gauss_2+gauss_3

def multi_gauss4(x,c0,c1,m1,s1,c2,m2,s2,c3,m3,s3,c4,m4,s4):
    gauss_1=c1*np.exp(-np.power((x-m1),2.)/(2.*s1**2))
    gauss_2=c2*np.exp(-np.power((x-m2),2.)/(2.*s2**2))
    gauss_3=c3*np.exp(-np.power((x-m3),2.)/(2.*s3**2))
    gauss_4=c4*np.exp(-np.power((x-m4),2.)/(2.*s4**2))
    return c0+gauss_1+gauss_2+gauss_3+gauss_4

