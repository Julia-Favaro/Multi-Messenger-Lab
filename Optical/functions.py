import numpy as np

def model_gauss(x,A,x0,sigma):
    return A*np.exp(-np.power((x-x0),2.)/(2.*sigma**2))

def line(x,q,m):
    return q+m*x

def gain(B1,err_B1,B2,err_B2,sigmaB12,err_B12,F1,err_F1,F2,err_F2,sigmaF12,err_F12):
    gain=((F1+F2)-(B1+B2))/(sigmaF12**2-sigmaB12**2)
    frac1=(1/(sigmaF12**2-sigmaB12**2))
    derF12=2*sigmaF12*(F1+F2-B1-B2)/(sigmaF12**2-sigmaB12**2)**2
    derB12=2*sigmaF12*(F1+F2-B1-B2)/(sigmaF12**2-sigmaB12**2)**2
    err_gain=np.sqrt(frac1**2*(err_B1**2+err_B2**2+err_F1**2+err_F2**2)+derF12**2*err_F12**2+derB12**2*err_B12**2)
    return gain,err_gain

def readout_noise(gain,err_gain,sigmaB12,err_B12):
    rnoise=gain*sigmaB12/np.sqrt(2)
    dergain=sigmaB12/np.sqrt(2)
    dersigma=gain/np.sqrt(2)
    err_rnoise=np.sqrt(dergain**2*err_gain**2+dersigma**2*err_B12**2)
    return rnoise,err_rnoise