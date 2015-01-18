import numpy as np
import sys
import matplotlib.pyplot as plt
import os
from scipy import stats
from lmfit.models import GaussianModel

# This is the place where global variable is defined:
infile = 'ZHIL1669.arw' # Input file name for Empower Raw data
bl_start = 1.
bl_end = 15. # Stand and End point for baseline correction
start = 17. # user-input term for the start of peak recognization
end = 33.   # user-input term for the end of peak recognization
amp_threshold = 5 # Peak threshold to identify as peak
# Function definition

def add_peak_info(p_loc, p_amp, p_num):
    '''
    This is the function to add the peak information to the p_list dictionary object
    Usage:
    peak = add_peak_info(p_loc, p_amp, p_left, p_right, index, p_num)
    '''
    peak_name = 'g' + str(p_num) + '_'
    peak_location = peak_name + 'center'
    peak_amp = peak_name + 'amplitude'
    peak = dict(zip([peak_name,peak_location,peak_amp],[p_num,p_loc,p_amp]))
    return peak
def gaussian (x,A,mu,sigma):
    '''
    This is the function to build the matrix to plot the individual peak,
    Usage: y = gaussian (x, A, mu, sigma)
    where x is the data (ndarray object); A, mu and sigma is amplitude, mean and standard deviation in 
    f(x; A, \\mu, \\sigma) = \\frac{A}{\\sigma\\sqrt{2\\pi}} e^{[{-{(x-\\mu)^2}/{{2\\sigma}^2}}]},
    '''
    res = A/(sigma*np.sqrt(2*np.pi))*np.exp(-(x-mu)**2/(2*sigma**2))
    return res

def calibration(time):
    '''
    This is the fuction to calculate the Log Mw from the Elution time
    Usage:
    Mp_Log = calibration(time)
    '''
    a = 2.27E1
    b = -1.77E0
    c = 6.21E-2
    d = -8.15E-4
    res = a + b * time + c * time**2 + d * time**3
    return res

def print_peak_info(peaks, num):
    '''
    This function is to pretify the peak information print out function
    Currently it is only used for terminal based function.
    The newer GUI interface would need a better layout
    Usage:
    print_peak_info(peaks, num)
    while peaks is the list of peaks. Each peak is a
    dictionary with peakM, Mn, Mw, PDI as keys
    '''
    for i in range(1,num):
        name = 'g' + str(i) + '_'
        print('Here is the information of Peak {0!s}'.format(i))
        print('Peak area is {0:.2f}'.format(peaks[i-1][name+'amplitude']))
        print('Peak Mn is {0:,.0f}'.format(peaks[i-1][name+'Mn']))
        print('Peak Mw is {0:,.0f}'.format(peaks[i-1][name+'Mw']))
        print('Peak Mp is {0:,.0f}'.format(peaks[i-1][name+'peakM']))
        print('Peak PDI is {0:.2f}'.format(peaks[i-1][name+'PDI']))
        print('{0:-<15}'.format('-'))

# os.chdir('C:\\Users\\zhil\\Documents\\ZQ_Programs\\Coding')
if (sys.platform != 'linux'):
    os.chdir('C:\\Users\\zhil\\Documents\\Coding')
else:
    os.chdir('/home/digitalpig/Coding/GPC_Fit')
os.getcwd()

# dat = np.loadtxt(open("ZHIL1669.arw","rb"),delimiter="\t",skiprows=1)
# Empower ARW file does not follow the line end convention in either Linux or Windows way, need to convert


outfile = infile + '-fixed'
with open(infile, 'rb') as f, open(outfile, 'wb') as g:
    content = f.read()
    g.write(content.replace(b'\r', b'\r\n'))
dat = np.loadtxt(open(outfile,"rb"),delimiter="\t")
x = dat[:, 0]
y = dat[:, 1]

# I need to do the linear regression first to get rid of baseline drifting.


ibl_start = np.where(x==bl_start)[0]
ibl_end = np.where(x==bl_end)[0]
slope, intercept, r_value, p_value, std_err = stats.linregress(x[ibl_start:ibl_end],y[ibl_start:ibl_end])
print('the R-sqaure of baseline fitting is: {0:.4}'.format(r_value))
y_drifting = slope * x + intercept
y = (y - y_drifting)

# Get the code from https://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-in-the-right-way

# The problem now is how to find the initial parameters for peak fit.
y1 = np.diff(y, 1)

# The next step is to find out the zero points for 1st and 2nd degree of differentiation...

x_start = np.where(x==start)[0]
x_end = np.where(x==end)[0]
peak_location = []
amp = []

previous = start
for i in range(x_start-1,x_end-1):
    if (y1[i] > sys.float_info.epsilon) and (y1[i+1] <sys.float_info.epsilon):  
    # Found that it is better not to smooth
        if (y[i] > amp_threshold):
            peak_location.append(x[i])
            amp.append(y[i])
            previous = x[i]
        continue
# Then, I need to clean the peak location data
peaks = []

num = 1
for i in range(len(peak_location)):
    peaks.append(add_peak_info(peak_location[i], amp[i],num))
    num += 1

# TODO: Here the A parameter considering the std. Manually input now but need 
# to be changed.


# Here we want to have the program to build the model by peak information.

# First, We will always have Peak 1. Also, need to do some initiation for peak
# 1.

gauss1  = GaussianModel(prefix='g1_')
pars = gauss1.make_params()

for i in range(1, num):
    name = 'g' + str(i) + '_'
    model = 'gauss' + str(i)
    x2 = peaks[i-1][name+'center']
    sigma = float(input("The sigma of peak {0} is".format(name)))
    y2 = peaks[i-1][name+'amplitude']/(sigma*np.sqrt(2*np.pi))
    plt.plot((x2,x2),(0,y2*sigma*np.sqrt(2*np.pi)),'-')
    if i == 1:
        pars[name+'center'].set(peaks[i-1][name+'center'])
        pars[name+'sigma'].set(sigma)
        pars[name+'amplitude'].set(peaks[i-1][name+'amplitude'])
        mod = gauss1
    else:
        model = GaussianModel(prefix=name)
        pars.update(model.make_params())
        pars[name+'center'].set(peaks[i-1][name+'center'])
        pars[name+'sigma'].set(sigma)
        pars[name+'amplitude'].set(peaks[i-1][name+'amplitude'])
        mod = mod + model

plt.plot(x[x_start:x_end],y[x_start:x_end])
plt.show()
print(peaks) 
   
out = mod.fit(y[x_start:x_end], pars, x=x[x_start:x_end])
print(out.fit_report(min_correl=0.5))

#comps = out.eval_components(x=x[x_start:x_end])
fitted_val = out.best_values

plt.plot(x[x_start:x_end],y[x_start:x_end])
plt.plot(x[x_start:x_end], out.best_fit, 'r-')
for i in range(1,num):
    name = 'g' + str(i) + '_'
    yp = gaussian(x[x_start:x_end],fitted_val[name+'amplitude'],fitted_val[name+'center'],fitted_val[name+'sigma'])
    plt.plot(x[x_start:x_end],yp,'--')

plt.show()

for i in range(1, num):
    name = 'g' + str(i) + '_'
    peaks[i-1][name+'center'] = fitted_val[name+'center']
    peaks[i-1][name+'area'] = fitted_val[name+'amplitude']
    peaks[i-1][name+'peakM'] = 10**(calibration(fitted_val[name+'center']))
    M_sigma = (calibration((peaks[i-1][name+'center']+fitted_val[name+'sigma']))
            - calibration((peaks[i-1][name+'center']-fitted_val[name+'sigma'])))
    M_sigma /= 2
    print(M_sigma)
    peaks[i-1][name+'Mn'] = peaks[i-1][name+'peakM']*np.exp(-(M_sigma**2)/2)
    peaks[i-1][name+'Mw'] = peaks[i-1][name+'peakM']*np.exp(+(M_sigma**2)/2)
    peaks[i-1][name+'PDI'] = np.exp(+(M_sigma**2))
    
print_peak_info(peaks,num)






