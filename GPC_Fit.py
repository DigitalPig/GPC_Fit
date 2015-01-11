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
amp_threshold = 20 # Peak threshold to identify as peak
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

for i in range(1, num):
    name = 'g' + str(i) + '_'
    x2 = peaks[i-1][name+'center']
    sigma = float(input("The sigma of peak {0} is".format(name)))
    y2 = peaks[i-1][name+'amplitude']/(sigma*np.sqrt(2*np.pi))
    plt.plot((x2,x2),(0,y2*sigma*np.sqrt(2*np.pi)),'-')

plt.plot(x[x_start:x_end],y[x_start:x_end])
plt.show()
print(peaks) 
   


gauss1  = GaussianModel(prefix='g1_')
pars = gauss1.make_params()

pars['g1_center'].set(24, min=22, max=26)
pars['g1_sigma'].set(3, min=1)
pars['g1_amplitude'].set(82, min=10)

gauss2  = GaussianModel(prefix='g2_')

pars.update(gauss2.make_params())

pars['g2_center'].set(26, min=24, max=28)
pars['g2_sigma'].set(3, min=1)
pars['g2_amplitude'].set(60, min=20)

gauss3  = GaussianModel(prefix='g3_')

pars.update(gauss3.make_params())
pars['g3_center'].set(27, min=25, max=29)
pars['g3_sigma'].set(0.5, min=0.1)
pars['g3_amplitude'].set(40, min=20)
mod = gauss1 + gauss2 + gauss3
out = mod.fit(y[x_start:x_end], pars, x=x[x_start:x_end])
print(out.fit_report(min_correl=0.5))

#comps = out.eval_components(x=x[x_start:x_end])
fitted_val = out.best_values
yp1 = gaussian(x[x_start:x_end],fitted_val['g1_amplitude'],fitted_val['g1_center'],fitted_val['g1_sigma'])
yp2 = gaussian(x[x_start:x_end],fitted_val['g2_amplitude'],fitted_val['g2_center'],fitted_val['g2_sigma'])
yp3 = gaussian(x[x_start:x_end],fitted_val['g3_amplitude'],fitted_val['g3_center'],fitted_val['g3_sigma'])



plt.plot(x[x_start:x_end],y[x_start:x_end])
plt.plot(x[x_start:x_end], out.best_fit, 'r-')
plt.plot(x[x_start:x_end],yp1,'--',x[x_start:x_end],yp2,'--',x[x_start:x_end],yp3,'--')
plt.show()

for i in range(1, num):
    name = 'g' + str(i) + '_'
    peaks[i-1][name+'center'] = fitted_val[name+'center']
    peaks[i-1][name+'area'] = fitted_val[name+'amplitude']
    peaks[i-1][name+'peakM'] = calibration(fitted_val[name+'center'])
    M_sigma = (calibration((peaks[i-1][name+'center']+fitted_val[name+'sigma']))
            - calibration((peaks[i-1][name+'center']-fitted_val[name+'sigma'])))
    M_sigma /= 2
    print(M_sigma)
    peaks[i-1][name+'Mn'] = peaks[i-1][name+'peakM']*np.exp(-(M_sigma**2)/2)
    peaks[i-1][name+'Mw'] = peaks[i-1][name+'peakM']*np.exp(+(M_sigma**2)/2)
    peaks[i-1][name+'PDI'] = np.exp(+(M_sigma**2))
    
print(peaks)    






