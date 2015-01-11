import numpy as np
import sys
import matplotlib.pyplot as plt
import os
from scipy import fftpack, stats
from lmfit.models import GaussianModel

# This is the place where global variable is defined:
infile = 'ZHIL1669.arw' # Input file name for Empower Raw data
bl_start = 1.
bl_end = 15. # Stand and End point for baseline correction
start = 17. # user-input term for the start of peak recognization
end = 32.   # user-input term for the end of peak recognization

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
y1_fft = fftpack.rfft(y1)
y1_freq = fftpack.rfftfreq(len(y1),x[3]-x[2])

spectrum = y1_fft **2
cutoff_idx = spectrum < (spectrum.max()/5)
y1_fft_filtered = y1_fft.copy()
y1_fft_filtered[cutoff_idx] = 0
y1_filtered = fftpack.irfft(y1_fft_filtered)



y2 = np.diff(y, 2)
y2_fft = fftpack.rfft(y2)
y2_freq = fftpack.rfftfreq(len(y2),x[3]-x[2])

spectrum = y2_fft **2
cutoff_idx = spectrum < (spectrum.max()/5)
y2_fft_filtered = y2_fft.copy()
y2_fft_filtered[cutoff_idx] = 0
y2_filtered = fftpack.irfft(y2_fft_filtered)



# The next step is to find out the zero points for 1st and 2nd degree of differentiation...

x_start = np.where(x==start)[0]
x_end = np.where(x==end)[0]
peak_location = []
peak_inflection_left = []
peak_inflection_right = []
amp = []
amp_threshold = 20
del_loc = []
previous = start
for i in range(x_start,x_end-1):
    if (y1[i]*y1[i+1] < sys.float_info.epsilon):  # Found that it is better not to smooth
        if (x[i]-previous > 1.) and (y[i] > amp_threshold):
            peak_location.append(x[i])
            amp.append(y[i])
            previous = x[i]
        continue
    if (y2_filtered[i-1]*y2_filtered[i] < sys.float_info.epsilon):
        if (y2_filtered[i-1] > sys.float_info.epsilon):
            peak_inflection_left.append(x[i+2])
        else:
            peak_inflection_right.append(x[i+2])
# Then, I need to clean the peak location data


print('peak locations are: {0}'.format(peak_location))
print('peak left inflections are: {0}'.format(peak_inflection_left))
print('peak right inflections are: {0}'.format(peak_inflection_right))
print('peak amplitudes are: {0}'.format(amp))
peaks = []


# Iterate for a peak location/width/amplitude
def add_peak_info(p_loc, p_amp, p_left, p_right, p_num):
    '''
    This is the function to add the peak information to the p_list dictionary object
    Usage:
    peak = add_peak_info(p_loc, p_amp, p_left, p_right, index, p_num)
    '''
    peak_name = 'g' + str(p_num) + '_'
    peak_location = peak_name + 'center'
    peak_width = peak_name + 'sigma'
    peak_amp = peak_name + 'amplitude'
    peak = dict(zip([peak_name,peak_location,peak_width,peak_amp],[p_num,p_loc,p_right-p_left,p_amp]))
    return peak


num = 0

#TODO: Here is the last break should let the circle go to the next i circle (looking for next left inflection)
try:
    for i in range(len(peak_inflection_left)):
        left = peak_inflection_left[i]
        for j in range(len(peak_inflection_right)):
            if (peak_inflection_right[j] > left):
                right = peak_inflection_right[j]
                for k in range(len(peak_location)):
                    if ((left <= peak_location[k]) and (peak_location[k] <= right)):
                        num += 1
                        print('num is {0}'.format(num))
                        peaks.append(add_peak_info(peak_location[k],amp[k],left,right,num)) # call the function to add peak info to the dictionary
                        peak_location.pop(k)
                        amp.pop(k)
                        peak_inflection_left.pop(i)
                        peak_inflection_right.pop(j)
                        break
                peak_inflection_right.pop(j)
        peak_inflection_right.pop(j)
    peak_inflection_left.pop(i)
except IndexError:
    pass
                    

for i in range(1, num+1):
    name = 'g' + str(i) + '_'
    x2 = peaks[i-1][name+'center']
    y2 = peaks[i-1][name+'amplitude']
    displace = peaks[i-1][name+'sigma']
    print('name is {0}'.format(name))
    plt.plot((x2,x2),(0,y2),'-')
#    plt.plot((x2-displace,x2-displace),(0,y2),'--')
#    plt.plot((x2+displace,x2+displace),(0,y2),'--')

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



#out = mod.fit(y, pars, x=x)
out = mod.fit(y[x_start:x_end], pars, x=x[x_start:x_end])
print(out.fit_report(min_correl=0.5))

comps = out.eval_components(x=x[x_start:x_end])

plt.plot(x[x_start:x_end],y[x_start:x_end])
plt.plot(x[x_start:x_end], out.best_fit, 'r-')
plt.plot(x[x_start:x_end],comps['g1_'],'--',x[x_start:x_end],comps['g2_'],'--',x[x_start:x_end],comps['g3_'],'--')
plt.show()
print(out.best_values)





