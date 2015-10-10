# GPC_Fit
GPC Fit is a small python script to decovolute your GPC spectrums. It reads the
data file as well as the calibration curve exported from Waters HPLC and then
process them using non-linear least square curve fitting to fit the chromatogram
with Gaussian peaks.

# Prerequisite
The script is written in Python 3. It also depends on numpy, scipy, matplotlib
and lmfit. All packages except lmfit are most likely already available in your
python environment if you install it using pre-compiled package like WinPython,
Anaconda etc. Lmfit package can be found
[here](https://lmfit.github.io/lmfit-py/).

# Usage
Currently, GPC_Fit is just a command line script which only interacts with you
in a terminal. It does the following steps:

1. Baseline correction
2. Identify the peak location by looking at the 2nd derivative of the
   chromatogram. It then ask for peak width by inputing the sigma of the
   possible peak. Note: This is just the *initial* value of the parameter for
   the curve fitting. You can guess it from your chromatogram by measuring the
   peak's half peak width and then divide by 2.
3. A plot showing the chromatogram and identified peaks are listed for you to
   confirm. You can close it after you feel the recognization is acceptable.
4. It then perform the deconvolution based on the peak locations and peak width
   obtained previously.
5. The peak information will then be printed out in the terminal. Feel free to
   save the deconvoluted images if you want.
6. After you close the image, the molecular weight information for each peak is
   then printed on the terminal. Those are calculated based on the calibration
   curve that you provide.

# Configuration

At current stage, the program is crude that you have to alter the source code in
order to change some parameters. Here is the brief list:

## Input file name; Baseline correction starting and ending position

```python
infile = 'ZHIL1669.arw' # Input file name for Empower Raw data
bl_start = 1.
bl_end = 15. # Stand and End point for baseline correction
start = 17. # user-input term for the start of peak recognization
end = 33.   # user-input term for the end of peak recognization
amp_threshold = 5 # Peak threshold to identify as peak

```
The code itself is quite self-explaintory. "bl_start" and "bl_end" tell it the
beginning and ending of baseline correction. "amp_threshold" tells the program
that only peaks higher than it will be recognized as peak.

## Calibration curve

```python

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
```

Here is the place for you to put the parameter of your calibration curve. Here
the example is a 3rd order calibration curve.

# Next Plan

I got a lot of work to finish this tool. Here is something I am thinking of or
working at right now:

* A GUI frontend to wrap the process and parameter adjustment
* Introducing Langevin fitting
* Get compatiblize with Agilent GPC data
