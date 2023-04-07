import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import math
import pandas as pd
from scipy.signal import find_peaks
import statistics

# Open file as CSV
with open ('test.csv', 'r') as csv_file:
    df_data = pd.read_csv(csv_file, sep=",", header=2)

# Required function for np.vectorize. Does the math to convert transmitence into absorbance
def logcalc(x):
    return -(math.log10(float(x) / 100))

# TODO: Put this in a function so you don't look like a fucking idiot
# Converts data into array N2-Dim Array
# Col 1: Frequency, Col 2: Transmittence
array_df = df_data.to_numpy(dtype = 'float32')
frequency = array_df[:,0]
frequency = frequency[::-1]
# Iterates through Col 2 and converts Transmittence into Absorbance
absorbance_calc = np.vectorize(logcalc)
intensity = absorbance_calc(array_df[:,1])

# Takes minimum value of intensity to remove from all points
noise =  min(intensity)
# Required function for np.vectorize. Does the math to remove noise from data
def denoise(x,noise):
    return x - noise
# Iterates through Col 2 and converts absorbance values into denoise-d absorbance values
denoise_intensity = np.vectorize(denoise)
absorbance = denoise_intensity(intensity, noise)
absorbance = absorbance[::-1]
intensity = intensity[::-1]

combined_array = np.column_stack((frequency,absorbance))
high_peak_locations = find_peaks(absorbance)
high_peak_locations = high_peak_locations[0].tolist()
low_peak_locations = find_peaks(-absorbance)
low_peak_locations = low_peak_locations[0].tolist()

combined_peak_locations = sorted(high_peak_locations + low_peak_locations)


combined_peak_freq = []
combined_peak_abs = []

for i in combined_peak_locations:
    combined_peak_freq.append(frequency[i])
    combined_peak_abs.append(absorbance[i])

x_values = combined_array[:,0].tolist()
y_values = combined_array[:,1].tolist()

combined_peak_frequencies = []
for i in combined_peak_locations:
    combined_peak_frequencies.append(x_values[i])


def window(range_over,combined_peak_locations,combined_array,combined_peak_frequencies):

    '''
    :param range_over: How many peaks you want to display
    :param combined_peak_locations: List of all peaks in order
    :param combined_array: Array of all frequency and intensity data in 2 Column Matrix
    :param combined_peak_frequencies: List of the peak indices for use in creating window bounds
    :return:
    - sh_freq: Full list of the windowed frequencies
    - sh_ab: Full list of the windowed intensities
    - sh_peaks_freq: Windowed peak frequencies
    - sh_peaks_abs: Window peak intensities
    '''

    #Established limit value - this is how many peaks you'll display
    limit = combined_peak_locations[range_over-1]

    # Separate the datastreams into two arrays
    sh_freq = combined_array[:, 0]
    sh_abs = combined_array[:, 1]
    # Break the freq and abs lists into shortened versions based on user
    # defined limits
    sh_freq = sh_freq[:limit].tolist()
    sh_abs = sh_abs[:limit].tolist()

    # Memory expensive, I know
    cpf = combined_peak_frequencies

    # Diagnostic print statements
    ## print('cpf: ', cpf)
    ## print('sh_freq: ', sh_freq)
    ## print('sh_abs: ',sh_abs)
    ## print('c_p_f: ', len(cpf))

    f_sbpk_h = [] # Full Subpeak Holder

    # TODO: Check to make sure this works in all instances
    # This is, quite possibly, the worst way to do this
    for i in range(len(cpf)):
        temp = []
        try:
            for j in sh_freq:
                if cpf[i] <= j <= cpf[i+1]: # checks to see if j is within peak range under study
                    temp.append(j)
            if sum(temp) != 0: # doesn't append anything if the list is empty
                f_sbpk_h.append(temp)
            else:
                pass
        except IndexError:
            break

    # Uses the index to represent the points that need to be iterated over
    # for later use by the gaussian function
    # I'm almost certain there is a 100X better way to do this
    sh_peaks_freq = []
    sh_peaks_abs = []
    sh_peaks_index = combined_peak_locations[0:range_over]
    for i in sh_peaks_index:
        sh_peaks_freq.append(combined_array[:,0][i])
        sh_peaks_abs.append(combined_array[:, 1][i])

    return sh_freq,sh_abs, sh_peaks_freq, sh_peaks_abs

window = window(2,combined_peak_locations,combined_array,combined_peak_frequencies)

# Decodes datastreams from window function
xdata = window[0]
ydata = window[1]
freq_peaks_window = window[2]
abs_peaks_window = window[3]

# Define the Gaussian function
# TODO: Make work for multiple gaussian functions

def gaussian(x,amp, cntr, stdev):
    return amp * np.exp(-(x - cntr) ** 2 / (2 * stdev ** 2))

# Define initial guesses for Gaussian parameters
amp_guess = statistics.mean(ydata)
mean_guess = statistics.mean(xdata)
stdev_guess = statistics.stdev(xdata)
p0 = [amp_guess,mean_guess,stdev_guess]

# Fit the curve to the data using non-linear least squares
# I'm going to be honest, I barely understand this.
# TODO: Iterate over all subpeaks
popt, pcov = curve_fit(gaussian, xdata, ydata, p0=p0)

# TODO: Add a main function so people don't think you're fucking insane
x = xdata
y = gaussian(x, *popt)
plt.plot(xdata, ydata, 'bo', label='data')
plt.plot(x, y, 'r-', label='fit')
plt.legend()
plt.show()
