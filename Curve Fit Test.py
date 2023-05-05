import numpy as np
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import math
import pandas as pd
from scipy.signal import find_peaks
import statistics

def file_import(ir_file):
# Open file as CSV
    with open (ir_file, 'r') as csv_file:
        df_data = pd.read_csv(csv_file, sep=",", header=2)

# Required function for np.vectorize. Does the math to convert transmitence into absorbance
    def logcalc(x):
        return -(math.log10(float(x) / 100))

    # Converts data into (number of datapoints) X 2 matrix
    # Col 1: Frequency, Col 2: Transmittance
    array_df = df_data.to_numpy(dtype = 'float32')
    frequency = array_df[:,0]
    frequency = frequency[::-1]

    # Iterates through Col 2 and converts transmittance into absorbance
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

    # Makes new combined array
    combined_array = np.column_stack((frequency,absorbance))

    # Finds high and low peaks
    high_peak_locations = find_peaks(absorbance)
    high_peak_locations = high_peak_locations[0].tolist()
    low_peak_locations = find_peaks(-absorbance)
    low_peak_locations = low_peak_locations[0].tolist()

    # Combines peak lists together
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

    return combined_array, combined_peak_locations, combined_peak_frequencies, combined_peak_abs

peak_info = file_import('test.csv')

# Decodes datastreams from file_import function
combined_array = peak_info[0]
combined_peak_locations = peak_info[1]
combined_peak_frequencies = peak_info[2]
combined_peak_abs = peak_info[3]
peak_freq_abs_array = np.vstack((combined_peak_frequencies,combined_peak_abs)).T

def window(block_peaks, block_size, window_start, window_end, combined_peak_locations,
           combined_array,peak_freq_abs_array):

    freq = combined_array[:,0].tolist()
    abs = combined_array[:,1].tolist()
    combined_array = combined_array.tolist()
    peak_blocks_index = []

    peak_freq_abs_dict = {}
    for j, k in peak_freq_abs_array:
        peak_freq_abs_dict[j] = k

    freq_abs_dict = {}
    for j,k in combined_array:
        freq_abs_dict[j] = k

    for i in range(block_peaks):
        peak_blocks_index.append(combined_peak_locations[0:block_size+1])
        combined_peak_locations = combined_peak_locations[block_size:]

    peak_blocks_freq_points = []
    for i in peak_blocks_index:
        peak_blocks_freq_points.append(freq[i[0]:i[-1]+1])

    peak_blocks_abs_points = []
    for block in peak_blocks_freq_points:
        temp = []
        for point in block:
            temp.append(freq_abs_dict[point])
        peak_blocks_abs_points.append(temp)

    pbfp = peak_blocks_freq_points[window_start:window_end+1]
    print(pbfp)
    pbap = peak_blocks_abs_points[window_start:window_end+1]


    pbpaf = [] # peak_blocks_peak_abs_freq
    for subpeak in pbfp:
        for point in subpeak:
            temp = []
            if point in peak_freq_abs_dict:
                temp += point, peak_freq_abs_dict[point]
                pbpaf.append(temp)

    return pbfp, pbap ,pbpaf

window = window(7,1,3,4,combined_peak_locations,combined_array,peak_freq_abs_array)
print(combined_peak_frequencies)
freq = window[0]
abs = window[1]
pbpaf = window[2]

# Takes the concatenated list creates by the window functions and
# creates a new list without all the sublist stuff. These data points can be used
# as the raw data for the optimization function.

concat_freq = [j for i in window[0] for j in i]
concat_abs = [j for i in window[1] for j in i]

mean = statistics.mean(concat_freq)
stdev = statistics.stdev(concat_freq)
print('Mean: ',mean)
print('Stdev: ',stdev)
print('Freq:',freq)
print('Abs:',abs)
print('Peak F:A Pairs', pbpaf)

# Create a normal distribution for data
def norm(x, mean, sd):
  norm = []
  for i in range(len(x)):
    norm += [1.0/(sd*np.sqrt(2*np.pi))*np.exp(-(x[i] - mean)**2/(2*sd**2))]
  return np.array(norm)

def plotter(concat_freq,concat_abs):

    normalized_block = []
    concat_freq,concat_abs = concat_freq,concat_abs
    mean = statistics.mean(concat_freq)
    stdev = statistics.stdev(concat_freq)

    print('Window Mean',mean,'Window Stddev',stdev)
    normal_dist = norm(concat_freq,mean,stdev)
    normalized_block.append(normal_dist)

    fig, axs = plt.subplots(2)
    axs[0].set_title('Original Data')
    axs[1].set_title('Normalized Data')

    axs[0].plot(concat_freq,concat_abs,color = 'red')
    axs[1].plot(concat_freq,normalized_block[0], color = 'blue')

    for i,(freq,abs) in enumerate(pbpaf):
        axs[0].scatter(freq,abs,color = 'orange')
        axs[1].scatter(freq,abs,color = 'orange')

    axs[0].scatter(pbpaf[0][0],pbpaf[0][1],color = 'orange', label = 'Peak Locations')
    axs[0].legend(loc = 'upper right')
    axs[1].scatter(pbpaf[0][0],pbpaf[0][1],color = 'orange', label = 'Peak Locations')
    axs[1].legend(loc = 'upper right')

    fig.tight_layout()
    plt.legend()
    plt.show()

# Non-Example Data:

# Decode + Create Means and Stdevs to use as initial guesses
actual_mean = mean
actual_mean1 = actual_mean + (actual_mean*0.01)
actual_stdev = stdev
actual_stdev1 = actual_stdev + (actual_stdev*0.01)

# Uses the actual datapoints as the raw inputs
actual_x = concat_freq
actual_y = concat_abs

# Shows the initial curve guesses
y_init = norm(actual_x, actual_mean, actual_stdev) + norm(actual_x, actual_mean1, actual_stdev1)

plt.show()
### Solving ###

# Initial Guesses
initial = [actual_mean,actual_mean1,actual_stdev,actual_stdev1]

def actual_res(initial, actual_y, actual_x):

    actual_y = np.array(actual_y)
    # Decode Data
    actual_mean,actual_mean1,actual_stdev,actual_stdev1 = initial

    # Create the fitting as two normal curves
    actual_y_fit = norm(actual_x, actual_mean, actual_stdev) + \
                 norm(actual_x, actual_mean1, actual_stdev1)

    # creates bivariate normal distribution by using the norm function which creates a normal curve
    # based on the x values and the initial guesses for the mean and stddev
    err = actual_y - actual_y_fit
    return err

plsq = leastsq(actual_res, initial, args = (actual_y, actual_x))
print(actual_res(initial, actual_y, actual_x))


print('plsq ',plsq)
actual_y_est = norm(actual_x, plsq[0][0], plsq[0][2]) + norm(actual_x, plsq[0][1], plsq[0][3])
plt.plot(actual_x, actual_y, label='Real Data')
plt.plot(actual_x, actual_y_est, 'r', label='Fitted')

plt.plot(actual_x,norm(actual_x, plsq[0][0], plsq[0][2]), label = 'Peak 1') # Peak 1
plt.plot(actual_x,norm(actual_x, plsq[0][1], plsq[0][3]), label = 'Peak 2') # Peak 2

# plt.plot(actual_x,norm(actual_x, actual_mean, actual_stdev),'black',label = 'Initial Guess Curve')
# plt.plot(actual_x,norm(actual_x, actual_mean1, actual_stdev1),'black')

# This just makes the graph larger
x = np.linspace(0, 1, 1000)
c =  [0] * 1000
plt.plot(x,c)

plt.legend()
plt.show()

### Recursive Shit ###

# Base Case: Error between actual_y_est and actual_y is less than .005% of actual_y

def average(actual_y, actual_y_est):
    for i in range(len(actual_y)):
        delta = 0
        delta += (actual_y_est[i] / actual_y[i])*100
        return delta/len(actual_y)

print('Average % Δ between fit and real data: ',average(actual_y,actual_y_est))


#########
# Example Data:
mean1, mean2 = 0, -2
std1, std2 = 0.5, 1
x = np.linspace(-20, 20, 500)
y_real = norm(x, mean1, std1) + norm(x, mean2, std2)

# Solving
m, dm, sd1, sd2 = [-2, 10, 1, 1]
p = [m, dm, sd1, sd2] # Initial guesses for leastsq
y_init = norm(x, m, sd1) + norm(x, m + dm, sd2) # For final comparison plot

# Pretty sure this creates the function that determines the difference between
# the calculated peak and the actual peak, which is then minimized as a function
# in scipy's leastsq

def res(p, y, x):
  m, dm, sd1, sd2 = p
  m1 = m
  m2 = m1 + dm
  y_fit = norm(x, m1, sd1) + norm(x, m2, sd2)
  # creates bivariate normal distribution by using the norm function which creates a normal curve
  # based on the x values and the initial guesses for the mean and stddev
  err = y - y_fit
  return err

plsq = leastsq(res, p, args = (y_real, x))

y_est = norm(x, plsq[0][0], plsq[0][2]) + norm(x, plsq[0][0] + plsq[0][1], plsq[0][3])

plt.plot(x, y_real, label='Real Data')
plt.plot(x, y_init, 'r.', label='Starting Guess')
plt.plot(x,norm(x, plsq[0][0], plsq[0][2]), label = 'Peak 1') # Peak 1
plt.plot(x,norm(x, plsq[0][0] + plsq[0][1], plsq[0][3]), label = 'Peak 2') # Peak 2
plt.plot(x, y_est, 'g.', label='Fitted')
plt.legend()
plt.show()