import numpy as np
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
import numpy as np
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

    return combined_array, combined_peak_locations, combined_peak_frequencies

peak_info = file_import('test.csv')

# Decodes datastreams from file_import function
combined_array = peak_info[0]
combined_peak_locations = peak_info[1]
combined_peak_frequencies = peak_info[2]

def window(block_peaks, block_size,combined_peak_locations,combined_array):

    freq = combined_array[:,0].tolist()
    abs = combined_array[:,1].tolist()
    combined_array = combined_array.tolist()
    peak_blocks_index = []

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

    return peak_blocks_freq_points, peak_blocks_abs_points

window = window(3,2,combined_peak_locations,combined_array)
freq = window[0]
print('Freq:',freq)
print(freq[0])
abs = window[1]
print('Abs:',abs)
print(abs[0])

# Create a normal distribution for data
def norm(x, mean, sd):
  norm = []
  for i in range(len(x)):
    norm += [1.0/(sd*np.sqrt(2*np.pi))*np.exp(-(x[i] - mean)**2/(2*sd**2))]
  return np.array(norm)

normalized_block = []

for i in abs:
    mean = statistics.mean(i)
    stdev = statistics.stdev(i)
    print('mean',mean,'stddev',stdev)
    normal_dist = norm(i,mean,stdev)
    normalized_block.append(normal_dist)

fig, axs = plt.subplots(2)
axs[0].plot(freq[0],abs[0],color = 'red')
axs[1].plot(freq[0],normalized_block[0], color = 'blue')

plt.show()

# for norm_abs in normalized_block:
#     for freq in freq:
#         plt.plot(norm_abs,freq,color = 'black')

print(normalized_block)

mean1, mean2 = 0, -2
std1, std2 = 0.5, 1

x = np.linspace(-20, 20, 500)
y_real = norm(x, mean1, std1) + norm(x, mean2, std2)

######################################
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
  err = y - y_fit
  return err

plsq = leastsq(res, p, args = (y_real, x))
print(plsq)

y_est = norm(x, plsq[0][0], plsq[0][2]) + norm(x, plsq[0][0] + plsq[0][1], plsq[0][3])

plt.plot(x, y_real, label='Real Data')
#plt.plot(x, y_init, 'r.', label='Starting Guess')
plt.plot(x,norm(x, plsq[0][0], plsq[0][2]), label = 'Peak 1') # Peak 1
plt.plot(x,norm(x, plsq[0][0] + plsq[0][1], plsq[0][3]), label = 'Peak 2') # Peak 2
plt.plot(x, y_est, 'g.', label='Fitted')
plt.legend()
plt.show()