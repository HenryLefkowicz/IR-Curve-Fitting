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

print('Freq:',freq)
print('Abs:',abs)
print('Peak F:A Pairs', pbpaf)

concat_freq = [j for i in window[0] for j in i]
concat_abs = [j for i in window[1] for j in i]

#plt.plot(concat_freq,concat_abs)

# Create a normal distribution for data
def norm(x, mean, sd):
  norm = []
  for i in range(len(x)):
    norm += [1.0/(sd*np.sqrt(2*np.pi))*np.exp(-(x[i] - mean)**2/(2*sd**2))]
  return np.array(norm)

normalized_block = []

mean = statistics.mean(concat_freq)
stdev = statistics.stdev(concat_freq)

print('mean',mean,'stddev',stdev)
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

mean1, mean2 = 0, -2
std1, std2 = 0.5, 1

x = np.linspace(-20, 20, 500)
y_real = norm(x, mean1, std1) + norm(x, mean2, std2)


print(y_real)

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
  # creates bivariate normal distribution by using the norm function which creates a normal curve
  # based on the x values and the initial guesses for the mean and stddev
  err = y - y_fit
  return err

plsq = leastsq(res, p, args = (y_real, x))
plsq1 = leastsq(res,p,args = ())
print(plsq)

y_est = norm(x, plsq[0][0], plsq[0][2]) + norm(x, plsq[0][0] + plsq[0][1], plsq[0][3])

plt.plot(x, y_real, label='Real Data')
# #plt.plot(x, y_init, 'r.', label='Starting Guess')
# plt.plot(x,norm(x, plsq[0][0], plsq[0][2]), label = 'Peak 1') # Peak 1
# plt.plot(x,norm(x, plsq[0][0] + plsq[0][1], plsq[0][3]), label = 'Peak 2') # Peak 2
# plt.plot(x, y_est, 'g.', label='Fitted')
plt.legend()
plt.show()