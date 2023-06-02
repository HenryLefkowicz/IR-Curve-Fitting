import math
import pandas as pd
from scipy.signal import find_peaks
import statistics
from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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

    '''

    :param block_peaks: How many blocks you want
    :param block_size: How many subpeaks in each block
    :param window_start: What subpeak the window starts at
    :param window_end: What subpeak the window ends at
    ## NOTE: The window ranges can override the block_peaks parameter.
    you can calculate N number of peaks, and it will give you that data points for all of those
    but you can still additionally truncate that with the window start:stop parameter
    :param combined_peak_locations: A list of where all the peaks are
    :param combined_array: All the freq:abs points for the entire data set
    :param peak_freq_abs_array: An array of the peak freq and abs values for the whole data set
    :return:Truncated range of datapoints based on parameters. Returns them in 3 separates lists.
    pbfp (frequency), pbap (absorbance), and pbpaf (peak frequency absorbance pairs)
    '''

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

    print('peak_blocks_freq_points'
          ,peak_blocks_freq_points)

    pbfp = peak_blocks_freq_points[window_start:window_end+1]
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

actual_x = concat_freq
actual_x_arr = np.array(actual_x)

actual_y = concat_abs
actual_y_arr = np.array(actual_y)


# Create a normal distribution for data
def normal(x, mean, sd):

    '''
    Makes normal curves used for just about everything else.

    :param x: Individual X Point
    :param mean: Mean of the data
    :param sd: Standard deviation
    :return: Returns a normal curve based on the inputted parameters.
    '''

    norm = []
    for i in range(len(x)):
        norm += [1.0/(sd*np.sqrt(2*np.pi))*np.exp(-(x[i] - mean)**2/(2*sd**2))]
    return np.array(norm)

def plotter(concat_freq,concat_abs):

    '''

    This guy just plots stuff. Useful mostly because it keeps thing tidy,
    but I'm not actually using it for anything right now.

    :param concat_freq: The actual frequency values in the range generated from the window function
    :param concat_abs: The actual abs values in the range generated from the window function
    :return:
    '''

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

def fun_rosenbrock(x):
    return np.array([10 * (x[1] - x[0]**2), (1 - x[0])])

p1 = [450,3]
p2 = [450,3,455,3]
p3 = [470,3,475,3,480,3]
p4 = [470,3,475,3,480,3,470,4]

def gaussian(x, mean, sd):
    return 1.0 / (sd * np.sqrt(2 * np.pi)) * np.exp(-(x - mean) ** 2 / (2 * sd ** 2))

def gaussian2(x,mean,sd,mean2,sd2):

     return 1.0 / (sd * np.sqrt(2 * np.pi)) * np.exp(-(x - mean) ** 2 / (2 * sd ** 2)) + \
        1.0 / (sd * np.sqrt(2 * np.pi)) * np.exp(-(x - mean2) ** 2 / (2 * sd2 ** 2))

def gaussian3(x, mean, sd, mean2, sd2, mean3, sd3):
    return 1.0 / (sd * np.sqrt(2 * np.pi)) * np.exp(-(x - mean) ** 2 / (2 * sd ** 2)) + \
        1.0 / (sd * np.sqrt(2 * np.pi)) * np.exp(-(x - mean2) ** 2 / (2 * sd2 ** 2)) + \
    1.0 / (sd * np.sqrt(2 * np.pi)) * np.exp(-(x - mean3) ** 2 / (2 * sd3 ** 2))

def gaussian4(x, mean, sd, mean2, sd2, mean3, sd3, mean4, sd4):
        return 1.0 / (sd * np.sqrt(2 * np.pi)) * np.exp(-(x - mean) ** 2 / (2 * sd ** 2)) + \
            1.0 / (sd * np.sqrt(2 * np.pi)) * np.exp(-(x - mean2) ** 2 / (2 * sd2 ** 2)) + \
        1.0 / (sd * np.sqrt(2 * np.pi)) * np.exp(-(x - mean3) ** 2 / (2 * sd3 ** 2)) + \
    1.0 / (sd * np.sqrt(2 * np.pi)) * np.exp(-(x - mean4) ** 2 / (2 * sd4 ** 2))



popt,pcov = curve_fit(gaussian2, actual_x_arr, actual_y_arr,p0=p2,bounds = (0,[500,10,500,10]))
print('Two Curves:', popt)

plt.figure('Two Curves')
plt.plot(actual_x_arr, actual_y_arr, '.b',label='Real Data')
plt.plot(actual_x_arr, gaussian2(actual_x_arr,popt[0],popt[1],popt[2],popt[3]))
plt.plot(actual_x_arr, gaussian(actual_x_arr,popt[0],popt[1]),label = 'Component 1',color = 'red')
plt.plot(actual_x_arr, gaussian(actual_x_arr,popt[2],popt[3]),label = 'Component 2',color = 'orange')

plt.legend()

popt,pcov = curve_fit(gaussian3, actual_x_arr, actual_y_arr,p0=p3,
        bounds = ((450, 0, 450, 0, 450, 0), (500, 10, 500, 10, 500, 10)))
print('Three Curves:', popt)
plt.figure('Three Curves')
plt.plot(actual_x_arr, actual_y_arr, '.b',label='Real Data')
plt.plot(actual_x_arr, gaussian3(actual_x_arr,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]))
plt.plot(actual_x_arr, gaussian(actual_x_arr,popt[0],popt[1]),label = 'Component 1', color = 'red' )
plt.plot(actual_x_arr, gaussian(actual_x_arr,popt[2],popt[3]),label = 'Component 2',color = 'orange')
plt.plot(actual_x_arr, gaussian(actual_x_arr,popt[4],popt[5]),label = 'Component 3',color = 'green')

popt,pcov = curve_fit(gaussian4, actual_x_arr, actual_y_arr,p0=p4,
        bounds = ((350, 0, 350, 0, 350, 0, 350, 0), (600, 10, 600, 10, 600, 10, 600, 10)), method = 'trf', maxfev = 10000)
print('Four Curves:', popt)
plt.figure('Four Curves')
plt.plot(actual_x_arr, actual_y_arr, '.b',label='Real Data')
plt.plot(actual_x_arr, gaussian4(actual_x_arr,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],popt[6],popt[7]))
plt.plot(actual_x_arr, gaussian(actual_x_arr,popt[0],popt[1]),label = 'Component 1', color = 'red' )
plt.plot(actual_x_arr, gaussian(actual_x_arr,popt[2],popt[3]),label = 'Component 2',color = 'orange')
plt.plot(actual_x_arr, gaussian(actual_x_arr,popt[4],popt[5]),label = 'Component 3',color = 'green')
plt.plot(actual_x_arr, gaussian(actual_x_arr,popt[6],popt[7]),label = 'Component 3',color = 'black')

plt.legend()
plt.show()



