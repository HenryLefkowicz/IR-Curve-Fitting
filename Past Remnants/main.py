import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from statistics import stdev
from collections import OrderedDict
import math
import pandas as pd


# TODO: Update documentation (7/03/2023)
# Why do I have a sinking feeling that this isn't going to fit different IRs well...
ir_dict = {}
ir_dict = OrderedDict()

# Used to knock the first X items off CSV list
# Changes due to IR software CSV orientation
n = 2
values = []
ir_int = []
ir_fr = []

# Open file as CSV
with open ('../test.csv', 'r') as csv_file:
    df_data = pd.read_csv(csv_file, sep=",", header=2)

# Required function for np.vectorize. Does the math to convert transmitence into absorbance
def logcalc(x):
    return -(math.log10(float(x) / 100))

# Converts data into array N2-Dim Array
# Col 1: Frequency, Col 2: Transmittence
array_df = df_data.to_numpy(dtype = 'float32')
frequency = array_df[:,0]
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


plt.scatter(frequency, absorbance, s=1, color='black')
plt.show()



# # Find peaks from intensity data
# # find_peaks returns the INDICES of the peaks, *NOT* the peaks themselves
# # this has to be done in two sets because find_peaks() can't find both at once
# # TODO: Identify flat regions / inflection points
# low_peaks, low_ = find_peaks(-inten_arr)
# high_peaks, high_ = find_peaks(inten_arr)
#
# # converts arrays (find_peaks output) into lists for comprehension
# low_peaks_list = low_peaks.tolist()
# high_peaks_list = high_peaks.tolist()
# # combines and sorts the two lists to make one master one.
# combined_peak_list = low_peaks_list + high_peaks_list
# combined_peak_list.sort()
#
# # converts lists into dicts
#     low_peaks_dict = {}
#     high_peaks_dict = {}
#     combined_peaks_dict = {}
#
#     for i in low_peaks:
#         low_peaks_dict[pointlist[i][0]] = pointlist[i][1] - noise
#     for i in high_peaks:
#         high_peaks_dict[pointlist[i][0]] = pointlist[i][1] - noise
#     for i in combined_peak_list:
#         combined_peaks_dict[pointlist[i][0]] = pointlist[i][1] - noise
#
#     low_peaks_INVERT_dict = high_peaks_dict
#     high_peaks_INVERT_dict = low_peaks_dict
#
#     # print('LPD', low_peaks_dict)
#     # print('HPD', high_peaks_dict)
#     # print('CPD', combined_peaks_dict)

def peak_calc(low_peaks_dict, high_peaks_dict):
    '''

    :param low_peaks_dict: Dict with K:V of Frequency:Intensity for low peaks
    :param high_peaks_dict: Dict with K:V of Frequency:Intensity for high peaks
    :return: 4 Arrays, 2 for the frequencies of the high and low peaks
    and 2 for the intensities of the high and low peaks
    '''

    # Pull keys and values into lists and converts them into arrays
    # for later processing
    # It's not actually that deep. Idk if it ~needs~ a function persay
    low_peaks_frequency = list(low_peaks_dict.keys())
    low_peaks_intensity = list(low_peaks_dict.values())

    high_peaks_frequency = list(high_peaks_dict.keys())
    high_peaks_intensity = list(high_peaks_dict.values())

    low_peaks_intensity_arr = np.array(low_peaks_intensity)
    low_peaks_frequency_arr = np.array(low_peaks_frequency)
    high_peaks_intensity_arr = np.array(high_peaks_intensity)
    high_peaks_frequency_arr = np.array(high_peaks_frequency)


    return (low_peaks_frequency_arr,high_peaks_frequency_arr,
            low_peaks_intensity_arr,high_peaks_intensity_arr)
def peak_break(peak_list_local,freq_arr):

    # kicks the first value of the list off because it just tells you
    # if the first value is a high or low point. It's important, but not here.
    peak_list_local.pop(0)

    local_peak_points = []
    sub = []

    print(peak_list_local)

    # goes through each peak range and then creates a list of all the integers
    # that are in that peak range, not just the bounds
    # and stores them in a sublist bc you need it for gaussian calculations and shit

    for i in range(len(peak_list_local)):
        for j in range(len(freq_arr)):
            if freq_arr[j] >= peak_list_local[i][0] and freq_arr[j] <= peak_list_local[i][-1]:
                sub.append(freq_arr[j])
        local_peak_points.append(sub)
        sub = []

    # Adds points to gaussian range calculation
    # TODO: Fix weird double plotting thing
    # Problem: It adds the points in a weird fucking order, I think it's a list slicing thing

    # Adding points to the range of the gaussian calculation makes it better
    # so this guy goes and adds extra points to each side of each sublist from the
    # sublists on either side to increase the range.

    # Yes, I recognize that from a time complexity standpoint this is a horrible
    # way to do this

    # Creates a reference version of local_peak_points that never changes size
    # so you can look at it to figure out what points you need to add
    # this way, local_peak_points is the only one that gets changed
    uc_local_peak_points = local_peak_points

    # for i in range(len(local_peak_points)):
    #     try:
    #         front = uc_local_peak_points[i+1][0:8]
    #         for j in range(len(front)):
    #             local_peak_points[i].append(front[j])
    #         if i != 0:
    #             back = uc_local_peak_points[i-1][-8:]
    #             #back.reverse()
    #             for k in range(len(back)):
    #                 local_peak_points[i].insert(0,back[k])
    #         # if you try to do this for an index of zero (the first sublist), it'll error out,
    #         # so the first peak just won't be as accurate. No way around that
    #         # without making up data.
    #         elif i == 0:
    #             pass
    #     except IndexError:
    #         break

    return local_peak_points
def min_max_min_id(peak_list_global):

    # Decode peak_list into individual data streams
    low_frequency_peaks_arr = peak_list_global[0]
    high_frequency_peaks_arr = peak_list_global[1]
    low_intensity_peaks_arr = peak_list_global[2]
    high_intensity_peaks_arr = peak_list_global[3]

    # It's important to know which one of the peaks comes first
    # so, this figures that out.

    if low_frequency_peaks_arr[0] < high_frequency_peaks_arr[0]:
        LIPf_Lowest = True
        HIPf_Lowest = False
    elif low_frequency_peaks_arr[0] > high_frequency_peaks_arr[0]:
        LIPf_Lowest = False
        HIPf_Lowest = True


# TODO: Fix the LIPf_Lowest so it puts the frequencies in the correct order
# TODO: Fix gaussian shift so it works for instances where the low peak comes first
    if LIPf_Lowest == True:
        holder = ['LIPf_Lowest']
        sublist = []
        for i in range(len(low_frequency_peaks_arr)):
            try:
                sublist.append(low_frequency_peaks_arr[i])
                sublist.append(high_frequency_peaks_arr[i])
                sublist.append(low_frequency_peaks_arr[i+1])
                holder.append(sublist)
                sublist = []
            except IndexError:
                break

    elif HIPf_Lowest == True:
        holder = ['HIPf_Lowest']
        sublist = []
        #sublist.append(freq_arr[0])
        # you don't want your calculated peak to start with the peak as it's first point,
        # so you have to add another point before it, but you have to convert it to a list first
        high_frequency_peaks_arr = high_frequency_peaks_arr.tolist()
        high_frequency_peaks_arr.insert(0,freq_arr[0])
        for i in range(len(low_frequency_peaks_arr)):
            try:
                # for literally the first term only, you want that very first, newly added term
                # then the first peak ([i+1]), then the peak from the other list, in this case
                # it's the low peaks bc this if for high peaks showing up first.
                if i == 0:
                    sublist.append(high_frequency_peaks_arr[i])
                    sublist.append(high_frequency_peaks_arr[i+1])
                    sublist.append(low_frequency_peaks_arr[i])
                    holder.append(sublist)
                    sublist = []
                # for all other times, you want to essentially pretend that first point didn't
                # exist, so you have to look backwards one ([i-1]) as your starting position for the peaks
                # you didn't add anything to and for the peak you added the extra value for,
                # you want to look one ahead of that so you can ignore it.
                elif i != 0:
                    sublist.append(low_frequency_peaks_arr[i-1])
                    sublist.append(high_frequency_peaks_arr[i+1])
                    sublist.append(low_frequency_peaks_arr[i])
                    holder.append(sublist)
                    sublist = []
            except IndexError:
                break

    return holder

def peak_break(peak_list_local,freq_arr):

    # kicks the first value of the list off because it just tells you
    # if the first value is a high or low point. It's important, but not here.
    peak_list_local.pop(0)

    local_peak_points = []
    sub = []

    print(peak_list_local)

    # goes through each peak range and then creates a list of all the integers
    # that are in that peak range, not just the bounds
    # and stores them in a sublist bc you need it for gaussian calculations and shit

    for i in range(len(peak_list_local)):
        for j in range(len(freq_arr)):
            if freq_arr[j] >= peak_list_local[i][0] and freq_arr[j] <= peak_list_local[i][-1]:
                sub.append(freq_arr[j])
        local_peak_points.append(sub)
        sub = []

    # Adds points to gaussian range calculation
    # TODO: Fix weird double plotting thing
    # Problem: It adds the points in a weird fucking order, I think it's a list slicing thing

    # Adding points to the range of the gaussian calculation makes it better
    # so this guy goes and adds extra points to each side of each sublist from the
    # sublists on either side to increase the range.

    # Yes, I recognize that from a time complexity standpoint this is a horrible
    # way to do this

    # Creates a reference version of local_peak_points that never changes size
    # so you can look at it to figure out what points you need to add
    # this way, local_peak_points is the only one that gets changed
    uc_local_peak_points = local_peak_points

    # for i in range(len(local_peak_points)):
    #     try:
    #         front = uc_local_peak_points[i+1][0:8]
    #         for j in range(len(front)):
    #             local_peak_points[i].append(front[j])
    #         if i != 0:
    #             back = uc_local_peak_points[i-1][-8:]
    #             #back.reverse()
    #             for k in range(len(back)):
    #                 local_peak_points[i].insert(0,back[k])
    #         # if you try to do this for an index of zero (the first sublist), it'll error out,
    #         # so the first peak just won't be as accurate. No way around that
    #         # without making up data.
    #         elif i == 0:
    #             pass
    #     except IndexError:
    #         break

    return local_peak_points

def gaussian_calc(full_peak_list, peak_list_local_bound, peak_list_global,intensity_ranges,flip_value):

    flip_value = flip_value

    '''
    :param full_peak_list: Entire spectra separated into their individual points. Each peak is
    represented in a sublist which contains all its frequency points.
    :param peak_list_local_bound: the highest point of the peak and it's two external bounds
    each is contained in a sublist
    :param peak_list_global: an dual array that represents each peak intensity and the frequency
    that peak is found at. Contains peaks for both the highest peaks and lowest peaks
    :param intensity_ranges: like full_peak_list, but for the intensities
    :param flip_value: tells the function if it's doing the inverse or non-inverse calculation
    :return: List of each calculated gaussian value for each subpeak
    '''

    # Required:
    # Stand. Dev for local peak range
    # Peak Height for local peak range
    # All other points

    gaussian_peaks = []

    if peak_list_local_bound[0] == 'LIPf_Lowest':
        e = 2.718281
        for subpeak in range(len(full_peak_list)):
            # Determines which subpeak we're iterating over
            subpeak_gaussian = []
            # this is just going to give you the original intensity value
            # that you're going to do math with (I think, maybe something different)
            a = peak_list_global[3][subpeak]
            # standard deviation of the subpeak req. for gaussian calc
            mu = stdev(full_peak_list[subpeak])
            # subpeak frequency that math is being done on
            # TODO: Change this so it works for LIPf_lowest
            b = peak_list_local_bound[subpeak+1][1]
            for point in full_peak_list[subpeak]:
                # does all the math bits
                over = (point-b)**2
                under = 2*(mu**2)
                exp = -(over/under)
                calc = a * e**(exp)
                # multiplies calculated value by flip_value depending on whether or not
                # it's inverted
                subpeak_gaussian.append(calc * flip_value)
            # Appends sublist to master list
            gaussian_peaks.append(subpeak_gaussian)

    # Yeah, um, I don't want to write the documentation twice but it's the same
    elif peak_list_local_bound[0] == 'HIPf_Lowest':
        e = 2.718281
        for subpeak in range(len(full_peak_list)):
            subpeak_gaussian = []
            a = peak_list_global[3][subpeak]
            mu = stdev(full_peak_list[subpeak])
            b = peak_list_local_bound[subpeak+1][1]
            for point in freq:
                over = (point-b)**2
                under = 2*(mu**2)
                exp = -(over/under)
                calc = a * e**(exp)
                #print(a)
                subpeak_gaussian.append(calc * flip_value)
            #print(len(subpeak_gaussian))
            gaussian_peaks.append(subpeak_gaussian)

# This adds all the points together, but you might need to do a convolution instead
# TODO: Check to make sure this is the right way to do this
    combined_gaussian = []
    for i in range(len(freq_arr)):
        holder = [item[i] for item in gaussian_peaks]
        print(len(holder))
        combined_gaussian.append(sum(holder))
        #print(holder)
    print(len(combined_gaussian))
    return combined_gaussian

def intensity_peaks(local_peak_points):

    intensity_ranges = []

    for subpeak in local_peak_points:
        sbpk = []
        for point in subpeak:
            #print(point)
            sbpk.append(ir_dict[int(point)])
        intensity_ranges.append(sbpk)

    return intensity_ranges

def optimizer():

    '''

    1. Import freq_arr

    :return:
    '''


    pass



def peak_print(peak_list_global,full_peak_list,calculated_peaks,color,flip_value):

    low_frequency_peaks_arr = peak_list_global[0]
    high_frequency_peaks_arr = peak_list_global[1]
    low_intensity_peaks_arr = peak_list_global[2]
    high_intensity_peaks_arr = peak_list_global[3]

    gaussian_peaks = calculated_peaks

    # Determines if the function is printing the inverted or the
    # non-inverted curves

    print(freq_arr)

    if flip_value == 1:
        plt.plot(freq_arr,gaussian_peaks, color = color)


    plt.plot(freq_arr,inten_arr, color = 'grey')
    plt.scatter(low_frequency_peaks_arr, low_intensity_peaks_arr, s=8, color='red')
    plt.scatter(high_frequency_peaks_arr, high_intensity_peaks_arr, s=8, color='black')
    #plt.show()

def main():

    color1 = 'green'

    peak_list_global = peak_calc(low_peaks_dict, high_peaks_dict)
    peak_list_local_bound = min_max_min_id(peak_list_global)
    full_peak_list = peak_break(peak_list_local_bound,freq_arr)
    intensity_ranges = intensity_peaks(full_peak_list)
    gaussian_peaks = gaussian_calc(full_peak_list, min_max_min_id(peak_list_global), peak_list_global,intensity_ranges,1)
    print(gaussian_peaks)
    peak_print(peak_list_global,full_peak_list,gaussian_peaks,color1,1)


    plt.show()

if __name__ == '__main__':
    main()