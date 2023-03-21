import numpy as np
from scipy.signal import savgol_filter, find_peaks
from constants import *
import scipy.special as spspec


def strip_analysis(data_array):

    channels = len(data_array)  # matrix horizontal dimension: number of columns, each column corresponding to one acquisition channel (127 per axis for nextQ, 256 per axis for Q360, 127 to be divided in the two axes for Qplus)
    x_coord_raw = np.arange(0, channels)
    x_coord = np.arange(-(channels/2), channels/2) #centering to zero
    sample_studied = np.asarray(data_array)
    ym = savgol_filter(sample_studied, 7, 4)
    max_data = np.amax(ym)  # maximum element in the sample to study
    maxIndex = np.argmax(ym)
    
    #retrieving the fit parameters
    mean = meanPos(ym[maxIndex-1:maxIndex+2], x_coord[maxIndex-1:maxIndex+2])
    calc_fwhm = calc_sides_withBS(ym, x_coord_raw, 0.5)
    manual_fwhm = calc_fwhm[1] - calc_fwhm[0]
    # manual calculation of peak width at 90 % of the peak level
    calc_sides90 = calc_sides_withBS(ym, x_coord_raw, 0.9)
    manual_sides90 = calc_sides90[1] - calc_sides90[0]
    #arrays conversion to lists
    rawdata_list = sample_studied.tolist()
    rawcoord_list = x_coord.tolist()  # converting arrays to lists
    ym_list = ym.tolist()
    
    results = {
        "mean":{"value":float(mean), "unit":"mm"},
        "sigma": {"value": abs(float(manual_fwhm/CONVERSION_SIGMA)), "unit": "mm"},
        "fwhm": {"value": abs(float(manual_fwhm)), "unit": "mm"},
        "peak_width": {"value": abs(float(manual_sides90)), "unit": "mm"},
        "coordinates_raw": rawcoord_list,
        "raw_data":rawdata_list,
        "coordinates_fit": rawcoord_list,
        "fit_data":ym_list,
    }
    return results

def mlic_analysis(data_array, bortfeld):
    channels = len(data_array)  # number of acquisition channels
    x_coord_we = np.arange(1, channels+1)*TO_WE # coordinates of channels converted to water equivalent thickness
    # coordinates of fit converted to water equivalent thickness (1 bin = 1/4 channels -> increased accuracy for fit)
    x_coord_fit_we = np.arange(1, channels, 0.1)*TO_WE
    sample_studied = np.asarray(data_array)  #selecting the complete sample (all channels)
    raw_max = np.amax(sample_studied)# maximum value of raw data array
    low_lim = sample_studied[0] #if sample_studied[0] > sample_studied[-1] else sample_studied[-1] # plateau entrance value (with check for reversed data)
    peaks, properties = find_peaks(sample_studied, prominence=(low_lim, raw_max), width=1)# Find peaks in data
    index_max_width = 0
    if(len(properties['widths'])):
        index_max_width = np.argmax(properties['widths'])
    else:
        index_max_width = np.argmax(sample_studied)
    #correction of prominent spikes
    count = 0
    if peaks.size > 1:
        for el in peaks:
            if count != index_max_width and properties['widths'][count] < 1.5:
                sample_studied[peaks[count]] = (sample_studied[peaks[count]-1]+sample_studied[peaks[count]+1]) / 2
            count+=1
    # norm_factor = sample_studied[0]
    norm_theo = (np.amax(sample_studied) - sample_studied[0])/2 + sample_studied[0]
    index_norm, norm_factor = find_value(sample_studied, norm_theo)
    #NORMALIZING DATA CURVE
    new_sample = []
    for el in sample_studied:
        new_sample.append(el/norm_factor)
    sample_studied = np.asarray(new_sample)
    # max_after_correction = np.amax(sample_studied)
    # norm_level_raw = (max_after_correction - sample_studied[0])/2 + sample_studied[0]
    #Smoothing corrected curve before fit
    # smoothed_curve = savgol_filter(sample_studied, 7, 4)
    # smoothed_curve = tls.multiply_data_points(sample_studied, x_coord_we, 10)
    smoothed_curve = savgol_filter(sample_studied, 9, 4)
    maxIndex_smooth = np.argmax(smoothed_curve)
    R0Index_smooth = find_right_perc(smoothed_curve, maxIndex_smooth, CLINICAL_RANGE_PERC)
    # startIndex_smooth = utls.find_left_perc(smoothed_curve, maxIndex_smooth, 0.9)
    stopIndex_smooth = R0Index_smooth + (R0Index_smooth-maxIndex_smooth)
    # bortfeld_proto = tls.createBortfeldCurves(x_coord_fit_we, np.average(sample_studied[0]))
    bortfeld_proto = createBortfeldCurves(x_coord_we, np.average(sample_studied[0]))
    #CONVOLUTE BORTFELD CURVES TO RETRIEVE DOSE-RANGE PROFILE
    count = 0
    final_bortfeld = np.zeros(len(smoothed_curve))
    for x in smoothed_curve:
        norm_bortfeld = (bortfeld_proto[count])*x
        if count <= stopIndex_smooth and count >= maxIndex_smooth :
            final_bortfeld = np.add(final_bortfeld, norm_bortfeld)
        count += 1
    #BORTFELD MODEL FITTING
    #fitting data with Bortfeld model [Bortfeld, MedPhys, 1997] - R0 constraint to be in the position of data max +- 2 bins, fluence between 0 and the value defined in PHI0
    # popt_fit, pcov_fit = curve_fit(utils._bortfeld, x_coord_fit_we, smoothed_curve, bounds=((0., (maxIndex_smooth-3)), (smoothed_curve[0]*110/100, (maxIndex_smooth+3))))
    # perr_fit = np.sqrt(np.diag(pcov_fit))
    # print(popt_fit) 
    # ym = utils._bortfeld(x_coord_fit_we, popt_fit[0], popt_fit[1])
    # ym = utls._bortfeld(x_coord_fit_we, 1, maxIndex_smooth*TO_WE)
    smoothed_curve = multiply_data_points(smoothed_curve, x_coord_we, 10)
    ym = multiply_data_points(final_bortfeld, x_coord_we, 10)
    #ym = final_bortfeld
    max_fit = np.amax(ym)
    # norm_level_fit = (max_fit-ym[0])/2 + ym[0]
    norm_index_fit = index_norm *10
    #norm_index_fit = index_norm
    norm_value_fit = ym[norm_index_fit]
    #norm_fit = max_fit/max_after_correction
    norm_fit = norm_value_fit
    ym_norm = []
    for el_fit in ym:
        ym_norm.append(el_fit/norm_fit)
    ym_norm = np.asarray(ym_norm)
    # ym_norm = ym
    #analysing curve based on smoothed version
    peak_index = 0
    peak_pos = 0
    plt_mean = 0
    peak_val = 0 
    cl_range = 0
    back_cl_range = 0
    if bortfeld:
        peak_val = np.amax(ym_norm) # peak value
        peak_index = np.where(ym_norm == peak_val)[0][0] 
        peak_pos = x_coord_fit_we[peak_index]
        #peak_pos = x_coord_we[peak_index]
        plt_mean = np.average(ym_norm[:3])
        cl_range = find_cl_range(ym_norm, x_coord_fit_we, peak_index, CLINICAL_RANGE_PERC)
        back_cl_range = find_cl_range_back(ym_norm, x_coord_fit_we, peak_index, CLINICAL_RANGE_PERC)
        #cl_range = utls.find_cl_range(ym_norm, x_coord_we, peak_index, 0.8)
        #back_cl_range = utls.find_cl_range_back(ym_norm, x_coord_we, peak_index, 0.8)
    else:
        peak_val = np.amax(smoothed_curve) # peak value
        peak_index = np.where(smoothed_curve == peak_val)[0][0]  # peak index in smoothed curve 
        peak_pos = x_coord_fit_we[peak_index]
        #peak_pos = x_coord_we[peak_index]
        plt_mean = np.average(smoothed_curve[:3])
        cl_range = find_cl_range(smoothed_curve, x_coord_fit_we, peak_index, CLINICAL_RANGE_PERC)
        back_cl_range = find_cl_range_back(smoothed_curve, x_coord_fit_we, peak_index, CLINICAL_RANGE_PERC)
        #cl_range = utls.find_cl_range(smoothed_curve, x_coord_we, peak_index, 0.8)
        #back_cl_range = utls.find_cl_range_back(smoothed_curve, x_coord_we, peak_index, 0.8)
    
    peak_plt_ratio = peak_val/plt_mean # peak to plateau ratio
    
    rawcoord_list = x_coord_we.tolist()  # converting arrays to lists
    rawdata_list = sample_studied.tolist()
    coord_list = x_coord_fit_we.tolist()
    #coord_list = x_coord_we.tolist()
    ym_list = ym_norm.tolist()
    smooth_list = smoothed_curve.tolist()
    coordinates_fit = []
    fit_data = []
    if bortfeld : 
        coordinates_fit = coord_list
        fit_data = ym_list
    else :
        coordinates_fit = coord_list
        fit_data = smooth_list
    # composing dictionary of results
    results = {
        "peak_pos":{"value":float(peak_pos), "unit":"cm w.e."},
        "pp_ratio":{"value":float(peak_plt_ratio),"unit":" "},
        "cl_range":{"value":float(cl_range),"unit":"cm w.e."},
        "peak_width":{"value":float(cl_range - back_cl_range),"unit":"cm w.e."},
        "coordinates_raw":rawcoord_list,
        "raw_data":rawdata_list,
        "coordinates_fit": coordinates_fit,
        "fit_data": fit_data,
    }
    return results


def meanPos(data, coord):
    mean_position = 0
    weights = 0
    for idx in range(len(data)):
        mean_position += data[idx]*coord[idx]
        weights += data[idx]
    result = mean_position/weights
    return result

def calc_sides_withBS(x, y, value):
    max_val = np.amax(x)
    max_index = np.where(x == max_val)[0][0]
    #calculating baseline
    bs = 0
    bs_sx = np.mean(x[:5])
    bs_dx = np.mean(x[-5:])
    bs_els = {"maxEl": max([bs_dx, bs_sx]), "minEl": min([bs_dx, bs_sx])}
    bsRatio = bs_els["maxEl"]/bs_els["minEl"]
    if bsRatio < 2:
        bs = np.mean([bs_dx, bs_sx])
    else:
        bs = bs_els["minEl"]
    level_max = value*(max_val-bs) + bs
    left_half = 0
    right_half = 0
    for cnt in range(0, max_index):
        val = x[cnt]
        if(val >= level_max):
            w1 = np.abs(value - val/max_val)
            w2 = np.abs(value - x[cnt - 1]/max_val)
            pos1 = y[cnt]
            pos2 = y[cnt - 1]
            left_half = ((pos1*w2) + (pos2*w1))/(w1+w2)
            break
    for cnt2 in range(max_index+1, len(x)-1):
        val2 = x[cnt2]
        if(val2 <= level_max):
            w3 = np.abs(value - val2/max_val)
            w4 = np.abs(value - x[cnt2 - 1]/max_val)
            pos3 = y[cnt2]
            pos4 = y[cnt2 - 1]
            right_half = ((pos3*w4) + (pos4*w3))/(w3+w4)
            break
    values = np.array([left_half, right_half])
    return values

def find_value(x, value):
    count = 0
    for val in x :
        if val >= value :
            return count, val
        else :
            count = count + 1

def find_right_perc(x, max_index, perc):
    count = 0
    max = x[max_index]
    max_perc = perc*max
    for val in x[max_index+1:] :
        if val <= max_perc :
            w1 = np.abs(perc - val/max)
            w2 = np.abs(perc - x[max_index + count]/max)
            pos1 = count + 1 + max_index
            pos2 = count + max_index
            point = ((pos1*w2) + (pos2*w1))/(w1+w2)
            return point
        else:
            count = count + 1

def createBortfeldCurves(ranges, fluence):
    ym = []
    for x in ranges:
        # print('Curve for range {}'.format(x))
        new_curve = _bortfeld(ranges, fluence, x)
        ym.append(np.nan_to_num(new_curve))
    return np.array(ym)

def _bortfeld(x, PHI0, R0):
    #R0 = constants.BRAGG_ALPHA*np.power(E0, constants.BRAGG_P)
    rho = 3.2 #g/cm3 aluminum/copper mix density
    #rho = 1 #g/cm3 aluminum
    #rho = 1 #g/cm3 water density
    #sigma for mono-energetic beam
    sigma_mono = BRAGG_SIGMA_MONO_FACT * np.power(R0,BRAGG_SIGMA_MONO_EXP)
    #sigma correction for energy spectrum
    sigma_E0 = BRAGG_SIGMA_E0_FACTOR*np.power(R0/BRAGG_ALPHA, 1/BRAGG_P)
    #cumulative sigma
    sigma = np.sqrt(sigma_mono**2 + sigma_E0**2 * BRAGG_ALPHA**2 * BRAGG_P**2 * np.power(R0/BRAGG_ALPHA, 2 - 2/BRAGG_P))
    #auxiliary variables
    zeta = (R0-x)/sigma
    inverse_p = 1/BRAGG_P
    numerator = np.exp(-(zeta**2)/4)*np.power(sigma, inverse_p)*BRAGG_GAMMA_FUNC
    denominator = np.sqrt(2*np.pi)*rho*BRAGG_P*np.power(BRAGG_ALPHA,inverse_p)*(1+BRAGG_BETA*R0)
    #parabolic cylinder function values
    pbdv_0 = spspec.pbdv(-1/BRAGG_P, -zeta)[0]
    pbdv_1 = spspec.pbdv(-1/BRAGG_P-1, -zeta)[0]
    factor_1 = pbdv_0/sigma
    factor_2 = pbdv_1*(BRAGG_BETA/BRAGG_P + BRAGG_BETA*BRAGG_GAMMA + BRAGG_EPSILON_MAX/R0)

    return PHI0*(numerator/denominator)*(factor_1 + factor_2) #Bortfeld fitting function


def multiply_data_points(points, coords, factor):
    new_array = np.asarray([])
    for x in range(len(points)-1):
        cc = [coords[x], coords[x+1]]
        pp = [points[x], points[x+1]]
        # Calculate the coefficients. This line answers the initial question. 
        coefficients = np.polyfit(cc, pp, 1)
        # Let's compute the values of the line...
        polynomial = np.poly1d(coefficients)
        x_axis = np.linspace(coords[x],coords[x+1],factor)
        y_axis = polynomial(x_axis)
        new_array = np.concatenate((new_array, y_axis), axis = 0)    
    return new_array

def find_cl_range(x, y, max_index, perc):
    #print('Searching clinical range ...')
    count = 0
    max = x[max_index]
    max_perc = perc*max
    for val in x[max_index+1:] :
        if val <= max_perc :
            w1 = np.abs(perc - val/max)
            w2 = np.abs(perc - x[max_index + count]/max)
            pos1 = y[count + 1 + max_index]
            pos2 = y[count + max_index]
            cl_range = ((pos1*w2) + (pos2*w1))/(w1+w2)
            return cl_range
        else:
            count = count + 1

def find_cl_range_back(x, y, max_index, perc):
    count = 0
    max = x[max_index]
    max_perc = perc*max
    for cnt in range(1, max_index) :
        val = x[max_index-cnt]
        if val <= max_perc :
            w1 = np.abs(perc - val/max)
            w2 = np.abs(perc - x[max_index - cnt + 1]/max)
            pos1 = y[max_index-cnt]
            pos2 = y[max_index-cnt + 1]
            back_perc = ((pos1*w2) + (pos2*w1))/(w1+w2)
            return back_perc