from __future__ import print_function


import numpy as np
import sys
import scipy.ndimage as ndi
import math
import matplotlib.pyplot as plt

import isinmrtbx.inout.writers as writers


AUTHOR = 'Tomas Psorn'
CONTACT = 'tomaspsorn@isibrno.cz'


def markValidPoints1(dataIn, metric='max', span=.05, bord=.3):
    """
    Creates a mask discarding noise suspect values

    Parameters
    ----------
    dataIn
    medOrMax
    ~ Switch between classifier functions
    trashold
    ~ trashold value used by classifier

    Returns
    -------
    validPoints
    ~ the mask
    """

    _dataIn = dataIn.copy()

    dim0 = _dataIn.shape[0]

    _dataIn += np.abs(np.amin(_dataIn))

    # plt.figure()
    # plt.subplot(2,3,2)
    # plt.title('+ abs(min(data))')
    # plt.plot(_dataIn)

    _dataIn[0 : int(dim0 * bord)] = 0.0
    _dataIn[int(dim0 * (1-bord)):] = 0.0

    # plt.subplot(2,3,3)
    # plt.title('crop)')
    # plt.plot(_dataIn)

    validPoints = np.ones(_dataIn.shape[0],dtype='i2')

    if metric == 'med':
        median = np.median(_dataIn)
        scale = median

    if metric == 'max':
        thrashold_up = np.amax(_dataIn)
        thrashold_down = np.amin(_dataIn[np.nonzero(_dataIn)])

    thrashold_up *= 1. -span # scale up the trashold
    thrashold_down *= 1. + span  # scale up the trashold

    # plt.subplot(2,3,4)
    # plt.title('scale + thrasholds)')
    # plt.plot(_dataIn)
    # plt.axhline(y=thrashold_up, xmax=dim0, xmin=0, color='black')
    # plt.axhline(y=thrashold_down, xmax=dim0, xmin=0, color='black')

    for _dim0 in range(0,dim0):
        if _dataIn[_dim0] < thrashold_down or _dataIn[_dim0] > thrashold_up:
            validPoints[_dim0] = 0

    # plt.subplot(2,3,5)
    # plt.title('result')
    # plt.plot(_dataIn)
    # plt.stem(_dataIn * validPoints)
    # plt.savefig('results/epi-srepi-vs-grepi/vp/' + time.strftime("%H:%M:%S") + '.png')
    # plt.close()

    return validPoints


def markValidPoints2(dataIn, medOrMax='med'):
    """
    2D version of markValidPoints iterates trough second dimension
    performing a markValidPoints1 over each line

    Parameters
    ----------
    dataIn
    medOrMax

    Returns
    -------
    validPoints
    """

    dim0 = dataIn.shape[0]

    for line in range(0,dim0):
        if line == 0:
            validPoints = markValidPoints1(dataIn[line,:])
        else:
            validPoints *= markValidPoints1(dataIn[line,:])

    return validPoints


def linearFitSafe(dataIn):

    dim0 = dataIn.shape[0]
    dim0_nonZero = np.count_nonzero(dataIn)

    x = np.zeros (dim0_nonZero, dtype='f4', order='F')
    y = np.zeros (dim0_nonZero, dtype='f4', order='F')

    count = 0

    for _dim0 in range(0,dim0):
        if dataIn[_dim0] != 0.0:
            x[count] = _dim0
            y[count] = dataIn[_dim0]
            count += 1

    if count >= 2:
        offset, slope, offset_sd, slope_sd = linearFit(x,y)
    else:
        offset, slope, offset_sd, slope_sd = 0.0, 0.0, 0.0, 0.0

    return offset, slope, offset_sd, slope_sd

def linearFit(x,y):

    dim0 = float(x.shape[0])

    sum_x, sum_y, sxx, syy, sxy, chi2 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    for i in range(0,int(dim0)):
        sum_x += x[i]
        sum_y += y[i]
        sxx += x[i] * x[i]
        syy += y[i] * y[i]
        sxy += x[i] * y[i]

    slope = (sxy - sum_x * sum_y / dim0) / (sxx - sum_x * sum_x / dim0)
    offset = ( sum_y - slope * sum_x) / dim0

    for i in range(0, int(dim0)):
        chi2 += ((y[i] - offset - slope * x[i])*(y[i]-offset-slope*x[i]))

    if dim0 > 2:
        offset_sd = np.sqrt(chi2/(dim0-2) * sxx/(dim0*sxx-sum_x*sum_y))
        slope_sd = np.sqrt(chi2/(dim0-2) *   dim0/(dim0*sxx-sum_x*sum_x))
    else:
        offset_sd = 0.0
        slope_sd = 0.0

    return offset, slope, offset_sd, slope_sd


def mirrorOddLines(dataIn, axis=None):

    if axis == None:
        print("Function: mirrorOddLines")
        print("Problem: axis not specified")
        sys.exit(1)

    dataDim = dataIn.shape
    dataOut = np.zeros_like(dataIn)

    if len(dataDim) == 2:
        if axis == 0:
            dataOut [:, :2] = dataIn[:, :2]
            dataOut [:, 1::2] = dataIn[::-1,1::2]
        elif axis == 1:
            dataOut [::2,:] = dataIn[::2,:]
            dataOut [1::2,:] = dataIn[1::2,::-1]
        else:
            print("Function: mirrorOddLines, problem: invalid axis")
            sys.exit(1)

    elif len(dataDim) == 3:
        if axis == 0:
            dataOut [:, :2, :] = dataIn[:, :2, :]
            dataOut [:, 1::2, :] = dataIn[::-1, 1::2, :]
        elif axis == 1:
            dataOut[::2, :, :] = dataIn[::2, :, :]
            dataOut[1::2, :, :] = dataIn[1::2, ::-1, :]
        elif axis == 2:
            print("To be implemented yet")
            sys.exit(1)
        else:
            print("Function: mirrorOddLines, problem: invalid axis")
            sys.exit(1)

    else:
        print("Function: mirrorOddLines, can only process 2D and 3D data")
        sys.exit(1)


    return dataOut

def fft1dInPlace (dataIn, axis):

    dim0 = dataIn.shape[0]
    dim1 = dataIn.shape[1]

    fftShiftFactor = np.exp(np.complex(0.0,np.pi))

    dataOut = np.zeros_like(dataIn)

    if axis == 0:
        dataIn[0:dim0:2, :] *= fftShiftFactor

        for _dim1 in range(0,dim1):
            dataOut[:,_dim1] = np.fft.fft(dataIn[:,_dim1])

        dataOut[0:dim0:2, :] *= fftShiftFactor

    elif axis == 1:
        dataIn[:, 0:dim1:2] *= fftShiftFactor

        for _dim0 in range(0, dim0):
            dataOut[_dim0,:] = np.fft.fft(dataIn[_dim0,:])

        dataOut[:,0:dim1:2] *= fftShiftFactor

    else:
        print (" prameter axis: ", axis, "not valid")
        sys.exit(1)

    return dataOut

def combineChannels (dataIn, mode="AdaptiveCombine"):
    """
    Description
    -----------
    Function to combine images reconstructed from the signal,
    from individual channels. It uses either adaptive combine,
    or sum of squares.
    Parameters
    ----------
    dataIn
        In format:
        [channels, rows, columns]
        or
        [channels, slices, rows, columns]
    mode
        Mode of combination:
        "AdaptiveCombine",  or it's equivalent 0
        or
        "SumOfSquares",     or it's equivalent 1
    Returns
    -------
    dataOut
        Combined data
        In format:
        [rows, columns]
        or
        [slices, rows, columns]
    """
    defaultMode = "AdaptiveCombine"

    dim = np.shape(dataIn)

    if len(dim) is 3:
        channels_nr = np.shape(dataIn)[0]
        slices_nr = 1
        rows_nr = np.shape(dataIn)[1]
        columns_nr = np.shape(dataIn)[2]
        temp = np.zeros((channels_nr, slices_nr, rows_nr, columns_nr),dataIn.dtype)
        temp[:,0,:,:] = dataIn[:,:,:]
        dataIn = temp
    elif len(dim) is 4:
        channels_nr = np.shape(dataIn)[0]
        slices_nr = np.shape(dataIn)[1]
        rows_nr = np.shape(dataIn)[2]
        columns_nr = np.shape(dataIn)[3]
    else:
        print("Data dimensions ", dim, " are not valid")
        sys.exit(1)

    dataOut = np.zeros((slices_nr, rows_nr, columns_nr),dataIn.dtype)
    midSlice = int(slices_nr/2)

    if mode != 0 \
            and mode != 1 \
            and mode != "AdaptiveCombine"\
            and mode != "SumOfSquares":
        mode = defaultMode

    #   Do adaptive combine
    if mode is "AdaptiveCombine" or mode is 0:
        correlationMatrix = np.zeros((channels_nr,channels_nr),dataIn.dtype)
        for _x in range(0,channels_nr):
            for _y in range(0,channels_nr):
                correlationMatrix[_x,_y] = sum(sum(dataIn[_x,midSlice,:,:] * np.conj(dataIn[_y,midSlice,:,:])))

        val, vec = np.linalg.eigh(correlationMatrix)
        abs_val = np.abs(val)
        max_index = np.argmax(abs_val)
        max_vec = vec[:,max_index]

        for _channel in range(0, channels_nr):
            dataOut = dataOut + dataIn[_channel,:,:,:] * np.conj(max_vec[_channel])

        print("end adaptive")

    #   Do sum of squares
    elif mode is "SumOfSquares" or mode is 1:
        for _slice in range (0, slices_nr):
            for _channel in range(0, channels_nr):
                dataOut[_slice,:,:] = dataOut[_slice,:,:] + dataIn[_channel,_slice,:,:] * np.conj(dataIn[_channel,_slice,:,:])
        dataOut = np.sqrt(dataOut)

        print("end squares")

    if len(dim) is 3:
        dataOut = np.squeeze(dataOut,0)

    return dataOut

def zeroFilling(dataIn, loc, size):
    """
    Description
    -----------
    Function to complete non-complete k-spaces with zeros,
    prior to 2D IFFT
    Parameters
    ----------
    dataIn
        In format:
        [rows, columns]
        or
        [channels, rows, columns]
        or
        [channels, slices, rows, columns]
    axis
        Determine an axis along which axis the zero filling is
        to be done. Options: 0,1
    Returns
    -------
    dataOut
        In format:
        [rows, columns]
        or
        [channels, rows, columns]
        or
        [channels, slices, rows, columns]
    """

    dim0_out = size[0]
    dim1_out = size[1]
    dim0_in = dataIn.shape[0]
    dim1_in = dataIn.shape[1]
    dim0_diff = abs(dim0_out - dim0_in)
    dim1_diff = abs(dim1_out - dim1_in)

    vert = loc[0]
    hor = loc[1]

    dataOut = np.zeros((dim0_out,dim1_out),dataIn.dtype)

    if vert == 'up' and hor == 'left':
        dataOut[0:-dim0_diff,0:-dim1_diff] = dataIn
    elif vert == 'up' and hor == 'mid':
        dataOut[0:-dim0_diff,int(dim1_out/2)-int(dim1_in/2):int(dim1_out/2)+int(dim1_in/2)] = dataIn
    elif vert == 'up' and hor == 'right':
        dataOut[0:-dim0_diff,-(dim1_out-dim1_diff):] = dataIn
    elif vert == 'mid' and hor == 'left':
        dataOut[int(dim0_out/2)-int(dim0_in/2):int(dim0_out/2)+int(dim0_in/2),0:-dim1_diff] = dataIn
    elif vert == 'mid' and hor == 'mid':
        dataOut[int(dim0_out/2)-int(dim0_in/2):int(dim0_out/2)+int(dim0_in/2),int(dim1_out/2)-int(dim1_in/2):int(dim1_out/2)+int(dim1_in/2)] = dataIn
    elif vert == 'mid' and hor == 'right':
        dataOut[int(dim0_out/2)-int(dim0_in/2):int(dim0_out/2)+int(dim0_in/2),-(dim1_out-dim1_diff):] = dataIn
    elif vert == 'down' and hor == 'left':
        dataOut[-(dim0_out-dim0_diff):,0:-dim1_diff] = dataIn
    elif vert == 'down' and hor == 'mid':
        dataOut[-(dim0_out-dim0_diff):,int(dim1_out/2)-int(dim1_in/2):int(dim1_out/2)+int(dim1_in/2)] = dataIn
    elif vert == 'down' and hor == 'right':
        dataOut[-(dim0_out-dim0_diff):,-(dim1_out-dim1_diff):] = dataIn
    else:
        print(zeroFilling.__doc__)

    return dataOut

def snr (dataIn, sigRoi, noiseRoi):

    signal = dataIn[sigRoi[0]:sigRoi[1],sigRoi[2]:sigRoi[3]]
    noise  = dataIn[noiseRoi[0]:noiseRoi[1],noiseRoi[2]:noiseRoi[3]]

    signal_p = np.sum(np.abs(signal)**2)/float(signal.size)
    noise_p = np.sum(np.abs(noise)**2)/float(noise.size)

    snr = 10*np.log10((signal_p-noise_p)/noise_p)

    return snr

def shoeLace(corners):

    n = corners.shape[0]
    area = 0

    for i in range(0,n):
        j = (i + 1) % n
        area += corners[i][0] * corners[j][1]
        area -= corners[j][0] * corners[i][1]


    area = abs(area) * 0.5

    return area


def calcLine (length, slope, offset, start=0, type='f4'):

    out = np.zeros(length, dtype=type)

    for x in range(start,length):
        out[x] = x * slope + offset

    return out

def shiftToPossitive(dataIn):

    minimum = np.amin(dataIn)

    if minimum < 0:
        minimum = np.abs(minimum)
        dataIn += minimum

    return dataIn

def zero_to_nan(values):
    """Replace every 0 with 'nan' and return a copy."""
    return [float('nan') if x==0 else x for x in values]