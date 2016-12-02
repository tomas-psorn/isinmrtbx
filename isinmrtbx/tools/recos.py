#   Compatibility with python 3.x
# from __future__ import print_function
# from __future__ import unicode_literals
# from __future__ import division

import math
import sys
import collections


import matplotlib.pyplot as plt
import numpy as np

from utils import combineChannels
import utils as utils

AUTHOR = 'Tomas Psorn'
CONTACT = 'tomaspsorn@isibrno.cz'


def flashReco(dataObject):

    kspace = dataObject.fid

    dataOut = np.fft.fftshift(np.fft.ifft2(kspace, axes=(0,1)), axes=(0,1))

    return dataOut

"""
EPI reco
"""

def epiReco(dataObject):
    """
    EPI Reco pipeline

    Working with navigator data
    ~ Mirror odd lines
    ~ Do FFT in read direction (step into hybrid space)
    ~ Get phase of adjacent lines in each segments
    ~ Unwrap and fit the phase to a line and get slope and offset

    Working with data
    ~ Every other line is corrected by subtracting the phase slope
    ~ All measured k-space data distributed to a zero, reco-sized matrix
    ~ Grappa reco is used to calculate missing lines if any
    ~ Partial Fourier reco is used to correct for k-space truncation if any
    ~ FFT in phase direction is done
    ~ Multifrequency reco is done to correct for B0 inhomogenities

    Parameters
    ----------
    dataObject

    Returns
    -------

    """


    dimXr = dataObject.method.PVM_EncMatrix[0]  # number of pixels in read direction in the resulting image
    dimYr = dataObject.method.PVM_EncMatrix[1]  # number of phase encoding lines in the raw data.
    dimXi = dataObject.method.PVM_Matrix[0]     # number of samples in read direction in the raw data
    dimYi = dataObject.method.PVM_Matrix[1]     # number of pixels in the PE direction in the resulting image
    dimZ = dataObject.acqp.NI   # number of slices
    dimRep = dataObject.acqp.NR     #number of repetitions
    dimSeg = dataObject.method.NSegments    # number of segments
    dimCh = dataObject.getSelectedReceivers()       # number of receivers
    dimYrSeg = dimYr // dimSeg

    phaseCorr = np.zeros((dimCh, dimYrSeg , dimXr, dimZ * dimSeg), dtype='f4') # An array to store the phase correction information of every data point.

    # if there are navigators measured, use them for phase correction, if not, use Bruker stuff
    try:
        omNavigators_nr = dataObject.method.OMNavigators    # how many navigators were measured
        dimRep -= omNavigators_nr   # navigators are stored in a dimension together with repetitions
    except:
        omNavigators_nr = 0

    dimNavFrames = dimZ * omNavigators_nr * dimSeg
    dimNoNavFrames = dataObject.fid.shape[3] - dimNavFrames

    navData = np.zeros((dimCh, dimYrSeg, dimXr, dimNavFrames), dtype=dataObject.fid.dtype)
    rawData = np.zeros((dimCh, dimYrSeg, dimXr, dimNoNavFrames), dtype=dataObject.fid.dtype)


    for _dimCh in range(0, dimCh):
        navData[_dimCh,:,:,:], rawData[_dimCh,:,:,:] = epi_splitEpiData(dataObject.fid[_dimCh, :, :,:], dimZ, dimSeg, omNavigators_nr)  # split fid to navigators and data
        phaseCorr[_dimCh,:,:,:] = epi_calcPhaseCor2(navData[_dimCh], dimXr, dimYr, dimZ, dimSeg)  # it'll just read the bruker slopes and offsets


    resultChan = np.zeros((dimCh, dimYi, dimXi, dimZ), dtype=rawData.dtype) # array to store reco results for each channel
    resultRep = np.zeros((dimYi, dimXi, dimZ, dimRep), dtype=rawData.dtype) # array to store reco results for each repetition

    # reco goes volume by volume, since some intra volume correction might be used eventually
    for _dimRep in range(0, dimRep):
        for _dimCh in range(0, dimCh):
            # get the particular volume and do the reco
            volume = rawData[_dimCh,:,:, _dimRep * dimZ * dimSeg : (_dimRep+1) * dimZ * dimSeg]
            resultChan[_dimCh, :, :, :] = epiRecoRaw (volume, dataObject, phaseCorr[_dimCh,:,:,:])

            plt.figure()
            plt.imshow(np.abs(resultChan[_dimCh,:,:,0]))
            plt.show()

        resultRep[:,:,:, _dimRep] = combineChannels(resultChan)

    return resultRep



def epiRecoRaw (kspace, dataObject, phaseCorr):

    XdimI = dataObject.method.PVM_Matrix[0] # image size
    YdimI = dataObject.method.PVM_Matrix[1]
    segments_nr = dataObject.method.NSegments

    kspace = utils.mirrorOddLines(kspace, axis=1)


    hybrid = np.fft.ifft2(kspace, axes=(1,))



    hybrid = np.fft.fftshift(hybrid, axes=1)


    hybrid =  phaseCorEpi(hybrid, phaseCorr)


    hybrid = sortEpi(hybrid, dataObject)


    kspace = np.fft.fft2(hybrid, axes=(1,))



    plt.figure()
    plt.imshow(np.log(np.abs(kspace[:,:,0]))**0.1, interpolation=None)
    plt.show()

    image = (np.fft.ifft2(hybrid, s=(YdimI,), axes=(0,)))

    image = np.fft.fftshift(image, axes=0)


    return image



def epi_calcPhaseCor(navData,dimXr, dimYr, dimZ, dimSeg):
    '''
    Function to get phase correction coefficients.
    This can be done in two ways.

    1/  Navigator scans are measured
    2/  If navigator scans are not a part of the dataset,
        then the Bruker calculated coefficients are used

    Parameters
    ----------
    navData
    dimXr
    dimYr
    dimZ
    dimSeg

    Returns
    -------
    phaseCorr
    '''

    dimHigh = dimZ * dimSeg
    dimYrSeg = dimYr // dimSeg

    slopes = np.zeros((dimYrSeg, dimHigh), dtype='f4')
    offsets = np.zeros((dimYrSeg, dimHigh), dtype='f4')

    # test if navigator scan is present
    if navData is None:
        mode = 'obj'
    else:
        mode = 'nav'

    # calculate phase correction coefficients
    if mode == 'nav':

        # mirror odd lines

        navData = utils.mirrorOddLines(navData,axis=1)

        # FFT in read dirrection
        navData = np.fft.fftshift(np.fft.ifft2(navData, axes=(1,)), axes=1)


        for _dimHigh in range(0, dimHigh):

            _dimSeg = _dimHigh % dimSeg

            for _dimYrSeg in range(0, dimYrSeg, 2):

                phaseEven = np.angle(navData[_dimYrSeg, :, _dimHigh])
                phaseOdd = np.angle(navData[_dimYrSeg + 1, :, _dimHigh])

                # get a reference for a intersegment correct

                # if first segment rewrite phase init reset phase diff to zero
                if _dimSeg == 0 and _dimYrSeg == 0 and _dimHigh == 0:
                    phaseFirst = phaseEven #todo just in case for now
                    phaseFirst_unwrap = np.unwrap(phaseFirst)
                    phaseFirst_unwrap *= utils.markValidPoints1(phaseFirst_unwrap)
                    offsetFirst, slopeFirst, offset_sd, slope_sd = utils.linearFitSafe(phaseFirst_unwrap)

                if _dimSeg == 0 and _dimYrSeg == 0:
                    phaseEven_unwrap = np.unwrap(phaseEven)
                    phaseEven_unwrap *= utils.markValidPoints1(phaseEven_unwrap)

                    offsetActuall, slopeActuall, offset_sd, slope_sd = utils.linearFitSafe(phaseEven_unwrap)

                    phaseDiffIntra = offsetFirst - offsetActuall

                    # get phase difference between adjacent lines
                phaseOddEven = phaseOdd - phaseEven

                # unwrap phase
                phaseOddEven = np.unwrap(phaseOddEven)

                phaseOddEven *= utils.markValidPoints1(phaseOddEven)

                offset, slope, offset_sd, slope_sd = utils.linearFitSafe(phaseOddEven)

                offsets[_dimYrSeg+1, _dimHigh] = offset - phaseDiffIntra
                slopes [_dimYrSeg+1, _dimHigh] = slope
                offsets[_dimYrSeg, _dimHigh] -= phaseDiffIntra



    #get bruker correction coefficients
    if mode == 'obj':
        print("EPI reco which would use bruker corr coefs is not developed yet")

    phaseCorr = np.zeros((dimYrSeg, dimXr, dimHigh),
                         dtype='f4')  # An array to store the phase correction information of every data point.

    for _dimHigh in range(0,dimHigh):
        for _dimYrSeg in range(0,dimYrSeg):
            phaseCorr[_dimYrSeg, :, _dimHigh] = utils.calcLine(dimXr, slopes[_dimYrSeg,_dimHigh], offsets[_dimYrSeg,_dimHigh])

    return phaseCorr

def epi_calcPhaseCor2(navData,dimXr, dimYr, dimZ, dimSeg):
    '''
    Function to get phase correction coefficients.
    This can be done in two ways.

    1/  Navigator scans are measured
    2/  If navigator scans are not a part of the dataset,
        then the Bruker calculated coefficients are used

    Parameters
    ----------
    navData
    dimXr
    dimYr
    dimZ
    dimSeg

    Returns
    -------
    phaseCorr
    '''

    dimHigh = dimZ * dimSeg
    dimYrSeg = dimYr // dimSeg


    # test if navigator scan is present
    if navData is None:
        mode = 'obj'
    else:
        mode = 'nav'

    # calculate phase correction coefficients
    if mode == 'nav':

        slopesOddEven = np.zeros((dimYrSeg, dimHigh), dtype='f4')
        offsetsOddEven = np.zeros((dimYrSeg, dimHigh), dtype='f4')
        slopesInterEven = np.zeros((dimYrSeg, dimHigh), dtype='f4')
        offsetInterEven = np.zeros((dimYrSeg, dimHigh), dtype='f4')
        slopesInterOdd = np.zeros((dimYrSeg, dimHigh), dtype='f4')
        offsetsInterOdd = np.zeros((dimYrSeg, dimHigh), dtype='f4')

        phaseFirstEven = np.zeros((dimYrSeg // 2, dimXr), dtype='f4')
        phaseFirstOdd = np.zeros((dimYrSeg // 2, dimXr), dtype='f4')

        # mirror odd lines

        navData = utils.mirrorOddLines(navData,axis=1)

        # FFT in read dirrection
        navData = np.fft.fftshift(np.fft.ifft2(navData, axes=(1,)), axes=1)


        for _dimHigh in range(0, dimHigh):

            _dimSeg = _dimHigh // dimSeg

            for _dimYrSeg in range(0, dimYrSeg, 2):

                phaseEven = np.angle(navData[_dimYrSeg, :, _dimHigh])
                phaseOdd = np.angle(navData[_dimYrSeg + 1, :, _dimHigh])

                # get a reference for a intersegment correct

                # if first segment rewrite phase init reset phase diff to zero
                if _dimSeg == 0:
                    phaseFirstEven[_dimYrSeg // 2, :] = phaseEven
                    phaseFirstOdd[_dimYrSeg // 2, :] = phaseOdd

                # get phase difference between adjacent lines
                phaseOddEven = phaseOdd - phaseEven
                phaseInterEven = phaseFirstEven[_dimYrSeg // 2, :] - phaseEven
                phaseInterOdd = phaseFirstOdd[_dimYrSeg // 2, :] - phaseOdd

                # unwrap phase
                phaseOddEven = np.unwrap(phaseOddEven)
                phaseInterEven = np.unwrap(phaseInterEven)
                phaseInterOdd = np.unwrap(phaseInterOdd)

                phaseOddEven *= utils.markValidPoints1(phaseOddEven)
                phaseInterEven *= utils.markValidPoints1(phaseInterEven)
                phaseInterOdd *= utils.markValidPoints1(phaseInterOdd)

                offset, slope, offset_sd, slope_sd = utils.linearFitSafe(phaseOddEven)
                offsetsOddEven[_dimYrSeg+1, _dimHigh] = offset
                slopesOddEven [_dimYrSeg+1, _dimHigh] = slope

                offset, slope, offset_sd, slope_sd = utils.linearFitSafe(phaseInterEven)
                offsetInterEven[_dimYrSeg, _dimHigh] = offset
                slopesInterEven [_dimYrSeg, _dimHigh] = slope

                offset, slope, offset_sd, slope_sd = utils.linearFitSafe(phaseInterOdd)
                offsetsInterOdd[_dimYrSeg+1, _dimHigh] = offset
                slopesInterOdd [_dimYrSeg+1, _dimHigh] = slope




    #get bruker correction coefficients
    if mode == 'obj':
        print("EPI reco which would use bruker corr coefs is not developed yet")

    phaseCorr = np.zeros((dimYrSeg, dimXr, dimHigh),
                         dtype='f4')  # An array to store the phase correction information of every data point.

    for _dimHigh in range(0,dimHigh):
        for _dimYrSeg in range(0,dimYrSeg):
            phaseCorr[_dimYrSeg, :, _dimHigh] = utils.calcLine(dimXr, slopesOddEven[_dimYrSeg,_dimHigh], offsetsOddEven[_dimYrSeg,_dimHigh])
            phaseCorr[_dimYrSeg, :, _dimHigh] += utils.calcLine(dimXr, slopesInterEven[_dimYrSeg, _dimHigh], offsetInterEven[_dimYrSeg, _dimHigh])
            phaseCorr[_dimYrSeg, :, _dimHigh] += utils.calcLine(dimXr, slopesInterOdd[_dimYrSeg, _dimHigh], offsetsInterOdd[_dimYrSeg, _dimHigh])

    return phaseCorr



def calcIntensityCorEpi(dataObject, **kwargs):

    NI = dataObject.acqp.NI
    NSegments = dataObject.method.NSegments

    intensityCor = np.zeros((NSegments, NI), dtype='f4')

    dim2 = NI * NSegments

    try:
        navigator = kwargs['navigator']
        mode = 'nav'
    except:
        mode = 'obj'
        pass



    if mode == 'nav':

        for _dim2 in range(0,dim2):
            _frame, _segment = divmod(_dim2, NSegments)

            segment = navigator[:,:, (_frame * NSegments) + _segment]

            intensityCor[_segment, _frame] = np.sum(np.abs(segment))

        for _frame in range(0, NI):
            intensityCor[:,_frame] /= np.amax(intensityCor[:,_frame])
            pass

    else:
        pass

    return intensityCor


def phaseCorEpi(dataIn, phaseCorr):

    linesPerSeg = dataIn.shape[0]
    dimXr = dataIn.shape[1]
    dimHigh = dataIn.shape[2]

    dataOut = np.zeros_like(dataIn)


    for _dimHigh in range(0, dimHigh):
        for _linesPerSeg in range(0, linesPerSeg):
            for _dimXr in range(0, dimXr):
                phase = np.complex(0.0, phaseCorr[_linesPerSeg , _dimXr , _dimHigh])

                dataOut[_linesPerSeg , _dimXr , _dimHigh] = dataIn[_linesPerSeg , _dimXr , _dimHigh] * np.exp(-1.0 * phase)

    return dataOut

def intensityCorEpi(dataIn, intensityCor):

    NSegments = intensityCor.shape[0]
    NI = intensityCor.shape[1]

    dim2 = dataIn.shape[2]

    for _dim2 in range(0, dim2):
        _frame, _segment = divmod(_dim2, NSegments)
        dataIn[:, :, (_frame * NSegments) + _segment] /= intensityCor[_segment,_frame]

    return dataIn

def epi_splitEpiData(dataIn, dimZ, dimSeg, dimNav ):
    '''
    Split epi data set to navigators and raw data, if there is no navigator data,
    it returns None + rawData
    Parameters
    ----------
    dataIn
    dimZ
    dimSeg
    dimNav

    Returns
    -------

    '''
    data = np.squeeze(dataIn)

    if dimNav > 0:
        navigator_data = data[ :, :, 0: dimZ * dimNav * dimSeg]
        raw_data = data[:,:, dimZ * dimNav * dimSeg:]
    else:
        navigator_data = None
        raw_data = data

    return navigator_data, raw_data

def prepareMultifrequencyReco(brukerObj, numberOfColumns, numberOfPElines, swh):

    # Get esential parameters from bruker data object
    fieldMap = brukerObj.data2dseq.data


    # Derive rest of parameters
    dwellTime_ms = 1.0 / float(swh)
    print (dwellTime_ms)
    acqTimeLine = dwellTime_ms * numberOfColumns
    minFieldMap = np.amin(fieldMap)
    maxFieldMap = np.amax(fieldMap)
    bandwidthPerVoxelPE = 1.0/acqTimeLine/numberOfPElines
    numberOfFrequencies = int(2.0 * (maxFieldMap - minFieldMap)/bandwidthPerVoxelPE)
    frequencyStep = (maxFieldMap - minFieldMap)/numberOfFrequencies

    # Alloc array to store all the offset frequencies
    frequenciesArray = np.zeros(frequencyStep, dtype=minFieldMap.dtype, order='F')

    # Fill frequenciesArray
    for i in range(0,numberOfFrequencies):
        frequenciesArray[i] = minFieldMap + i * frequencyStep + .5 * frequencyStep

    return frequenciesArray

def sortEpi(dataIn, dataObject):

    dimSeg = dataObject.method.NSegments
    dimZ = dataObject.acqp.NI

    PVM_EncSteps1 = np.asarray(dataObject.method.PVM_EncSteps1)
    PVM_EncSteps1 += np.abs(np.amin(PVM_EncSteps1))

    dimYrSeg = dataIn.shape[0]
    dimXr = dataIn.shape[1]
    dimHigh = dataIn.shape[2]

    dimYr = dimYrSeg * dimSeg

    dataOut = np.zeros((dimYr,dimXr, dimZ), dtype=dataIn.dtype)

    for _dimHigh in range(0, dimHigh):

        _dimZ = int(_dimHigh / dimSeg)
        _dimSeg = _dimHigh % dimSeg

        for _dimYrSeg in range(0, dimYrSeg):

            _PVM_EncSteps1 = _dimSeg *  dimYrSeg + _dimYrSeg
            _dataOut = PVM_EncSteps1[_PVM_EncSteps1]
            dataOut[_dataOut, :, _dimZ] = dataIn[_dimYrSeg, :, _dimHigh]

    return dataOut

"""
Gridding stuff
"""

def gridding_dcf1(traj, fog, kerWidth, iter_nr=5):

    out = np.ones_like(traj,order='F')

    for _iter in range(0,iter_nr):

        grid_tmp = gridding1_1(out, traj, fog, kerWidth, trim=True, doDcf=False)

        weights_tmp = degridding1_1(grid_tmp, traj, fog, kerWidth)

        for i in range(0,len(weights_tmp)):
            if weights_tmp[i] == 0.: out[i] = 0.
            else: out[i] /= weights_tmp[i]


    return out


def gridding1_1 (dataIn, trajIn, fog, kernelWid_n, trim=True, doDcf = True, iterations=4):
    """
    _n - points
    _k - kspace
    _g - gridspace
    """

    if doDcf:
        densityCompensationFactors = gridding_dcf1(trajIn, fog, kernelWid_n, iterations)
        dataIn *= densityCompensationFactors

    kernelSamples_nr = 1000
    beta = np.pi*np.sqrt( math.pow(kernelWid_n,2)/math.pow(fog,2)*math.pow((fog-0.5),2)-0.8 )
    kernel = np.kaiser(2 * kernelSamples_nr, beta)[-kernelSamples_nr:]

    dim0_n = len(dataIn)
    dim0_g = dim0_n * fog

    kernelWid_k = kernelWid_n / float(dim0_n)
    kernelWid_k_div2 = kernelWid_k / 2.

    ddim0_k = 1.0 / float(dim0_n)
    ddim0_g = 1.0/ float(dim0_g)

    kmax = 0.5 + kernelWid_k_div2
    overlap = np.ceil(kernelWid_k_div2 / ddim0_g)
    dataOut = np.zeros(dim0_g + 2 * overlap, dtype=dataIn.dtype, order='F')

    for data_x in range(0, dim0_n):
        data = dataIn[data_x]
        kx = trajIn[data_x]

        for grid_x in range(0, len(dataOut)):

            gx = grid_x * ddim0_g - kmax
            dx = abs(kx - gx)

            if dx < kernelWid_k_div2:
                ind = np.floor(kernelSamples_nr* dx / kernelWid_k_div2)
                ker = kernel[ind]
                dataOut[grid_x] += data * ker

    if trim:
        # Add overlaping parts
        dataOut[-2*overlap:-overlap] += dataOut[0:overlap]
        dataOut[overlap: 2*overlap] += dataOut[-overlap:]
        # Trim
        dataOut = dataOut[overlap:-overlap]

    return dataOut

def degridding1_1 (dataIn, trajIn, fog, kernelWid_n):
    """
    From grid coordinates to traj coordinates
    """

    kernelSamples_nr = 1000

    beta = np.pi * np.sqrt(math.pow(kernelWid_n, 2) / math.pow(fog, 2) * math.pow((fog - 0.5), 2) - 0.8)
    kernel = np.kaiser(2 * kernelSamples_nr, beta)[-kernelSamples_nr:]

    dim0_n = len(trajIn)
    dim0_g = len(dataIn)

    kernelWid_k_n = kernelWid_n / float(dim0_n)
    kernelWid_k_div2_n = kernelWid_k_n / 2.

    ddim0_g = (1 + kernelWid_k_n) / float(dim0_g)
    ddim0_k = 1 / float(dim0_n)

    kmax = 0.5 + kernelWid_k_div2_n

    dataOut = np.zeros(dim0_n, dtype=dataIn.dtype, order='F')



    for grid_x in range(0, dim0_g):

        gx = grid_x * ddim0_g - kmax
        data_g = dataIn[grid_x]

        for data_x in range(0, len(trajIn)):

            kx = trajIn[data_x]
            dx = abs(kx - gx)

            if dx < kernelWid_k_div2_n:
                ind = np.floor(kernelSamples_nr * dx / kernelWid_k_div2_n)
                ker = kernel[ind]
                dataOut[data_x] += data_g * ker

    return dataOut


def gridMinMax(x, maximum, radius):

    out = collections.namedtuple('MinMax',['min','max'])
    out.min = int(np.ceil(x-radius))
    out.max = int(np.floor(x+radius))
    if out.min < 0: out.min  = 0
    if out.max > maximum: out.max  = maximum

    return out

def gridMinMax1_1(x, maximum, radius):

    out = collections.namedtuple('MinMax',['min','max'])
    out.min = int(np.ceil(x-radius))
    out.max = int(np.floor(x+radius))
    if out.min < 0: out.min  = 0
    if out.max > maximum: out.max  = maximum

    return out

def _poly_sdc_kern_0lobes (r):

    POLY_ORDER = 5
    FIT_LEN = 394
    SPECTRAL_LEN = 25600
    FOV = 63

    x = SPECTRAL_LEN * r / float(FOV)

    poly = np.array ((-1.1469041640943728E-13,
                      8.5313956268885989E-11,
                      1.3282009203652969E-08,
                      -1.7986635886194154E-05,
                      3.4511129626832091E-05,
                      0.99992359966186584))

    out = poly[5]

    for i in range(1,POLY_ORDER+1):
        out += pow(x,i) * poly[POLY_ORDER-i]

    if out < 0:
        out = 0

    return out

def loadKernelTable(length):

    rfp = 0.96960938
    out = np.zeros(length,order='F')

    for i in range(0,length):
        out[i] = _poly_sdc_kern_0lobes(math.sqrt(rfp**2 * float(i) / float(length-1)))

    return out


def demo_interSegmentPhaseCorr (navData, rawData, dataObject):


    navData = np.squeeze(navData)

    dimHigh = navData.shape[2]
    dimZ = dataObject.acqp.NI
    dimSeg = dataObject.method.NSegments
    dimYrSeg = navData.shape[0]

    navData = utils.mirrorOddLines(navData,1)

    navData = np.fft.fftshift(np.fft.ifft2(navData, axes=(1,)), axes=1)



    # navData = np.unwrap(np.angle(navData), axis=1)

    phase = np.zeros(navData.shape, 'f4')

    for _dimHigh in range (dimHigh):

        _dimZ = _dimHigh // dimSeg
        _dimSeg = _dimHigh % dimSeg

        if _dimHigh == _dimSeg:
            pass

        for _dimYrSeg in range(dimYrSeg):
            validPoints = utils.markValidPoints1(navData[_dimYrSeg,:,_dimHigh])
            phase[_dimYrSeg,:,_dimHigh] = np.unwrap(np.angle(navData[_dimYrSeg,:,_dimHigh] * validPoints)) * validPoints


    plt.figure()

    plt.suptitle('Segment to segment phase change as it goes trough the readout', fontsize=24)

    line = 5
    ind_image1 = range(0,10)
    ind_image2 = range(10, 20)

    for i in range(28, 68, 5):
        plt.subplot(1, 2, 1)
        plt.plot((phase[line, i, ind_image1]), label='sample = %s' % i)
        plt.subplot(1, 2, 2)
        plt.plot((phase[line, i, ind_image2]), label='sample = %s' % i)

    plt.subplot(1, 2, 1)
    plt.title("Slice nr. 1", fontsize=24)
    plt.xlabel('Segment nr.', fontsize=24)
    plt.ylabel('Phase [rad]', fontsize=24)

    plt.subplot(1, 2, 2)
    plt.title("Slice nr. 2", fontsize=24)
    plt.xlabel('Segment nr.', fontsize=24)
    plt.ylabel('Phase [rad]', fontsize=24)

    plt.legend(loc = 'best', fontsize=24 )

    plt.figure()
    for i in ind_image2:
        plt.plot(phase[line, :, i], label='segment = %s' % (i % dimSeg))

    plt.xlabel('x', fontsize=24)
    plt.ylabel('Phase [rad]', fontsize=24)
    plt.suptitle('Segment to segment phase change in different stages of readout', fontsize=24)
    plt.legend(loc = 'upper left', fontsize=24)

    plt.show()

    return