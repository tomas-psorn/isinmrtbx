import numpy as np
#from scipy import fft, arange
import os
from math import sqrt, exp
import matplotlib
import matplotlib.pyplot as plt

AUTHOR = 'Ross Wiliams'
CONTACT = 'tomaspsorn@isibrno.cz'


def rawDataDict(rawdataobject, searchTerms):
    """
    This is a function to take a file (acqp) and a list of seachTerms and return a dictionary of the values
    """
    fdata = open(rawdataobject, 'r').read()
    tempDict = {}
    comment_start_str = '##$'
    lines = fdata.split('\n')
    desired_lines = []
    track = 0
    for line in lines:
        track += 1
        for item in searchTerms:
            if line.startswith(comment_start_str) and item in line:
                if '(' not in line:
                    desired_lines.append(tuple(line.replace("##$", "").split('=')))
                elif 'Yes' not in lines[track]:
                    tempDict[item] = map(int, lines[track].split(' '))
                else:
                    tempDict[item] = lines[track].split(' ')
    tempDict.update(dict(desired_lines))
    return tempDict


def readBrukerPartialFidFile(location, partialFileName = "rawdata.job0"):
    """
    function to read partial spectroscopic data adapted for Tomas Psorn
    Parameters
    ----------
    path: This selects the path to the file to read
    partialFileName: This is set to default to "rawdata.job0" however can be changed dependant upon needs.

    Returns
    -------
    returns an ndarray of 32bit complex numbers
    -------
    """

    acqp = location + 'acqp'
    pm = ['GO_raw_data_format', 'BYTORDA', 'NI', 'NR', 'ACQ_size', 'GO_data_save', 'GO_block_size', 'AQ_mod',
          'ACQ_experiment_mode']
    pmDict = rawDataDict(acqp, pm)
    pathToTemp = location + partialFileName
    #print (pmDict)
    foo = int(list(pmDict.get('ACQ_size'))[0])

    NI = int(pmDict.get('NI'))
    NR = int(pmDict.get('NR'))
    NA = 5 #This is hardcoded as we need to change the dictionary reader from acessing acqp to accessing the method file

    numDataHighDim = 1
    numSelectedReceivers = 1

    #   get data type and number of bits

    format = np.dtype('i4').newbyteorder('<')
    bits = 32
    blockSize = np.ceil(float(foo) * numSelectedReceivers)
    endFileSize = blockSize * numDataHighDim * NI * NR * NA
    print (endFileSize)
    fidFile = open(pathToTemp, "rb")
    currentFidData = np.fromfile(fidFile, dtype=format, count=-1)
    fidFile.close()
    start = 0
    end = len(currentFidData)
    fidData = np.zeros(endFileSize, dtype='i4')
    fidData[start:end] = currentFidData

    if len(fidData) != blockSize * numDataHighDim * NI * NR * NA:
        print('Missmatch')

    fidData = np.reshape(fidData, (blockSize, numDataHighDim * NI * NR * NA), order='F')

    if blockSize != endFileSize * numSelectedReceivers:
        fidData = np.transpose(fidData, (1, 0))
        fidData = fidData[:, :(endFileSize * numSelectedReceivers)]
        fidData = np.reshape(fidData, (numDataHighDim * NI * NR * NA, foo, numSelectedReceivers), order='F')
        fidData = np.transpose(fidData, (2, 1, 0))
    else:
        fidData = np.reshape(fidData, (foo, numSelectedReceivers, numDataHighDim * NI * NR * NA), order='F')
        fidData = np.transpose(fidData, (1, 0, 2))

    fidData = fidData[:, 0::2, :] + 1j * fidData[:, 1::2, :]
    return fidData[0]


def bruker_getSelectedReceivers(path):
    """
    pvtools original function to determine number of channels used for acquisition

    Parameters
    ----------
    rawdataobject

    Returns
    -------
    Number of channels
    """
    numSelectedReceivers = 0
    acqpFile = path + 'acqp'
    pm = ['ACQ_experiment_mode', 'GO_ReceiverSelect']
    rawdataobject = rawDataDict(acqpFile, pm)
    if rawdataobject.get('ACQ_experiment_mode') == 'ParallelExperiment':
        if hasattr(rawdataobject.acqp, 'GO_ReceiverSelect'):
            if rawdataobject.acqp.GO_ReceiverSelect[0].isalpha():
                numSelectedReceivers = 0
                for channel in rawdataobject.acqp.GO_ReceiverSelect:
                    if channel == 'Yes':
                        numSelectedReceivers += 1
        elif hasattr(rawdataobject.acqp, 'ACQ_ReceiverSelect'):
            if rawdataobject.acqp.ACQ_ReceiverSelect[0].isalpha():
                numSelectedReceivers = 0
                for channel in rawdataobject.acqp.ACQ_ReceiverSelect:
                    if channel == 'Yes':
                        numSelectedReceivers += 1
        else:
            print('Information about number of receivers is unknown, check your acqp file')
    else:
        numSelectedReceivers = 1
    return numSelectedReceivers


def removeBrukerLeadingZeroes(data):
    """
    function that strips leading 0s from Bruker data
    """

    total = 0
    for n in range(len(data)):
        point = data[0,n,0]
        if point == 0j:
            total += 1
        else:
            break
    return data[:,total:,:]


def plotter(data, Fs=400, phi=0):
    """
    This is a function that will take the data and return a plot of both the time domain and the frequency domain signals

    """
    y = removeBrukerLeadingZeroes(data)
    y = data.flatten()
    x = range(len(y))  # length of the signal
    real, imag = removeGroupDelay(data)
    real = phaseCorrect(real, phi)
    imag = phaseCorrect(imag, phi)
    fftdata = [z.real * exp(1j * phi) for z in np.fft.fft(y)][::-1]
    x1 = np.array(range(len(real)))
    x2 = np.array(range(len(imag)))
    plt.figure()
    plt.subplot(3, 1, 1)

    plt.plot(x, y)
    plt.xlabel('Time Domain')

    plt.subplot(3, 1, 2)
    plt.plot(x1, real)
    plt.xlabel('Real Domain')

    plt.subplot(3, 1, 3)
    plt.plot(x2, imag)
    plt.xlabel('Imaginary data')

    plt.savefig('displayFig.png')

    plt.show()

    print ('finished')
    return real


def removeGroupDelay(data):
    """
    This is a function that takes NON fourier transform data as an input and corrects the group delay artefact.

    Returns
    ----------------------------
    two nparrays one with the group delay corrected FOURIER TRANSFORMED real data and one with imaginary
    ----------------------------
    """
    """
    data = data.flatten()
    reals = [x.real for x in data]
    imags = [x.imag for x in data]
    maxReal = np.argmax(reals)
    maxImag = np.argmax(imags)
    realsGDC = reals[maxReal:]  # group delay corrected real data
    imagsGDC = imags[maxImag:]  # gdc corrected imag data
    realsGDCFFT = np.fft.rfft(realsGDC)  # realfft the reals
    imagsGDCFFT = np.fft.rfft(imagsGDC)  # realfft the imags

    return realsGDCFFT, imagsGDCFFT

    """
    return data[np.argmax(np.abs(data)):]


def phaseCorrect ( data, phi=0):
    """

    :param data: input data to be phase corrected (ndarray)
    :param phi: angle of phase correction
    :return: phase corrected data
    """
    return [x*exp(1j*phi) for x in data]