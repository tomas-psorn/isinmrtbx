from __future__ import print_function

import sys
import operator


import numpy as np
import numpy.matlib as npml
import matplotlib.pyplot as plt

import isinmrtbx.tools.mrs as mrs
from isinmrtbx.inout.jcampdx import *

"""
Temp imports
"""
import matplotlib.pyplot as plt



AUTHOR = 'Tomas Psorn'
CONTACT = 'psorn@isibrno.cz'

def readBrukerParamFile(path, dataobject):
    """
    Function to read Bruker parameter file in JCAMPDX format. It stores all attributes and
    their values found in a file into particular object parameter.

    Parameters
    ----------
    path (string)
    imagedataobject

    Returns
    -------
    -
    """

    fileType = path.split('/')
    fileType = fileType[len(fileType)-1]

    with open(path) as paramFile:
        paramData = paramFile.read()

    paramDataLines, comments = strip_comments(paramData)

    ldr_sep, ldr_usr_sep, ldr_dict_sep = '##', '$', '='

    ldrs = [ldr for ldr in paramDataLines.split(ldr_sep) if ldr]

    for ldr in ldrs:
        try:
            ldr_name, ldr_val = ldr.split(ldr_dict_sep)
        except:
            print('Parameter file line is immposible to process, program continues to run')
            print(ldr_dict_sep)
            continue
        ldr_name=ldr_name[1:]
        ldr_val = parse_record(ldr_val)

        command = 'setattr(''dataobject'+'.' + fileType + ', ldr_name, ldr_val)'
        eval(command)

    return

def readBrukerParamFile2(path, dataobject):
    """
    Function to read Bruker parameter file in JCAMPDX format. It stores all attributes and
    their values found in a file into particular object parameter.

    Parameters
    ----------
    path (string)
    imagedataobject

    Returns
    -------
    -
    """

    fileType = path.split('/')
    fileType = fileType[len(fileType)-1]

    with open(path) as paramFile:
        paramData = paramFile.read()

    ldr_sep, ldr_usr_sep, ldr_dict_sep = '##', '$', '='

    ldrs = [ldr for ldr in paramData.split(ldr_sep) if ldr]

    for ldr in ldrs:
        try:
            ldr_name, ldr_val = ldr.split(ldr_dict_sep)
        except:
            print('Parameter file line is immposible to process, program continues to run')
            print(ldr_dict_sep)
            continue

        ldr_val = parse_record(ldr_val)


        setattr(eval('dataobject.' + fileType + '2'),ldr_name,ldr_val)

    return

def readBrukerTrajFile(path,rawdataobject):
    """
    Function to read trajectory data. It saves trajectory data into rawdataobject.traj.data.

    Parameters
    ----------
    path
    rawdataobject

    Returns
    -------
    -
    """
    dimensions = 2
    projections_nr = rawdataobject.method.PVM_TrajIntAll
    samples_nr = rawdataobject.method.PVM_TrajSamples

    with open(path,"rb") as trajFile:
        rawTrajData = np.fromfile(trajFile,dtype='float64',count=-1)

    trajData = np.zeros((samples_nr, projections_nr, dimensions), dtype=rawTrajData.dtype)
    # x coordinates
    trajData [:,:,0] = np.reshape(rawTrajData[::2],(samples_nr, projections_nr),order='F')
     # y coordinates
    trajData [:,:,1] = np.reshape(rawTrajData[1::2],(samples_nr, projections_nr),order='F')

    return trajData

def readBrukerFidFile(path, dataObject):
    """
    Pvtools original algorithm to read and reshape fid data is used. It lacks some of the advanced functionality
    of the original code.

    The reading process is encapsulated so that it can be used for online data analysis.

    Parameters
    ----------
    path - path  to fid file
    dataObject - bruker scan object

    Returns
    -------
    dataOut - fid data in folloving shape [dimCh, dimAcq0. dimHigh]
    """

    # test the presence of all required parameters, get read essential parameters
    format, bits, dimBlock, dimZ, dimR, dimAcq0, dimAcqHigh, dimCh, dimA = fidBasicInfo(dataObject)

    # reading part
    with open(path,"rb") as fidFile:
        fidData = fidBasicRead(fidFile, format)

    # pv-tools like reshape
    dataOut = fidBasicReshape(fidData, dimBlock, dimZ, dimR, dimAcq0, dimAcqHigh, dimCh, dimA)

    return dataOut


def fidBasicInfo(dataObject):

    minCondition = ('GO_raw_data_format','BYTORDA','NI','NR','ACQ_size','GO_data_save','GO_block_size', 'AQ_mod');
    all_here = bruker_requires(dataObject, minCondition, 'acqp')
    if not all_here:
        print("ERROR: visu_pars file does not provide enough info")
        sys.exit(1)

    dimZ = dataObject.acqp.NI
    dimR = dataObject.acqp.NR
    dimAcqHigh = int(np.prod(dataObject.acqp.ACQ_size[1:]))
    dimAcq0 = dataObject.acqp.ACQ_size[0]
    dimCh = dataObject.getSelectedReceivers()
    acqp_BYTORDA = dataObject.acqp.BYTORDA
    acqp_GO_raw_data_format = dataObject.acqp.GO_raw_data_format
    acqp_GO_block_size = dataObject.acqp.GO_block_size
    method_Method = dataObject.method.Method
    # because of reading jobfiles
    if method_Method == 'PRESS':
        dimA = dataObject.method.PVM_NAverages
    else:
        dimA = 1

    #   get data type and number of bits
    if acqp_GO_raw_data_format == 'GO_32BIT_SGN_INT' and acqp_BYTORDA == 'little':
        format = np.dtype('i4').newbyteorder('<')
        bits = 32
    elif acqp_GO_raw_data_format == 'GO_16BIT_SGN_INT' and acqp_BYTORDA == 'little':
        format = np.dtype('i').newbyteorder('<')
        bits = 16
    elif acqp_GO_raw_data_format == 'GO_32BIT_FLOAT' and acqp_BYTORDA == 'little':
        format = np.dtype('f4').newbyteorder('<')
        bits = 32
    elif acqp_GO_raw_data_format == 'GO_32BIT_SGN_INT' and acqp_BYTORDA == 'big':
        format = np.dtype('i4').newbyteorder('>')
        bits = 32
    elif acqp_GO_raw_data_format == 'GO_16BIT_SGN_INT' and acqp_BYTORDA == 'big':
        format = np.dtype('i').newbyteorder('>')
        bits = 16
    elif acqp_GO_raw_data_format == 'GO_32BIT_FLOAT' and acqp_BYTORDA == 'big':
        format = np.dtype('f4').newbyteorder('>')
        bits = 32
    else:
        format = np.dtype('i4').newbyteorder('<')
        print('Data format not specified correctly, set to int32, little endian')
        bits = 32

    if acqp_GO_block_size == 'Standard_KBlock_Format':
        dimBlock = int(np.ceil(float(dimAcq0)*dimCh*(bits/8)/1024)*1024/(bits/8))
    else:
        dimBlock = dimAcq0 * dimCh

    return format, bits, dimBlock, dimZ, dimR, dimAcq0, dimAcqHigh, dimCh, dimA


def fidBasicReshape(fidData, dimBlock, dimZ, dimR, dimAcq0, dimAcqHigh, dimCh, dimA):

    if len(fidData) != dimBlock * dimAcqHigh * dimZ * dimR * dimA:
        print('Missmatch')

    fidData = np.reshape(fidData, (dimBlock, dimAcqHigh * dimZ * dimR));

    if dimBlock != dimAcq0 * dimCh:
        fidData = np.transpose(fidData, (1, 0))
        fidData = fidData[:, :(dimAcq0 * dimCh)]
        fidData = np.reshape(fidData, (dimAcqHigh * dimZ * dimR * dimA, dimAcq0, dimCh), order='F')
        fidData = np.transpose(fidData, (2, 1, 0))
    else:
        fidData = np.reshape(fidData, (dimAcq0, dimCh, dimAcqHigh * dimZ * dimR * dimA), order='F')
        fidData = np.transpose(fidData, (1, 0, 2))

    fidDataOut = fidData[:, 0::2, :] + 1j * fidData[:, 1::2, :]

    return fidDataOut


def fidBasicRead(fileObject, format, startRead=0, dimRead=-1):


    fileObject.seek(startRead)  # find the start position in the file

    fidData = np.fromfile(fileObject, dtype=format, count=dimRead) # read

    if dimRead == -1:
        return fidData
    else:
        endRead = fileObject.tell()
        return fidData, endRead




def readBruker2dseq(path, imagedataobject):
    """
    pvtools original function to read 2dseq data. It lacks some of the advanced functionality
    of the original function. At the end it saves read data into imagedataobject.data2seq.data.

    Parameters
    ----------
    path
    imagedataobject

    Returns
    -------

    """

    #   Test the presence of all required parameters in visu_pars
    minCondition = ('VisuCoreWordType', 'VisuCoreByteOrder', 'VisuCoreSize', 'VisuCoreFrameCount',
            'VisuCoreDataSlope', 'VisuCoreDataOffs','VisuCoreFrameType', 'VisuCoreDim', 'VisuCoreDimDesc')
    all_here = bruker_requires(imagedataobject, minCondition, 'visu_pars')
    if not all_here:
        print("ERROR: visu_pars file does not provide enough info")
        sys.exit(1)

    #   Transform used visu_pars variables to local variables
    VisuCoreWordType = imagedataobject.visu_pars.VisuCoreWordType
    VisuCoreByteOrder = imagedataobject.visu_pars.VisuCoreByteOrder
    VisuCoreSize = imagedataobject.visu_pars.VisuCoreSize
    VisuCoreDataSlope = imagedataobject.visu_pars.VisuCoreDataSlope
    VisuCoreDataOffs = imagedataobject.visu_pars.VisuCoreDataOffs
    VisuCoreFrameType = imagedataobject.visu_pars.VisuCoreFrameType
    VisuCoreDim = imagedataobject.visu_pars.VisuCoreDim
    VisuCoreDimDesc = imagedataobject.visu_pars.VisuCoreDimDesc
    NF = imagedataobject.visu_pars.VisuCoreFrameCount

    #   Get an idea about higher dim meaning
    if any('VisuFGOrderDesc' in s for s in dir(imagedataobject.visu_pars)) and \
        any('VisuFGOrderDescDim' in s for s in dir(imagedataobject.visu_pars)):
        VisuFGOrderDesc = imagedataobject.visu_pars.VisuFGOrderDesc
        VisuFGOrderDescDim = imagedataobject.visu_pars.VisuFGOrderDescDim

        for i in range (0,VisuFGOrderDescDim):
            if i is 0:
                dim5_n = np.zeros(1, dtype='i4', order='f')
                dim5_n[0] = VisuFGOrderDesc[0]
            else:
                dim5_n = np.append(dim5_n,VisuFGOrderDesc[5*i])

    #   get data type and number of bits
    if VisuCoreWordType == '_32BIT_SGN_INT':
        format = np.dtype('int32')
    elif VisuCoreWordType == '_16BIT_SGN_INT':
        format = np.dtype('int16')
    elif VisuCoreWordType == '_32BIT_FLOAT':
        format = np.dtype('float32')
    elif VisuCoreWordType == '_8BIT_USGN_INT':
        format = np.dtype('uint8')
    else:
        print('Data format not specified correctly!')

    if VisuCoreByteOrder == 'littleEndian':
        format = format.newbyteorder('L')
    elif VisuCoreWordType == 'bigEndian':
        format = format.newbyteorder('B')
    else:
        print('Byte order not specified correctly!')

    #   Read 2seq file ADD catch
    with open(path,"rb") as twodseqFile:
        twodseqData = np.fromfile(twodseqFile, dtype=format, count=-1)

    twodseqData = twodseqData.astype('f4') # so that it can be multiplyed by float slope
    twodseqData = np.reshape(twodseqData,(VisuCoreSize[0],-1),order='F')

    # Mapping data
    # to get the real values, we have to multiply with the slope-factor and
    # add the offset to the data
    slope = npml.repmat(VisuCoreDataSlope,1,np.prod(VisuCoreSize))
    slope = np.reshape(slope, twodseqData.shape, order='F')
    offset = npml.repmat(VisuCoreDataOffs,1,np.prod(VisuCoreSize))
    offset = np.reshape(offset, twodseqData.shape, order='F')

    twodseqData *= slope.astype('f4')
    twodseqData += offset.astype('f4')

    #   Pad VisuCoresize
    if VisuCoreDim is 1:    VisuCoreSize = np.append(VisuCoreSize,[1, 1, 1])
    elif VisuCoreDim is 2:  VisuCoreSize = np.append(VisuCoreSize,[1, 1])
    elif VisuCoreDim is 3:  VisuCoreSize = np.append(VisuCoreSize, 1)
    else:                   pass


    if 'dim5_n' in locals():
        twodseqData = np.reshape(twodseqData,np.append(VisuCoreSize, dim5_n),order='F')
    else:
        twodseqData = np.reshape(twodseqData,np.append(VisuCoreSize, NF),order='F')

    twodseqData  = np.transpose(twodseqData, (1,0,2,3,4)) #todo this is rather heuristic, but seems to work

    return twodseqData


def bruker_requires(dataobject,minCondition, fileType ):
    """
    pvtools original function to control the presence of an essential parameters
    in a dataobject's parameter defined by fileType

    Parameters
    ----------
    dataobject
    minCondition (tuple of strings)
    fileType (string)

    Returns
    -------
    all_here (bool)
    """
    all_here = True
    for conditionElement in minCondition:
        condition = 'dataobject'+'.'+fileType+'.'+ conditionElement
        try:
            eval(condition)
        except AttributeError:
            print('ERROR: ', fileType,' file does not contain essential parameter: ',conditionElement)
            all_here = False
    return all_here


def fidHandle_UTE(rawdataobject):
    minCondition = ('NI',)
    all_here = bruker_requires(rawdataobject, minCondition, 'acqp')
    if not all_here:
        print("ERROR: acqp file does not provide enough info")
        sys.exit(1)
    slices_nr = rawdataobject.acqp.NI
    return

def fidHandle_FAIR_RARE(rawdataobject):
    """
    Description
    -----------
    Simple function to reshape fid data using prior knowledge about data storage.
    For FAIR_RARE sequence in this case. It is just a draft version so far, some
    improvements are to be made.

    Parameters
    ----------
    rawdataobject

    Returns
    -------
    -
    """

    minCondition = ('NR','ACQ_rare_factor')
    paramFile = 'acqp'
    all_here = bruker_requires(rawdataobject, minCondition, paramFile)
    if not all_here:
        print("ERROR: ", paramFile, " file does not provide enough info")
        sys.exit(1)

    # dataIn = rawdataobject.fid.data

    dataIn = rawdataobject.fid.data
    dataIn_type = dataIn.dtype  # so the data format is preserved

    #   Get dataOut dimension parameters, create dataOut np.array
    channels_nr = dataIn.shape[0]
    repetitions_nr = rawdataobject.acqp.NR     # just for the readability
    images_nr = rawdataobject.acqp.NI
    views_nr = rawdataobject.acqp.ACQ_rare_factor
    samples_nr = dataIn.shape[1]
    dataOut = np.zeros((channels_nr,repetitions_nr,views_nr,samples_nr),dataIn_type)
    if repetitions_nr > 1:
        for _repetition in range(0,repetitions_nr):
            dataOut[:,_repetition,:,:] = dataIn[:,:,_repetition*views_nr:(_repetition+1)*views_nr].transpose((0,2,1))

    rawdataobject.fid.data = dataOut    # replace
    return

def fidHandle_Epi(rawdataobject):
    """
    Ii is supposed to rearrange all fid data acquired by EPI based sequence,
    such as
    Parameters
    ----------
    rawdataobject

    Returns
    -------

    """
    # dataIn = rawdataobject.fid.data #todo delete once the new structures are implemented
    dataIn = rawdataobject.fid

    minCondition = ('NR','NI')
    all_here = bruker_requires(rawdataobject, minCondition, 'acqp')
    if not all_here:
        print("ERROR: acqp file does not provide enough info")
        sys.exit(1)

    minCondition = ('PVM_EncMatrix',)
    all_here = bruker_requires(rawdataobject, minCondition, 'method')
    if not all_here:
        print("ERROR: acqp file does not provide enough info")
        sys.exit(1)

    segments_nr = rawdataobject.method.NSegments

    try:
        lines_nr = rawdataobject.method.PVM_EpiMatrix[1]
        samples_nr = rawdataobject.method.PVM_EpiMatrix[0]
    except:
        lines_nr = rawdataobject.method.PVM_EncMatrix[1]
        samples_nr = rawdataobject.method.PVM_EncMatrix[0]

    linesSegment_nr = lines_nr // segments_nr
    numDataHighDim = dataIn.shape[2]
    receivers_nr = rawdataobject.getSelectedReceivers()

    # dataIn is supposed to be in the shape of [channels_nr, acq size[0], numhigh dim]
    dataIn = dataIn[:, -linesSegment_nr*samples_nr:, :] # to get rid of navigator data, if present

    dataOut = np.reshape(dataIn,(receivers_nr,linesSegment_nr, samples_nr, numDataHighDim),order='C') #todo why order='C'?

    del rawdataobject.fid
    setattr(rawdataobject,'fid',dataOut)

    return

def fidHandle_Flash(scan):

    dataIn = scan.fid

    dataDim = scan.method.PVM_Matrix
    NI = scan.acqp.NI
    NR = scan.acqp.NR
    numHighDim = NI*NR

    dataOut = np.reshape(dataIn,(dataDim[0], dataDim[1], numHighDim), order='C')

    scan.fid = None
    scan.fid = dataOut

    return

def fidHandle_PRESS(scan):

    dataIn = scan.fid
    scan.fid = None

    dims = dataIn.shape # [channels, xDimR, highDim]
    dataOut = np.zeros_like(dataIn)

    for _dim0 in range(dims[0]):
        for _dim2 in range(dims[2]):
            fid = mrs.phaseCorrect(mrs.removeGroupDelay(dataIn[_dim0, :, _dim2]))
            dataOut [_dim0, 0:len(fid) , _dim2] = fid

    scan.fid = dataOut

    return

