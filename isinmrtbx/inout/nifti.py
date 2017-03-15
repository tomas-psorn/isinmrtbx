from __future__ import print_function

import numpy as np
import struct
import sys


import matplotlib.pyplot as plt

from isinmrtbx.inout.codeXchange import codeXchange

AUTHOR = 'Tomas Psorn'
CONTACT = 'tomaspsorn@isibrno.cz'

def readNiftiImg(path, niftiObject):

    fileType = path[-4:]

    if fileType == '.nii':
        offset = niftiObject.hdr.vox_offset
    elif fileType == '.img':
        offset = 0
    else:
        print("File: ", path, " is not .nii, nor .img")
        sys.exit(1)

    bytes = niftiObject.hdr.bitpix/8
    image_size = np.prod(niftiObject.hdr.dim[1:4])


    dataType = codeXchange(niftiObject.hdr.datatype,'niftiCode2pythonDtype')


    unpackCode = codeXchange(dataType,'pyDtype2pyUnpack')

    d1 = niftiObject.hdr.dim[1]
    d2 = niftiObject.hdr.dim[2]
    d3 = niftiObject.hdr.dim[3]
    d4 = niftiObject.hdr.dim[4]
    d5 = niftiObject.hdr.dim[5]
    d6 = niftiObject.hdr.dim[6]
    d7 = niftiObject.hdr.dim[7]

    if d1<1: d1 = 1
    if d2<1: d2 = 1
    if d3<1: d3 = 1
    if d4<1: d4 = 1
    if d5<1: d5 = 1
    if d6<1: d6 = 1
    if d7<1: d7 = 1

    fid = open(path,"rb")

    img_idx = range(0,d4)
    dim5_idx = range(0,d5)
    dim6_idx = range(0,d6)
    dim7_idx = range(0,d7)

    image = np.zeros((image_size, len(img_idx)* len(dim5_idx) * len(dim6_idx) * len(dim7_idx)),dtype=dataType)
    currentIndex = 0

    for i7 in dim7_idx:
        for i6 in dim6_idx:
            for i5 in dim5_idx:
                for t in img_idx:
                    inFilePointer = np.ravel_multi_index((0, 0, 0, img_idx[t], dim5_idx[i5],dim6_idx[i6],dim7_idx[i7]), dims=(d1, d2, d3, d4, d5, d6, d7), order='F')
                    inFilePointer *= bytes
                    fid.seek(offset + inFilePointer,0)
                    fid.seek(offset + inFilePointer,0)
                    image[:,currentIndex] = struct.unpack(unpackCode*image_size, fid.read(image_size * bytes))
                    currentIndex +=1

    fid.close()


    image = np.reshape(image,(d1,d2,d3,d4,d5,d6,d7),order='F')

    niftiObject.img = image

    return


def readNiftiHdr(path, niftiObject):

    fid = open(path,"rb")

    if fid.tell() != 0:
        fid.seek(0)

    # Header key substructure
    niftiObject.hdr.sizeof_hdr =    struct.unpack('i', fid.read(4))[0]
    niftiObject.hdr.data_type =     fid.read(10)
    niftiObject.hdr.db_name =       fid.read(18)
    niftiObject.hdr.extents =       struct.unpack('i', fid.read(4))[0]
    niftiObject.hdr.session_error = struct.unpack('h', fid.read(2))[0]
    niftiObject.hdr.regular =       struct.unpack('c', fid.read(1))[0]
    niftiObject.hdr.dim_info =      struct.unpack('c', fid.read(1))[0]

    # Image dim
    niftiObject.hdr.dim =           struct.unpack('h'*8, fid.read(16))
    niftiObject.hdr.intent_p1 =     struct.unpack('f', fid.read(4))[0]
    niftiObject.hdr.intent_p2 =     struct.unpack('f', fid.read(4))[0]
    niftiObject.hdr.intent_p3 =     struct.unpack('f', fid.read(4))[0]
    niftiObject.hdr.intent_code =   struct.unpack('h', fid.read(2))[0]
    niftiObject.hdr.datatype =      struct.unpack('h', fid.read(2))[0]
    niftiObject.hdr.bitpix =        struct.unpack('h', fid.read(2))[0]
    niftiObject.hdr.slice_start =   struct.unpack('h', fid.read(2))[0]
    niftiObject.hdr.pixdim =           struct.unpack('f'*8, fid.read(32))
    niftiObject.hdr.vox_offset =      struct.unpack('f', fid.read(4))[0]
    niftiObject.hdr.scl_slope =        struct.unpack('f', fid.read(4))[0]
    niftiObject.hdr.scl_inter =   struct.unpack('f', fid.read(4))[0]
    niftiObject.hdr.slice_end =        struct.unpack('h', fid.read(2))[0]
    niftiObject.hdr.slice_code =       struct.unpack('c', fid.read(1))[0]
    niftiObject.hdr.xyzt_units =      struct.unpack('c', fid.read(1))[0]
    niftiObject.hdr.cal_max =     struct.unpack('f', fid.read(4))[0]
    niftiObject.hdr.cal_min =     struct.unpack('f', fid.read(4))[0]
    niftiObject.hdr.slice_duration =     struct.unpack('f', fid.read(4))[0]
    niftiObject.hdr.toffset =     struct.unpack('f', fid.read(4))[0]
    niftiObject.hdr.glmax =     struct.unpack('i', fid.read(4))[0]
    niftiObject.hdr.glmin =     struct.unpack('i', fid.read(4))[0]

    # File history substructure
    niftiObject.hdr.descrip = fid.read(80)
    niftiObject.hdr.aux_file = fid.read(24)
    niftiObject.hdr.sform_code = struct.unpack('h', fid.read(2))[0]
    niftiObject.hdr.qform_code = struct.unpack('h', fid.read(2))[0]
    niftiObject.hdr.quatern_b = struct.unpack('f', fid.read(4))[0]
    niftiObject.hdr.quatern_c = struct.unpack('f', fid.read(4))[0]
    niftiObject.hdr.quatern_d = struct.unpack('f', fid.read(4))[0]
    niftiObject.hdr.qoffset_x = struct.unpack('f', fid.read(4))[0]
    niftiObject.hdr.qoffset_y = struct.unpack('f', fid.read(4))[0]
    niftiObject.hdr.qoffset_z = struct.unpack('f', fid.read(4))[0]
    niftiObject.hdr.srow_x = struct.unpack('f'*4, fid.read(16))
    niftiObject.hdr.srow_y = struct.unpack('f'*4, fid.read(16))
    niftiObject.hdr.srow_z = struct.unpack('f'*4, fid.read(16))
    niftiObject.hdr.intent_name = fid.read(16)
    niftiObject.hdr.magic = fid.read(4)

    fid.close()

    return


def writeNiftiHdr(fid, niftiObject):

    # Header key substructure
    fid.write(struct.pack('i', niftiObject.hdr.sizeof_hdr))
    fid.write(struct.pack('10s', niftiObject.hdr.data_type.encode('UTF-8')))
    fid.write(struct.pack('18s', niftiObject.hdr.db_name.encode('UTF-8')))
    fid.write(struct.pack('i', niftiObject.hdr.extents))
    fid.write(struct.pack('h', niftiObject.hdr.session_error))
    fid.write(struct.pack('s', niftiObject.hdr.regular.encode('UTF-8')))
    fid.write(struct.pack('s', niftiObject.hdr.dim_info.encode('UTF-8')))

    # Image dim
    fid.write(struct.pack('8h', *niftiObject.hdr.dim))
    fid.write(struct.pack('f', niftiObject.hdr.intent_p1))
    fid.write(struct.pack('f', niftiObject.hdr.intent_p2))
    fid.write(struct.pack('f', niftiObject.hdr.intent_p3))
    fid.write(struct.pack('h', niftiObject.hdr.intent_code))
    fid.write(struct.pack('h', niftiObject.hdr.datatype))
    fid.write(struct.pack('h', niftiObject.hdr.bitpix))
    fid.write(struct.pack('h', niftiObject.hdr.slice_start))
    fid.write(struct.pack('8f', *niftiObject.hdr.pixdim))
    fid.write(struct.pack('f', niftiObject.hdr.vox_offset))
    fid.write(struct.pack('f', niftiObject.hdr.scl_slope))
    fid.write(struct.pack('f', niftiObject.hdr.scl_inter))
    fid.write(struct.pack('h', niftiObject.hdr.slice_end))
    fid.write(struct.pack('c', niftiObject.hdr.slice_code.encode('UTF-8')))
    fid.write(struct.pack('c', niftiObject.hdr.xyzt_units.encode('UTF-8')))
    fid.write(struct.pack('f', niftiObject.hdr.cal_max))
    fid.write(struct.pack('f', niftiObject.hdr.cal_min))
    fid.write(struct.pack('f', niftiObject.hdr.slice_duration))
    fid.write(struct.pack('f', niftiObject.hdr.toffset))
    fid.write(struct.pack('i', niftiObject.hdr.glmax))
    fid.write(struct.pack('i', niftiObject.hdr.glmin))

    # File history substructure
    fid.write(struct.pack('80s', niftiObject.hdr.descrip.encode('UTF-8')))
    fid.write(struct.pack('24s', niftiObject.hdr.aux_file.encode('UTF-8')))
    fid.write(struct.pack('h', niftiObject.hdr.sform_code))
    fid.write(struct.pack('h', niftiObject.hdr.qform_code))
    fid.write(struct.pack('f', niftiObject.hdr.quatern_b))
    fid.write(struct.pack('f', niftiObject.hdr.quatern_c))
    fid.write(struct.pack('f', niftiObject.hdr.quatern_d))
    fid.write(struct.pack('f', niftiObject.hdr.qoffset_x))
    fid.write(struct.pack('f', niftiObject.hdr.qoffset_y))
    fid.write(struct.pack('f', niftiObject.hdr.qoffset_z))
    fid.write(struct.pack('4f', *niftiObject.hdr.srow_x))
    fid.write(struct.pack('4f', *niftiObject.hdr.srow_y))
    fid.write(struct.pack('4f', *niftiObject.hdr.srow_z))
    fid.write(struct.pack('16s', niftiObject.hdr.intent_name.encode('UTF-8')))
    fid.write(struct.pack('4s', niftiObject.hdr.magic.encode('UTF-8')))

    return

def writeNiftiImg(fid, niftiObject):

    skip_bytes = niftiObject.hdr.vox_offset - 348

    if skip_bytes:
        fid.write(struct.pack('B' * int(skip_bytes), *np.zeros(skip_bytes,dtype='i1', order='F')))

    im = np.reshape(niftiObject.img,(1, -1),order='F')
    im.tofile(fid)

    return

def bruker2nii(bruker, nii):

    # todo why?
    # bruker.data2dseq = bruker.data2dseq.astype('i',casting='safe')

    print('here')

    bruker2hdr(bruker, nii)
    bruker2img(bruker, nii)


    return

def bruker2hdr(bruker,nii):

    data_dims = np.asarray(bruker.data2dseq.shape)
    datatype = codeXchange(bruker.data2dseq.dtype, 'datatype')
    fov = bruker.visu_pars.VisuCoreExtent
    matrixSize = bruker.visu_pars.VisuCoreSize
    frameThickness = bruker.visu_pars.VisuCoreFrameThickness
    bitpix = bruker.data2dseq.dtype.itemsize * 8
    dataSlope = bruker.visu_pars.VisuCoreDataSlope[0]
    dataOffset = bruker.visu_pars.VisuCoreDataOffs[0]
    cal_max = np.amax(dataSlope * bruker.data2dseq + dataOffset)
    cal_min = np.amin(dataSlope * bruker.data2dseq + dataOffset)


    data_dims = data_dims[data_dims != 1]
    dims = np.ones(8, dtype='i2', order='F')
    dims[0] = len(data_dims)
    dims[1:1 + len(data_dims)] = data_dims


    voxel_size = np.ones(8, dtype='f4', order='F')
    voxel_size[0] = 0.0
    voxel_size[1] = fov[0]/float(matrixSize[0])
    voxel_size[2] = fov[1] / float(matrixSize[1])
    voxel_size[3] = frameThickness



    # Header key substructure
    nii.hdr.sizeof_hdr      = 384
    nii.hdr.data_type       = 10 * ' '
    nii.hdr.db_name         = 18 * ' '
    nii.hdr.extents         = 0
    nii.hdr.session_error   = 0
    nii.hdr.regular         = 'r'
    nii.hdr.dim_info        = ' '

    # Image dim
    nii.hdr.dim = dims
    nii.hdr.intent_p1 = 0
    nii.hdr.intent_p2 = 0
    nii.hdr.intent_p3 = 0
    nii.hdr.intent_code = 0
    nii.hdr.datatype = datatype
    nii.hdr.bitpix = bitpix
    nii.hdr.slice_start = 0
    nii.hdr.pixdim = voxel_size
    nii.hdr.vox_offset = 348
    nii.hdr.scl_slope = dataSlope
    nii.hdr.scl_inter = dataOffset
    nii.hdr.slice_end = 0
    nii.hdr.slice_code = ' '
    nii.hdr.xyzt_units = '0'
    nii.hdr.cal_max = cal_max
    nii.hdr.cal_min = cal_min
    nii.hdr.slice_duration = 0
    nii.hdr.toffset = 0
    nii.hdr.glmax = 0
    nii.hdr.glmin = 0

    # File history substructure
    nii.hdr.descrip = 80 * ' '
    nii.hdr.aux_file = 24 * ' '
    nii.hdr.sform_code = 0
    nii.hdr.qform_code = 0
    nii.hdr.quatern_b = 0
    nii.hdr.quatern_c = 0
    nii.hdr.quatern_d = 0
    nii.hdr.qoffset_x = 0
    nii.hdr.qoffset_y = 0
    nii.hdr.qoffset_z = 0
    nii.hdr.srow_x = np.zeros(4, dtype='float32', order='F')
    nii.hdr.srow_y = np.zeros(4, dtype='float32', order='F')
    nii.hdr.srow_z = np.zeros(4, dtype='float32', order='F')
    nii.hdr.intent_name = 16 * ' '
    nii.hdr.magic = 'n+1\0'

    return

def bruker2img(bruker, nii):

    dataOut = bruker.data2dseq.data
    nii.img = dataOut

    return

def data2nii(data):
    print ("Function data2nii is not implemented yet")
    sys.exit(1)