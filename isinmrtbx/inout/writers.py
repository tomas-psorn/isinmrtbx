from __future__ import print_function
from __future__ import unicode_literals

import sys

import matplotlib.image as mp_img
import matplotlib.pyplot as plt
import matplotlib.cm as mp_cm

import numpy as np
from isinmrtbx.inout.codeXchange import codeXchange

AUTHOR = 'Tomas Psorn'
CONTACT = 'tomaspsorn@isibrno.cz'

"""
def writeComplex2DArrayToNifty(data, dataObject,nameOfNiftyFile):
    if nameOfNiftyFile[-4:] != '.nii':
        nameOfNiftyFile += '.nii'

    image = nib.Nifti1Image(data,None, extra=None)

    fovDim0_cm = 0
    fovDim0_cm = 0
    fovDim1_cm = 0
    fovDim0_px = 0
    fovDim1_px = 0
    NI = 0
    NR = 0

    try:
        fovDim0_cm = dataObject.acqp.ACQ_fov[0]
        fovDim1_cm = dataObject.acqp.ACQ_fov[1]
        fovDim0_px = dataObject.method.PVM_Matrix[0]
        fovDim1_px = dataObject.method.PVM_Matrix[1]
        NI = dataObject.acqp.NI
        NR = dataObject.acqp.NR
    except:
        fovDim0_cm = dataObject.visu_pars.VisuCoreExtent[0] * 0.1
        fovDim1_cm = dataObject.visu_pars.VisuCoreExtent[1] * 0.1
        fovDim0_px = dataObject.visu_pars.VisuCoreSize[0]
        fovDim1_px = dataObject.visu_pars.VisuCoreSize[1]
        NI = dataObject.reco.RecoObjectsPerRepetition
        NR = dataObject.reco.RecoNumRepetitions

    #   datatype
    image.header['data_type'] = data.dtype.name     # just analyze 7.5 stuff?
    image.header['datatype'] = codeXchange(data.dtype.name,'datatype')

    #   bitpix
    image.header['bitpix'] = data.dtype.itemsize*8


    #   pixdim
    qfac = -1.0
    image.header['pixdim'][0] = qfac
    image.header['pixdim'][1] = 10*float(fovDim0_cm) / float(fovDim0_px)
    image.header['pixdim'][2] = 10*float(fovDim1_cm) / float(fovDim1_px)
    image.header['pixdim'][3] = NI
    image.header['pixdim'][4] = NR

    #   xyzt_units
    image.header['xyzt_units']= codeXchange('NIFTI_UNITS_MM','xyzt_units') +  codeXchange('NIFTI_UNITS_SEC','xyzt_units')

    #   description
    lenghtOfDescription = 80
    image.header['descrip'] = ('Data owner: ' + dataObject.acqp.ACQ_operator + ' NMR group ISI Brno' )[0:lenghtOfDescription-1]

    #   qform & sform
    image.header['qform_code'] = codeXchange('NIFTI_XFORM_SCANNER_ANAT','xform_code')
    image.header['sform_code'] = codeXchange('NIFTI_XFORM_SCANNER_ANAT','xform_code')

    nib.save(image,nameOfNiftyFile)
    return

"""

def writeJpeg(image,name):
    mp_img.imsave(name,image)
    return

def writePng(name, image):
    mp_img.imsave(name,image, None, None, mp_cm.Greys_r)
    return

def writeFigure(name, image, formatStr):
    """

Examples of format strings:
============================
formatStr = 'normalize-colorbar-jet'
------------------------------------
    Saves figure with normalized image, using jet colormap,
    including colorbar

formatStr = 'normalize-colorbar-gray'
------------------------------------
    Saves figure with normalized image, using grayscale colormap,
    including colorbar

    """

    if formatStr.find('normalize') != -1:
        image = image/np.amax(image)

    plt.figure()

    if formatStr.find('gray') != -1:
        plt.imshow(image, mp_cm.Greys_r)
    elif formatStr.find('hsv') != -1:
        plt.imshow(image, mp_cm.hsv)
    elif formatStr.find('jet') != -1:
        plt.imshow(image, mp_cm.jet)
    else:
        plt.imshow(image, mp_cm.Greys_r)

    if formatStr.find('colorbar') != -1:
        if formatStr.find('colorbarlabeled') != -1:
            plt.colorbar(ticks=np.linspace(np.amin(image),np.amax(image),10))
        else:
            plt.colorbar()



    plt.savefig(name)

    return


