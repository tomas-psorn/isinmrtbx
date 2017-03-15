#   Compatibility with python 3.x
from __future__ import print_function
from __future__ import unicode_literals

AUTHOR = 'Tomas Psorn'
CONTACT = 'tomaspsorn@isibrno.cz'

import sys

import numpy as np

from isinmrtbx.datatypes.GeneralDataTypes import ParameterClass
import isinmrtbx.inout.nifti as niftiIO


class NiftiHdr(ParameterClass):
  """
  """
  def __init__(self):
    # Header key substructure
    self.sizeof_hdr = 384
    self.data_type = 10 * ' '
    self.db_name = 18 * ' '
    self.extents = 0
    self.session_error = 0
    self.regular = 'r'
    self.dim_info = 0

    # Image dim
    self.dim = None
    self.intent_p1 = 0
    self.intent_p2 = 0
    self.intent_p3 = 0
    self.intent_code = 0
    self.datatype = None
    self.bitpix = None
    self.slice_start = 0
    self.pixdim = None
    self.vox_offset = 348
    self.scl_slope = None
    self.scl_inter = None
    self.slice_end = 0
    self.slice_code = ' '
    self.xyzt_units = '0'
    self.cal_max = None
    self.cal_min = None
    self.slice_duration = 0.0
    self.toffset = 0
    self.glmax = 0
    self.glmin = 0

    # File history substructure
    self.descrip = 80 * ' '
    self.aux_file = 24 * ' '
    self.sform_code = 0
    self.qform_code = 0
    self.quatern_b = 0
    self.quatern_c = 0
    self.quatern_d = 0
    self.qoffset_x = 0
    self.qoffset_y = 0
    self.qoffset_z = 0
    self.srow_x = np.zeros(4, dtype='float32', order='F')
    self.srow_y = np.zeros(4, dtype='float32', order='F')
    self.srow_z = np.zeros(4, dtype='float32', order='F')
    self.intent_name = 16 * ' '
    self.magic = 'n+1\0'
    pass


class NiftiObject:
    """

    """
    def __init__(self, **kwargs):

        self.hdr = None
        self.img = None

        self.path = None
        self.bruker = None
        self.data = None

        if 'path' in kwargs:
            self.path = kwargs['path']
            niftiIO.readNiftiHdr(self.path,self)
            if self.path [-4:] == '.hdr':
                path = self.path.replace('.hdr' ,'.img')
                niftiIO.readNiftiImg(path,self)

        elif 'bruker' in kwargs:
            self.bruker = kwargs['bruker']
            niftiIO.bruker2nii(self.bruker,self)

        elif 'hdr' and 'data' in kwargs:
            self.hdr = kwargs['hdr']
            self.img = kwargs['data']
        else:
            print ("You've created an empty nifiti object")


    def write(self, path):

      if path[-4:] != '.nii':
        print("Invalid name of nifti file: ", path)
        sys.exit(1)

      fid = open(path, "wb")

      niftiIO.writeNiftiHdr(fid, NiftiObject)
      niftiIO.writeNiftiImg(fid, NiftiObject)

      fid.close()


