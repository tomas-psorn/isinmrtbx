#   Compatibility with python 3.x
from __future__ import print_function
from __future__ import unicode_literals

AUTHOR = 'Tomas Psorn'
CONTACT = 'tomaspsorn@isibrno.cz'

from isinmrtbx.datatypes.GeneralDataTypes import ParameterClass
import isinmrtbx.inout.nifti as niftiIO

class NiftiImage:
    pass

class NiftiObject:
    """

    """
    def __init__(self, **kwargs):

        self.hdr = ParameterClass()
        self.img = NiftiImage()

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

        elif 'data' in kwargs:
            self.data = kwargs['data']
            niftiIO.data2nii(self.data)
        else:
            print ("You've created an empty nifiti object")

