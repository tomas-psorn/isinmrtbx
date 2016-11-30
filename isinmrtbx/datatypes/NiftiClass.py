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

    def __init__(self, **kwargs):

        self.hdr = ParameterClass()
        self.img = NiftiImage()

        path = None
        bruker = None
        data = None

        try:
            path = kwargs['path']
        except:
            pass

        try:
            bruker = kwargs['bruker']
        except:
            pass

        try:
            data = kwargs['data']
        except:
            pass

        if path:
            niftiIO.readNiftiHdr(path,self)
            if path [-4:] == '.hdr':
                path = path.replace('.hdr' ,'.img')
                niftiIO.readNiftiImg(path,self)

        elif bruker:
            niftiIO.bruker2nii(bruker,self)

        elif data:
            niftiIO.data2nii(data)

        else:
            print ("You've created an empty nifiti object")

