from __future__ import print_function
from __future__ import unicode_literals

import sys

def codeXchange(querry,querryType):
    """
    You have a code of, for instance, a Bruker data type,
    but you could use its python equivalent. Say no more!

    Parameters
    ----------
    querry
    querryType

    Returns
    -------
    -
    """
    querry = str(querry)

    if querryType == 'datatype':
        codes = {
            'int8':            2,
            'int16':           4,       #DT_INT16 or DT_SIGNED_SHORT
            'int32':           8,
            'float32':        16,       #DT_FLOAT32
            'float64':        64,
            'complex128':   1792,       #DT_COMPLEX128
        }
        code = codes.get(querry)

    elif querryType == 'xform_code':
        codes = {
            'NIFTI_XFORM_UNKNOWN':      0,
            'NIFTI_XFORM_SCANNER_ANAT': 1,
            'NIFTI_XFORM_ALIGNED_ANAT': 2,
            'NIFTI_XFORM_TALAIRACH':    3,
            'NIFTI_XFORM_MNI_152':      4
        }
        code = codes.get(querry)

    elif querryType == 'xyzt_units':
        codes = {
            'NIFTI_UNITS_UNKNOWN':  0,
            'NIFTI_UNITS_METER':    1,
            'NIFTI_UNITS_MM':       2,
            'NIFTI_UNITS_MICRON':   3,
            'NIFTI_UNITS_SEC':      8,
            'NIFTI_UNITS_MSEC':     16,
            'NIFTI_UNITS_USEC':     24,
            'NIFTI_UNITS_HZ':       32,
            'NIFTI_UNITS_PPM':      40,
            'NIFTI_UNITS_RADS':     48,
        }
        code = codes.get(querry)

    elif querryType == 'bruker2nifti_units':
        codes = {
            'NIFTI_UNITS_MM': 2
        }
        code = codes.get(querry)


    elif querryType == 'niftiCode2pythonDtype':
        codes = {
        4:      'i2',
        8:      'i4',
        16:     'f4',
        }
        code = codes.get(querry)

    elif querryType == 'pythonDtype2niftiCode':
        codes = {
        'i2':      4,
        'i4':      8,
        'f4':     16,
        }
        code = codes.get(querry)


    elif querryType == 'pyDtype2pyUnpack':
        codes = {
        'i2':     'H',
        'i4':     'i',
        'f4':     'f',
        'f8':     'd',
        }
        code = codes.get(querry)

    elif querryType == 'brukerLessAcqp':
        code = ['PULPROG', 'ACQ_size','NI', 'NA','NR','ACQ_grad_str_Y']

    elif querryType == 'brukerLessMethod':
        code = ['EchoTime', 'NSegments','OMNavigators','NAverages', 'NRepetitions','PVM_SignalType',
                'PVM_EncZfRead', 'PVM_EncPpiAccel1', 'PVM_EncPftAccel1', 'PVM_EncZfAccel1',
                'PVM_EncMatrix', 'PVM_EncSteps1', 'PVM_EpiEchoPosition', 'PVM_Matrix']

    else:
        print("Unknown querry type: ", querryType)
        sys.exit(1)

    if code != None:
        return code
    else:
        print("codeXchange ERROR: Code: ", querry, " of querry type: ", querryType, "not found")
        sys.exit(1)
