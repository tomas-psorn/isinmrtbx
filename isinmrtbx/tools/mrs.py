import numpy as np
#from scipy import fft, arange
import os
from math import sqrt, exp
import matplotlib
import matplotlib.pyplot as plt

"""
Based on onlineSpectraDisp.py by Ross Wiliams
"""

AUTHOR = 'Tomas Psorn'
CONTACT = 'tomaspsorn@isibrno.cz'


def removeGroupDelay(data):
    """

    """

    return data[np.argmax(np.abs(data)):]

def phaseCorrect (dataIn):
    """

    """
    return dataIn * np.exp(-1 * np.complex(0,np.angle(dataIn[0])))