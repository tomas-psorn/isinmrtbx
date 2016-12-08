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


    plt.figure()
    plt.plot(np.imag(dataIn[0:20]))


    phi = np.complex(0,np.angle(dataIn[0]))
    ex = np.exp(-1*phi)
    dataOut = dataIn * ex

    plt.figure()
    plt.plot(np.imag(dataOut[0:20]))
    plt.show()

    return dataOut