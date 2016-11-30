"""
Subject
Scan
Reco
"""
from __future__ import print_function
from __future__ import unicode_literals

AUTHOR = 'Tomas Psorn'
CONTACT = 'tomaspsorn@isibrno.cz'

import sys

import numpy as np

from os import listdir
from os import walk

# import xlwt

from isinmrtbx.datatypes.GeneralDataTypes import ParameterClass
from isinmrtbx.inout import readers
from isinmrtbx.tools import recos
from isinmrtbx.inout import codeXchange
from isinmrtbx.datatypes.NiftiClass import NiftiObject
from isinmrtbx.inout import nifti

class Scan():

    def __init__(self,path, readFid=True, readReco=False):

        self.path = None
        self.scanContent = None
        self.isAcqp = None
        self.isMethod = None
        self.isFid = None
        self.isTraj = None
        self.isPdata = None
        self.recoList = None
        self.fid = None
        self.reco = None
        self.isScan = None

        if path[-1] != "/":
            self.path = path + '/'
        else:
            self.path = path

        # Get idea about what is contained in the scan
        self.browseScan()

        if self.isAcqp is False or self.isMethod is False or self.isFid is False:
            self.isScan = False
            return None
        else:
            self.isScan = True

        # Create instances of parameter and data classes
        if self.isAcqp:
            self.acqp = ParameterClass()
            readers.readBrukerParamFile(self.path + 'acqp',self)

        if self.isMethod:
            self.method = ParameterClass()
            readers.readBrukerParamFile(self.path + 'method',self)

        if self.isFid and readFid:
            self.fid = readers.readBrukerFidFile (self.path + 'fid',    self)
            self.methodBasedReshape()

        # if at least one complete reco is present and reading of recos is demanded, read them
        if self.isPdata and readReco:
            for recoDir in self.recoList:
                if not self.reco:
                    self.reco = [Reco(path = self.path + 'pdata/' + recoDir + '/')]
                else:
                    self.reco.append(Reco(path=self.path + 'pdata/' + recoDir + '/'))

        self.checkOverflow()


    def browseScan(self):
        # Get idea about what is contained in a scan
        self.scanContent = listdir(self.path)
        self.isAcqp = 'acqp' in self.scanContent
        self.isMethod = 'method' in self.scanContent
        self.isFid = 'fid' in self.scanContent
        self.isTraj = 'traj' in self.scanContent
        self.isPdata = 'pdata' in self.scanContent

        if self.isPdata:
            scanRecos = listdir(self.path + 'pdata/')
            self.recoList = []
            for file in scanRecos: #todo not the best criteria
                if file.isdigit():
                    self.recoList.append(file)

    def methodBasedReshape(self):
        """
        Description
        -----------
        Call a particular function, according to the pulse program, to reshape fid data.
        If no function for handling of the fid of a given pulse program (sequence)
        is not specified, nothing happens.

        To add a new function for a certain sequence, just add an another elif line.

        Parameters
        ----------
        rawdataobject

        Returns
        -------

        -
        """

        if self.acqp.PULPROG in ['UTE.ppg']:
            readers.fidHandle_UTE(self)
        elif self.acqp.PULPROG in ['FAIR_RARE.ppg']:
            readers.fidHandle_FAIR_RARE(self)
        elif self.acqp.PULPROG in ['DtiEpi.ppg', 'EPI.ppg', 'navigatorEPI_OM.ppg']:
            readers.fidHandle_Epi(self)
        elif self.acqp.PULPROG in ['FLASH.ppg']:
            readers.fidHandleFlash(self)
        else:
            print('Function to reshape fid data of' + self.acqp.PULPROG + ' sequence is not developed yet')

        return

    def reconstruct(self, **kwargs):
        """
        Description
        -----------
        Call a particular function, according to the pulse program, to reconstruct fid data.

        Parameters
        ----------
        rawdataobject

        Returns
        -------

        -
        """
        # clean & init

        if self.acqp.PULPROG in ['UTE.ppg']:
            result = recos.uteReco(self, **kwargs)
        elif self.acqp.PULPROG in ['DtiEpi.ppg', 'EPI.ppg', 'navigatorEPI_OM.ppg']:
            result = recos.epiReco(self, **kwargs)
        elif self.acqp.PULPROG in ['FLASH.ppg']:
            result = recos.flashReco(self, **kwargs)
        else:
            print('Function to reconstruct fid data of' + self.acqp.PULPROG + ' sequence is not developed yet')

        if not self.reco:   # if the list of recos of this scan is empty
            self.reco = [Reco()]    # create an empty reco object
        else:
            self.append(Reco())

        myRecoIndex = len(self.reco) - 1

        # fill the data and meta data about reco to the structure
        self.reco[myRecoIndex].data2dseq = result


        return

    def getSelectedReceivers(self):
        """
        pvtools original function to determine number of channels used for acquisition

        Parameters
        ----------
        rawdataobject

        Returns
        -------
        Number of channels
        """
        if self.acqp.ACQ_experiment_mode == 'ParallelExperiment':
            if hasattr(self.acqp, 'GO_ReceiverSelect'):
                if self.acqp.GO_ReceiverSelect[0].isalpha():
                    numSelectedReceivers = 0
                    for channel in self.acqp.GO_ReceiverSelect:
                        if channel == 'Yes':
                            numSelectedReceivers += 1
            elif hasattr(self.acqp, 'ACQ_ReceiverSelect'):
                if self.acqp.ACQ_ReceiverSelect[0].isalpha():
                    numSelectedReceivers = 0
                    for channel in self.acqp.ACQ_ReceiverSelect:
                        if channel == 'Yes':
                            numSelectedReceivers += 1
            else:
                print('Information about number of receivers is unknown, check your acqp file')
        else:
            numSelectedReceivers = 1
        return numSelectedReceivers


    def getGrappaCoef(self):

        try:
            dataIn = self.method.PVM_EpiGrappaCoefficients
        except:
            print("Scan does not posses PVM_EpiGrappaCoefficients parameter")
            sys.exit(1)

        recievers_nr = readers.bruker_getSelectedReceivers(self)
        slices_nr = self.acqp.NI

        dataOut = np.reshape(dataIn,(slices_nr,recievers_nr,-1),order='F')

        return dataOut

    def checkOverflow(self):

        try:
            overflow = self.acqp.ACQ_adc_overflow
            if 'Yes' in overflow:
                print("There was an overflow in one of the channels")

        except:
            print ("Unable to chceck wehther there was an ACQ overflow")

        pass

    def export(self, attribute, less=False, mode='print'):
        """
        Export parameters from given parameter file.
        Parameters
        ----------
        attribute
        less

        Returns
        -------

        """
        # check whether the attribute is valid
        valid = ['acqp', 'method', 'all']
        if attribute not in valid:
            print('Function:', __name__, '.Scan.export')
            print('Problem: atribute: ', attribute, ' is not avalible, please chose one of [acqp, method]')
            sys.exit(1)


        if attribute == 'all':
            avalible = valid[:-1]
        else:
            avalible = [attribute,]

        namesGlobal = list()
        valuesGlobal = list()

        for _attribute in avalible:
            # get all avalible parameters within the parameter group of the instance
            all = eval('dir(' + 'self' + '.' + str(_attribute) + ')')

            if less:
                # get the list of what is considered less, defined in codeXchange
                if _attribute == 'acqp':
                    listOfLess = codeXchange.codeXchange(None,'brukerLessAcqp')
                elif _attribute == 'method':
                    listOfLess = codeXchange.codeXchange(None, 'brukerLessMethod')
                # make an intersection of what's avalible and what's defined to be the less set
                names = list(set(all).intersection(listOfLess))
            else:
                names = all

            for name in names:
                namesGlobal.append(name)
                valuesGlobal.append(eval('getattr(' + 'self' + '.' + str(_attribute) + ',name)'))

            # print names and values
        if mode == 'print':
            for name in namesGlobal:
               print(name, valuesGlobal[namesGlobal.index(name)])
            return

        elif mode == 'lists':
            return namesGlobal, valuesGlobal
        else:
            print('Mode invalid')
            return



class Reco():

    # Create empty reco
    def __init__(self, path=None):

        self.data2dseq = None
        self.isVisu = None
        self.isReco = None
        self.is2dseq = None
        self.visu_pars = None
        self.reco = None


        # if path was passed
        if path:
            self.path = path
            self.data2dseq = Reco()
            self.browseReco()

            if self.isVisu:
                self.visu_pars = ParameterClass()
                readers.readBrukerParamFile(self.path + 'visu_pars', self)

            if self.reco:
                self.reco = ParameterClass()
                readers.readBrukerParamFile(self.path + 'reco', self)

            if self.is2dseq and self.isReco and self.isVisu:
                self.data2dseq = readers.readBruker2dseq(self.path + '2dseq', self)
                # setattr(self.data2dseq,'data',data)
        else:
            pass # create an empty reco


    def browseReco(self):
        recoContent = listdir(self.path)
        self.isVisu = 'visu_pars' in recoContent
        self.isReco = 'reco' in recoContent
        self.is2dseq = '2dseq' in recoContent


    def methodBasedReshape(self):
        # reshape reco data once they've ben reconstructed
        pass

#todo revise this

    def tonifiti(self, path):
        """
        Write image data to nifti file
        Parameters
        ----------
        path - path to nifti file

        Returns
        -------
        -
        """

        nii = NiftiObject(bruker=self)
        nifti.writeNifti(path, nii)
        del nii

        return

    def parameters2dict(self):
        """

        Returns
        -------
        visu and reco parameters as a dictionary so that it can be easily diplayed
        """

        visu = {}
        reco = {}

        for key in dir(self.visu_pars):
            visu[key] = eval("self.visu_pars." + key)

        for key in dir(self.reco):
            reco[key] = eval("self.reco." + key)

        return visu, reco



class Study():
    def __init__(self, path=None, readData = False):

        if path[-1] != "/":
            self.path = path + '/'
        else:
            self.path = path


        # get a list of study folder's subfolders
        self.folders = walk(self.path).next()[1]

        for folder in self.folders:
            # create a scan from a folder
            obj = Scan(self.path + folder, readData=readData)

            # if the folder is a scan folder, set it as an attribute
            if obj.isScan:
                setattr(self, 'scan_' + folder, obj)

            # clean
            obj = None

    def export(self, mode='print', path=None):
        scans = dir(self)

        if mode == 'xls':
            # book = xlwt.Workbook(encoding="utf-8")
            print("XLS export is not supported at the moment")


        #todo keep, or not?
        '''
        for scan in scans:
            if scan.startswith('scan_'):

                if mode == 'xls':
                    exec ('sheet_' + scan + ' = book.add_sheet(\''+ scan +'\')')

                names, values = eval('self.' + scan + '.export(\'all\', less=True, mode=\'lists\')')

                names.insert(0,'File name')
                values.insert(0,scan.replace('scan_',''))

                if mode == 'xls':
                    row = 0

                    for name in names:
                        stringName = 'sheet_' + str(scan) + '.write('+ str(row) + ',0 , \'' + str(name) + '\')'
                        stringValue = 'sheet_' + str(scan) + '.write('+ str(row) + ' ,1, \'' + str(values[row]) + ' \')'
                        try:
                            exec (stringValue)
                            exec (stringName)
                            row += 1
                        except:
                            pass



        if mode == 'xls':
            book.save(path)
        '''

        return