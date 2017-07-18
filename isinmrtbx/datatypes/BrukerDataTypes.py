"""
Subject
Scan
Reco
"""
from __future__ import print_function
from __future__ import unicode_literals

import sys

import numpy as np

from os import listdir
from os import walk

import matplotlib.pyplot as plt

import xlwt

from isinmrtbx.datatypes.GeneralDataTypes import ParameterClass
from isinmrtbx.inout import readers
from isinmrtbx.tools import recos
from isinmrtbx.inout.codeXchange import codeXchange
from isinmrtbx.datatypes.NiftiDataTypes import NiftiObject, NiftiHdr

AUTHOR = 'Tomas Psorn'
CONTACT = 'tomaspsorn@isibrno.cz'


class Scan():

	def __init__(self, path, readFid=True, readReco=False, readTraj=True):

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
			self.fid = readers.readBrukerFidFile (self.path + 'fid', self)
			self.methodBasedReshape()

		if self.isTraj and readTraj:
			self.traj = readers.readBrukerTrajFile (self.path + 'traj',    self)

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

		if self.acqp.PULPROG in ['UTE.ppg','macicekGAUTE.ppg']:
			readers.fidHandle_UTE(self)
		elif self.acqp.PULPROG in ['FAIR_RARE.ppg']:
			readers.fidHandle_FAIR_RARE(self)
		elif self.acqp.PULPROG in ['DtiEpi.ppg', 'EPI.ppg', 'navigatorEPI_OM.ppg']:
			readers.fidHandle_Epi(self)
		elif self.acqp.PULPROG in ['FLASH.ppg']:
			readers.fidHandle_Flash(self)
		elif self.acqp.PULPROG in ['PRESS.ppg']:
			readers.fidHandle_PRESS(self)
		else:
			print('Function to reshape fid data of ' + self.acqp.PULPROG + ' sequence is not developed yet')

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
		Pvtools original function to determine number of channels used for acquisition

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

	def showTrajectory(self, projStep = 1, sampleStep = 1):
		try:
			traj = self.traj
			print("Showing trajectory of file: " + self.path)
			plt.figure()
			plt.axis('equal')
			for projection in np.arange(0, traj.shape[1], projStep):
				plt.stem(traj[0::sampleStep, projection, 0], traj[0::sampleStep, projection, 1], "+")
			plt.show()
		except:
			print("No traj info in: " + self.path)


	def export(self, attribute, mode='print'):
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
			if _attribute == 'acqp':
				for _name in dir(self.acqp):
					if _name[0] is not "_":
						namesGlobal.append(_name)
						_value = str(getattr(self.acqp, _name))
						valuesGlobal.append(_value)

			elif _attribute == 'method':
				for _name in dir(self.method):
					if _name[0] is not "_":
						namesGlobal.append(_name)
						_value = str(getattr(self.method, _name))
						valuesGlobal.append(_value)

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

			if path[-1] != "/":
				self.path = path + '/'
			else:
				self.path = path

			self.data2dseq = Reco()
			self.browseReco()

			if self.isVisu:
				self.visu_pars = ParameterClass()
				readers.readBrukerParamFile(self.path + 'visu_pars', self)

			if self.isReco:
				self.reco = ParameterClass()
				readers.readBrukerParamFile(self.path + 'reco', self)

			if self.is2dseq and self.isReco and self.isVisu:
				self.data2dseq = readers.readBruker2dseq(self.path + '2dseq', self)
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

	def getNiiHdr(self):

		hdr = NiftiHdr()

		# hdr.sizeof_hdr = DEFAULT

		# hdr.data_type = DEFAULT

		# hdr.self.db_name = DEFAULT

		# hdr.self.extents = DEFAULT

		# hdr.session_error = DEFAULT

		# hdr.regular = DEFAULT

		# hdr.dim_info = DEFAULT

		# IMAGE DIM

		# hdr.dim
		_data_dims = np.asarray(self.data2dseq.shape)
		_data_dims = _data_dims[_data_dims != 1] 	# only populated dimensions
		dim = np.ones(8, dtype='i2', order='F')
		dim[0] = len(_data_dims)
		dim[1:1 + len(_data_dims)] = _data_dims
		hdr.dim = dim

		# hdr.intent_p1 = DEFAULT

		# hdr.intent_p2 = DEFAULT

		# hdr.intent_p3 = DEFAULT

		#hdr.datatype
		hdr.datatype = codeXchange('int32', 'datatype')
		# hdr.datatype = codeXchange(self.data2dseq.dtype, 'datatype')

		# hdr.bitpix
		hdr.bitpix = self.data2dseq.dtype.itemsize * 8

		# hdr.slice_start
		hdr.slice_start = 0

		# hdr.pixdim
		fov = self.visu_pars.VisuCoreExtent
		matrixSize = self.visu_pars.VisuCoreSize
		frameThickness = self.visu_pars.VisuCoreFrameThickness
		_pixdim = np.ones(8, dtype='f4', order='F')
		_pixdim[0] = 0.0
		_pixdim[1] = fov[0] / float(matrixSize[0])
		_pixdim[2] = fov[1] / float(matrixSize[1])
		_pixdim[3] = frameThickness
		hdr.pixdim = _pixdim

		# hdr.vox_offset = DEFAULT

		# hdr.scl_slope
		hdr.scl_slope = self.visu_pars.VisuCoreDataSlope[0]

		# hdr.scl_inter
		hdr.scl_inter = self.visu_pars.VisuCoreDataOffs[0]

		# hdr.slice_end
		hdr.slice_end = hdr.dim[3]

		# hdr.slice_code
		hdr.slice_code = '1'

		# hdr.xyzt_units mm
		hdr.xyzt_units = '2'

		# hdr.cal_max
		hdr.cal_max = np.amax(hdr.scl_slope * self.data2dseq + hdr.scl_inter)

		# hdr.cal_min
		hdr.cal_min = np.amin(hdr.scl_slope * self.data2dseq + hdr.scl_inter)

		# hdr.slice_duration = DEFAULT

		# hdr.toffset = DEFAULT

		# hdr.glmax = DEFAULT

		# hdr.glmin = DEFAULT

		# FILE HISTORY

		# hdr.descrip = DEFAULT

		# hdr.aux_file = DEFAULT

		# hdr.sform_code - scanner_anat
		hdr.sform_code = 1

		# hdr.qform_code
		hdr.qform_code = 1

		# hdr.quatern
		rot_x = 0.
		rot_y = 0.
		rot_z = 0.
		R_X = np.array([[1., 0., 0.], [0., np.cos(rot_x), np.sin(rot_x)], [0., np.sin(rot_x), np.cos(rot_x)]])
		R_Y = np.array([[np.cos(rot_y), np.sin(rot_y), 0.], [0., 1., 0.], [np.sin(rot_y), np.cos(rot_y), 0.]])
		R_Z = np.array([[np.cos(rot_z), - np.sin(rot_z), 0.], [np.sin(rot_z), np.cos(rot_z), 0.], [0., 0., 1.]])
		R = R_X * R_Y * R_Z

		a = 0.5 * np.sqrt(1. + R[0,0] + R[1,1] + R[2,2])
		hdr.quatern_b = 0.25 * (R[2,1] - R[1,2]) / a
		hdr.quatern_c = 0.25 * (R[0, 2] - R[2, 0]) / a
		hdr.quatern_d = 0.25 * (R[1, 0] - R[0, 1]) / a

		# hdr.qoffset
		hdr.qoffset_x = 0.
		hdr.qoffset_y = 0.
		hdr.qoffset_z = 0.

		hdr.srow_x = np.zeros(4, dtype='float32', order='F')
		hdr.srow_y = np.zeros(4, dtype='float32', order='F')
		hdr.srow_z = np.zeros(4, dtype='float32', order='F')

		return hdr



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
		hdr = self.getNiiHdr()
		nii = NiftiObject(hdr=hdr, data=self.data2dseq)
		nii.write(path)
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

		self.folders = next(walk(self.path))[1]

		for folder in self.folders:
			# create a scan from a folder
			obj = Scan(self.path + folder, readFid=readData, readTraj=readData)

			# if the folder is a scan folder, set it as an attribute
			if obj.isScan:
				setattr(self, 'scan_' + folder, obj)

			# clean
			obj = None

	def export(self, mode='print', path=None):

		scans = dir(self)

		if mode == 'xls':
			book = xlwt.Workbook(encoding="utf-8")

		for scan in scans:
			if scan.startswith('scan_'):

				if mode == 'xls':
					exec ('sheet_' + scan + ' = book.add_sheet(\''+ scan +'\')')

				names, values = eval('self.' + scan + '.export(\'all\', mode=\'lists\')')

				names.insert(0,'File name')
				values.insert(0,scan.replace('scan_',''))

				if mode == 'xls':
					row = 0

					for name in names:
						stringName = "sheet_" + str(scan) + ".write("+ str(row) + ",0 , \"" + str(name) + "\")"
						stringValue = "sheet_" + str(scan) + ".write("+ str(row) + " ,1 , \"" + str(values[row]) + " \")"

						try:
							exec (stringValue)
							exec(stringName)
							row += 1
						except:
							exec(stringName)
							stringValue = "sheet_" + str(scan) + ".write("+ str(row) + " ,1 , \" Value too long \")"
							row += 1



		if mode == 'xls':
			book.save(path)

		return
