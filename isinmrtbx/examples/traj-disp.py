from isinmrtbx.datatypes.BrukerDataTypes import Scan, Study
from os import listdir

PATH = "/home/psorn/Data/20170313_090733_UTE_relaxometrie_MnCl2_guma_1_2/19/"


# paths = listdir(PATH)
#
# for path in sorted(paths):
#
#  	scan = Scan(PATH + path)
#  	scan.showTrajectory(projStep=5, sampleStep=10)

scan = Scan(PATH, readFid=False)
scan.showTrajectory(projStep=2, sampleStep=5)



	# try:
	# 	scan = Scan(PATH + path)
	# 	scan.showTrajectory(projStep=5, sampleStep=10)
	# except:
	# 	pass

