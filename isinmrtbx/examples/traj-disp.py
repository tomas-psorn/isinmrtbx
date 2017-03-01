from isinmrtbx.datatypes.BrukerDataTypes import Scan
from os import listdir

PATH = "/home/psorn/Data/27jan17_signalnejsi_guma_ute_1_1/"

paths = listdir(PATH)

for path in sorted(paths):
	try:
		scan = Scan(PATH + path)
		scan.showTrajectory(projStep=5, sampleStep=10)
	except:
		pass