from isinmrtbx.datatypes.BrukerDataTypes import Scan
from os import listdir

PATH = "FILE with multiple UTE scans"

paths = listdir(PATH)

for path in sorted(paths):
	try:
		scan = Scan(PATH + path)
		scan.showTrajectory(projStep=5, sampleStep=10)
	except:
		pass