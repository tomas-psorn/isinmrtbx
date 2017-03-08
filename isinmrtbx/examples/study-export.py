import sys
from isinmrtbx.datatypes.BrukerDataTypes import Study

PATH = sys.argv[1]

Study(path=PATH).export(mode='xls', path= PATH + "/study-export.xls" )
