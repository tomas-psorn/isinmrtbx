import sys
from isinmrtbx.datatypes.BrukerDataTypes import Study

PATH = sys.argv[1]

# PATH = "/home/psorn/Data/20170307_155607_Strukturni_bily_fantom_mrizka_1_2/"

Study(path=PATH).export(mode='xls', path= PATH + "/study-export.xls" )
