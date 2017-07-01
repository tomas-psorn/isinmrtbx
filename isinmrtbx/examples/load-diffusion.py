import sys

from isinmrtbx.datatypes.BrukerDataTypes import Scan, Reco

# create a scan object
scan = Scan("/media/tomas/OS/work_local/upt/data/rot/pl_16Apr20_ROT35_DKI.Bc1/5")

# get parameter from method
DwBMat = scan.method.PVM_DwBMat
print(DwBMat)

# get parameter from acqp
NI = scan.acqp.NI
print(NI)

# get kspace
kspace = scan.fid

# create a reco object
reco = Reco("/media/tomas/OS/work_local/upt/data/rot/pl_16Apr20_ROT35_DKI.Bc1/5/pdata/1")

# get parameter from reco
RECO_size = reco.reco.RECO_size
print(RECO_size)

# get parameter from visu_pars
VisuCoreFrameCount = reco.visu_pars.VisuCoreFrameCount
print(VisuCoreFrameCount)

# get image
image = reco.data2dseq

