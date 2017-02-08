"""
Script demonstrating the difference between reconstruction performed using measured and computed trajectory
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

from isinmrtbx.datatypes.BrukerDataTypes import Scan, Reco

path_example = "data/ute"

scan = Scan(path_example)

projStep = 30
sampleStep = 2
x = 0
y = 1

plt.figure()

plt.axis('equal')
for projection in np.arange(0, scan.traj.shape[1], projStep):
	plt.stem(scan.traj[0::sampleStep, projection, x], scan.traj[0::sampleStep, projection, y], "+")

plt.title("Trajektorie sekvence "+ str(scan.acqp.PULPROG) + ", zobrazena každá " + str(projStep) + "-tá projekce a každý " + str(sampleStep) + "-tý vzorek")

# sio.savemat('/home/psorn/data/LEGOPHANTOM_BRUKER_1_3/29/traj.mat', {'traj':dataset.traj})

reco = Reco(path_example + "/pdata/1/")

plt.figure()

plt.title("reco")
plt.imshow(reco.data2dseq[:,:,0,0,0])

plt.show()



