from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from isinmrtbx.tools.recos import gridding1_1


#defaults
alfa = 5
W = 5
data_N = 1000

# sampling function

traj = np.linspace(0, 0.5, data_N/2) ** 2
traj = np.append(traj[::-1] * -1, traj) * 2 #merge halves and stretch it back to 0.5
signal = np.sin(np.pi* 4. * traj)

interp = gridding1_1(signal, traj, alfa, W)

grid = np.linspace(-.5, .5, interp.shape[0])

plt.figure()
plt.plot(traj, signal)
plt.plot(grid, interp, 'r')
plt.show()



