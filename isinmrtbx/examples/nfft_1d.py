from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from isinmrtbx.tools.recos import gridding1_1


#defaults
alfa = 2.
W = 2.5
data_N = 40

# sampling function

traj = np.linspace(0, 0.5, data_N/2) ** 2
traj = np.append(traj[::-1] * -1, traj)
signal = np.sin(np.pi* 4. * traj)

interp = gridding1_1(signal, traj, alfa, W)

plt.figure()
plt.stem(traj, signal)
plt.show()



