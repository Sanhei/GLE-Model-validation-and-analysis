import matplotlib.pyplot as plt
import matplotlib
import numpy as np
matplotlib.use('TkAgg')
traj = np.loadtxt("test2.txt")
figure = plt.Figure()
plt.plot(traj)
plt.show()

