import matplotlib.pyplot as plt
import numpy as np

trajectory = np.loadtxt("0.000500.txt")
traj = []
for i in range((len(trajectory)/10):
	traj.append(trajectory[10*i])
x = [4*i*0.0005 for i in range(len(traj))]
plt.plot(x,trajectory[:100000000] )
plt.xlabel("time($\\tau_D$)")
plt.ylabel("position")
plt.title("time interval = 0.0005, tau=0.01")
plt.savefig("0.0005_traj.svg")
