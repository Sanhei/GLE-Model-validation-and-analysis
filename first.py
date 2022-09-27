import numpy as np
import matplotlib.pyplot as plt
import matplotlib

with open('timestep_changing/parameter.h') as f:
    for line in f:
        print(line)
        splitByComma=line.split('=')
        if len(splitByComma)>1:
            
            if(splitByComma[0]=="const double tau   "):
                tau = float(splitByComma[1].strip().replace(';', ''))
                print("tau = ", tau)
                break
            if(splitByComma[0]=="const double tau_D "):
                tau_D = float(splitByComma[1].strip().replace(';', ''))
                print("tau_D = ", tau_D)    
# To get the simulation time step
with open('timestep_changing/traj.cpp') as f:
    for line in f:
        #print(line)
        splitByComma=line.split('=')            
        if splitByComma[0]=="    double h":
            time = splitByComma[1].strip()
            print("time interval = ", time)
h=0
h = float(time[-2])*10**(-len(time)+3)

hist_MFPT = np.loadtxt("./figure/hist_MFPT.txt")



plt.scatter(hist_MFPT[:, 0], hist_MFPT[:, 1], label="First-first passage time")
plt.title("Passage time, $\\tau_D$="+ str(tau_D)+"$\\tau_\Gamma$="+str(tau))
plt.xscale("log")
plt.xlabel("Passage time ")
plt.ylabel("Distribution probability")
plt.title("MFFPT")
plt.legend()
plt.savefig("./figure/MFPT.png")
plt.clf()


hist_numb = np.loadtxt("./figure/hist_num_cross.txt")
hist_numb = hist_numb[:200000]
plt.title("Crossing barrier distribution, $\\tau_D$="+ str(tau_D)+"$\\tau_\Gamma$="+str(tau))
plt.scatter(hist_numb[:, 0], hist_numb[:, 1])
plt.xlabel("Crossing barrier times")
plt.ylabel("Probability")
plt.yscale("log")
plt.xscale("log")
plt.savefig("./figure/crossing_barrier.png")

plt.clf()


print("Fitting beginning")
print("tau_D = ", tau_D)
print("tau = ", tau)
from lmfit import Model
#Use package 
#Input function
def func(x, b):
    return b*np.exp(x*(-b))
kernel = np.loadtxt("kernel.txt").T
print("Value is", kernel[0][0])
expect = [func(i, tau_D/tau) for i in kernel[0]]
print("Test is", expect[0])
model = Model(func)

params = model.make_params(b=10)
#params['a'].min = -2
#params['b'].max = 0.00030
result = model.fit(kernel[1], params, x = kernel[0])
print(result.fit_report())
plt.plot(kernel[0], expect, label = 'expectation')
plt.plot(kernel[0], kernel[1], 'go', markersize=1, label = 'data point')
plt.plot(kernel[0], result.best_fit, 'r--', label = 'fit:')
plt.legend()
plt.xlabel("time($\\tau_D$)")
plt.xscale("log")
#plt.yscale("log")
plt.ylabel(r"Memory kernel")
plt.savefig("./figure/Memorykernel.svg")
plt.clf()



with open('timestep_changing/traj.cpp') as f:
    for line in f:
        #print(line)
        splitByComma=line.split('=')            
        if splitByComma[0]=="    double h":
            time = splitByComma[1].strip()
            print("time interval = ", time)
h=0
h = float(time[-2])*10**(-len(time)+3)
# reduce the size, every 30 steps record the trajectory.
record_interval = 30;
print("trajectory record interval is ", record_interval)

trajectory = np.loadtxt("./figure/small_traj.txt")
trajectory = trajectory[::record_interval]
print("Read trajectory finished")
time = []
for i in range(len(trajectory)):
    time.append(i*h*record_interval)
plt.plot(time, trajectory)
plt.title("trajectory, $\\tau_D$="+ str(tau_D)+"$\\tau_\Gamma$="+str(tau))
plt.xlabel("time($\\tau_D$)")
plt.ylabel("Position")
plt.savefig("./figure/trajectory.png")
plt.clf()
