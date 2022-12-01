#!/usr/bin/env python
# coding: utf-8

# In[9]:


import numpy as np
import matplotlib.pyplot as plt
import math
L = 2000000
time_interval = 0.005
trajs_number = 37


# In[10]:


import random
class Runge_kutta():
    """
    Class handles Runge_Kutta two dimensions trajectory
    """
    def __init__(self, L, x0, kB=1.38e-23, T=300, tau_m=10, tau_D=1000, tau_L=1, h=0.001 ):
        """
        The initial state for x, and R
        L: int
            Length of the trajectory
        x0: float
            initial position.
        R0: float
            initial R
        h: float
            the step size of Runge Kutta.
        U0: float
            set 3*kB*T
        """

        self.L = L
        self.x0 = x0
        self.U0 = 3*kB*T
        self.tau_m = tau_m
        self.tau_D = tau_D
        self.kB = kB
        self.tau_L = tau_L
        self.h = h
        self.T = T
        #self.fig = plt.figure()
        return
    
    def random_term(self):
        x = random.gauss(0.0, 2*self.tau_D)
        return x
    def f0(self, y):
        #x_dot = y
        return y
    def f1(self, x, R):
        #y_dot = f(x, R)
        return 1/self.tau_m*(4*self.U0/(self.kB*self.T)*x*(1-x**2)/self.tau_D + R)
    def f2(self, R, y):
        #R_dot = f(x_dot, R)
        return -1/self.tau_L*(R+y-self.random_term())

    def summation_RG(self, x, k1, k2, k3, k4):
        return x + k1/6 + k2/3 + k3/3 + k4/6 

    def stepfunction(self, x, y, R):
        # K1 = hf(x,y)
        # K2 = hf(x+h/2, y+k1/2)
        # k3 = hf(x+h/2, y+k2/2)
        # k4 = hf(x+h/2, y+k3)
        # y(n+1) = y(n)+ k1/6 + k2/3 +k3/3 +k4/6
        #first order
        K_01 = self.h*self.f0(y)
        #x_dot = y
        K_11 = self.h*self.f1(x, R)
        #y_dot = f(x,R)
        K_21 = self.h*self.f2(y, R)
        #R_dot = g(y,R)
        #second order
        K_02 = self.h*self.f0(y+K_11/2)
        K_12 = self.h*self.f1(x+K_01/2, R+K_21/2)
        K_22 = self.h*self.f2(y+K_11/2, R+K_21/2)
        #third order
        K_03 = self.h*self.f0(y+K_12/2)
        K_13 = self.h*self.f1(x+K_02/2, R+K_22/2)
        K_23 = self.h*self.f2(y+K_12/2, R+K_22/2)
        #fourth order
        K_04 = self.h*self.f0(y+K_13)
        K_14 = self.h*self.f1(x+K_03, R+K_23)
        K_24 = self.h*self.f2(y+K_13, R+K_23)
        
        x_new = self.summation_RG(x, K_01, K_02, K_03, K_04)
        y_new = self.summation_RG(y, K_11, K_12, K_13, K_14)
        R_new = self.summation_RG(R, K_21, K_22, K_23, K_24)

        return x_new, y_new, R_new
        

 
    def run(self):
        """
        h is the step size
        """
        x_list = []
        #y_list = []
        #R_list = []
        t_list = [i*self.h/self.tau_D for i in range(0, self.L+1)]
        
        R = random.gauss(0.0, 1.0/(self.tau_D*self.tau_L))
        x = self.x0
        y = random.gauss(0.0, 1.0/(self.tau_D*self.tau_m))
       # print("Initial value: x:", x)
        x_list.append(x)
        #y_list.append(y)
        #R_list.append(R)

        for i in range(self.L):
            x, y, R = self.stepfunction(x, y, R)
            x_list.append(x)
            #y_list.append(y)
            #R_list.append(R)
        # calculate mass by tau_m and tau_D
        mass = self.kB*self.T*self.tau_m*self.tau_D/self.L**2
        
        assert len(t_list) == len(x_list), "Length not cooperate"
        return t_list, x_list, mass


# In[11]:


    
#    plt.title(r'$ \tau_\Gamma/ \tau_D= $ %d' %c)
#    plt.plot(t, x)
#    plt.xlabel(r'$\frac{t}{\tau_D}$')
#    plt.show()
    
#np.savetxt("traj_2e7_memory.txt", [position[:-4], mass[:-4]], delimiter =", ")


# # Memory Kernel
# Our potential is here:
# $$U(x) = U_0[(\frac{x}{L})^2-1]^2.$$
# And derivative is:
# $$\triangledown U(x) =  2U_0[\frac{x^5}{L^3}-\frac{2x}{L}]$$

# In[12]:


class Memory_Kernel():
    def __init__(self, traj, time_interval, mass):
        """
        The potential function we must already know!
        traj: list
            trajectory
        time_interval: float 
            time interval of the trajectory.
        mass: float.
            particle mass.
        """
        self.traj = traj
        self.time_interval = time_interval
        self.mass = mass
        return
        
        
        
    def velocity(self, x, t):
        v = []
        for i in range(len(x)-2):
            v.append((x[i+2]-x[i])/(2*t))
        return v

    def acceleration(self, x, t):
        a = []
        for i in range(len(x)-2):
            a.append((x[i]+x[i+2]-2*x[i+1])/t**2)
        return a

    def Gradpotential(self, x):
        U0 = 1.0
        L = 1.0
        Gp = [2*U0*(2*pow(i, 3)/pow(L,3)-2*i/L) for i in x]
        return Gp

    def correlation(self, x, y):
        x = np.fft.fft(x)
        y = np.fft.fft(y)
        cor = []
        for i in range(len(x)):
            cor.append(x[i]*np.conj(y[i])/len(x))
        cor = np.fft.ifft(cor)
        return cor.real


    def run(self):
        v = self.velocity(self.traj, self.time_interval)
        a = self.acceleration(self.traj, self.time_interval)
        G_u = self.Gradpotential(self.traj)[1:-1]
        corvv = self.correlation(v, v)
        corva = self.correlation(v, a)
        corvgu = self.correlation(v, G_u)
        coraa = self.correlation(a,a)
        aa = coraa[0]
        coragu = self.correlation(a, G_u)
        b = coragu[0]
        Gamma_0 = (self.mass*aa-b)/corvv[0]
        Gamma_1 = -2/(self.time_interval*corvv[0])*(self.mass*corva[1]+self.time_interval/2*Gamma_0*corvv[1]+corvgu[1])
        temp_gamma = 0.0
        memory_kernel = [Gamma_0, Gamma_1]
        temp = 0.0
        pro_parameter = self.time_interval*corvv[0]
        length = int(10/time_interval)
        for i in range(2, length):
            for j in range(1, i):
                temp = temp + self.time_interval*memory_kernel[j]*corvv[i-j] 
            temp_gamma = -2/pro_parameter*(temp + self.mass*corva[i] + self.time_interval/2*Gamma_0*corvv[i]+ corvgu[i])
            memory_kernel.append(temp_gamma)
            temp_gamma = 0.0
            temp = 0.0
        
        
        print(aa)
        print(b)
        return memory_kernel


# In[13]:

    


# In[8]:


def FMPT(x, L, time_interval):
    """
    For double well:
    """
    upertime=0.0
    downtime=0.0
    jumpup_time = 0
    jumpdown_time = 0
    state = 0
    #state means the particle is in which well.
    for i in x:
        if state == 0 and i < L:
            downtime = downtime+time_interval
        if state == 0 and i > L:
            state = 1
            jumpup_time += 1
        if state == 1 and i > -1*L:
            upertime = upertime + time_interval
        if state == 1 and i < -1*L:
            state = 0
            jumpdown_time +=1
    print("Down time is", downtime)
    mean_tFirst_passage_time = downtime/(jumpup_time + 1)
    return mean_tFirst_passage_time


# In[5]:


import math
time = []
position = []
mass = []
tm_td = []
MFPT_list = []
for i in np.linspace(-2.7, 3, trajs_number):
    c = math.pow(10, i)/10
    
    test = Runge_kutta(L, x0=-1, tau_m=0.1, tau_D=10.0, tau_L=math.pow(10, i), h=time_interval)
    t, x, m= test.run()
    test = Memory_Kernel(x, time_interval, m)
    memory_kernel = test.run()
    time_list = [i*time_interval/10 for i in range(len(memory_kernel))]
    mid = []
    for j in range(1, round(len(memory_kernel)/2)):
        mid.append((memory_kernel[2*j]+memory_kernel[2*j+1])/2)
    plt.plot(mid)
#mem.append(memory_kernel)
    plt.savefig(f'D_10_m_100_Gamma_m2.5_to_p3_%f.png' %c)
    plt.cla()
    MFPT_list.append(FMPT(x, 1, time_interval))
    tm_td.append(c)

plt.scatter(tm_td,MFPT_list)
plt.xscale('log')
plt.yscale('log')
plt.savefig('D_10_m_10_Gamma_m2.5_to_p3_MFPT.png')

# In[ ]:




