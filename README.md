# Markovian Embedding model in Double-well potential.
Two main parameters: $\tau_m=\frac{m}{\tau_D}$, $\tau_\Gamma$. 
$$\Gamma=\frac{\gamma}{\tau_\Gamma}e^{-\frac{t}{\tau_\Gamma}}$$
$\tau_D=\frac{L^2\gamma}{k_BT}$ is the friction-related diffusion time, the time taken to diffuse a distance $L$. $k_BT$ is thermal energy. $\gamma$ is the friction of system.

By varying these two parameters, you can change the inertia and memory effect in the system.


## Trajectory making
This will be in the document named (timestep_changing)
We put a particle in a double well potenital(potential.cpp), 
![alt text](https://github.com/Sanhei/A-naive-memory-extraction/blob/main/potential.png?raw=true)

using the Markovian embedding system, which the discretizing equation is
$$x(t_{j+1})=2x(t_j)-\frac{\tau_D}{\tau_m}F(x(t))t^2-\frac{\tau_D^2}{\tau_m\gamma}(y(t_j)-x(t_j))\Delta t^2.$$
$$y(t_{j+1})=y(t_j)-\tau_D/\gamma(y(t_j)-x(t_j))\Delta t+\sqrt{2\Delta t}\xi$$

In the document parameter.h, it saved all the parameter we need; the integrator function is stored in the v_verlet.cpp (see the aboved function). And in the traj.cpp we just do the iteration.

And it will generate a .txt file which record each step of time position. And a single example:
![alt text](https://github.com/Sanhei/A-naive-memory-extraction/blob/main/0.005_traj.svg?raw=true)
## Trajectory analyse.
We need to extract the memory kernel, which is 
$$\gamma = \int^\infty_0\Gamma(t)$$.
To calculate this part, we need to get the correlation function of $C_{vv}$ and $C_{\triangledown U x}$.
### Potential gradient $\triangledown U$
For calculating $\triangledown U$, we use free energy by calculating the distribution of the position. Therefore, potential.h gives three function.
1. Histogram calculation, which will record the distribution of position in the trajectory. And $$F = -k_BTln(P(x))$$, so this can directly transfer to potential.
2. Gradient of potential.
3. Calculate a small area around boundary, because the behavior is not like stochastic behavior. So will print out how much influence on the whole calculation. And add a correction, which may a little bit shift due to numerical error.



### Correlation calculation
This part, we do correlation function. First get the velocity from the time derivative of trajectory. Then we get everything to calculate the correlation function. For calculating the correlation, we will do three steps for that:
1. Do the Fast Fourier Tranfer(FFT) to the variables. Here in our case, they are velocity, position and derivative of potential.
2. Get the convolution product of the variables(which got from FFT).
3. Reverse FFT to the convolution product, then we will get correlation function.

## Memory kernel Calculation.
Here we use the equation from paper(Non-markov protien folding).

$$G_n = \frac{2}{\Delta t C^{vv}_0} (C^{\triangledown U q}_n -\frac{C^{\triangledown U q}}{C^{vv}_0} - \Delta t\sum^{n-1}_{i=1}G_{n-i}C^{vv}_i)$$. 
![alt text](https://github.com/Sanhei/A-naive-memory-extraction/blob/main/Gplot.png?raw=true)

(Somehow it cannot show the right format of the equation.)Then all the calculation done.

## Program explain.
In our program, we can call matplotlib function in python to directly plot the Potential and Correlaiton. And in the end, after the calculation, program can automatically call the fitting.py, which uses lmfit package to do the fitting. Here we may need the python package path to include this environment.
## Runing the program.
Environment:
GCC 14++ or higher version.
FFTW(boost/math for spline, akima spline need to set the boundary, but you can delete it). 
For fitting part and plotting, we use python. There we may use the package: lmfit, matplotlib, numpy

```bash
export PYTHONPATH="~/python_path/lib/python3.x/site-packages"
```
This is to include the path of python, but this step may not necessary, since fourier transfer can be very time consuming, therefore, we want to see if the parameter setting is proper. So we may want to see if the free energy profile is same as we expected(a double well potential showed above), before the Fourier transfer.
By default, we comment all the code of python. The simulation is all done by C++ code. Setting the environment can be frustrated. Another thing you need is FFTW.

```bash
cd timestep_changing
g++ -o traj traj.cpp
./traj
```
This part is for the simulation, create a .txt(named by your discretised time step) file as your trajectory.
```bash
cd ..
make
./memoryextraction
```
This file will generate a kernel.txt which is the memory kernel. Then the correlation data and all figures will in the "figure" director. If you want to change the figure results, you can do it by changing fitting.py. 

## Fitting Result
The red curve is exponential fitting curve and blue solid line is input memory kernel. Here, we proved our model is self-consistant.
![alt text](https://github.com/Sanhei/GLE-Model-validation-and-analysis/blob/main/Memorykernel.svg?raw=true)
## Space and Memory
For now, we use 2GB trajectory, for several parameter settings, it works. However it is not enough, so maybe we want to try 8GB trajectory. As for memory occupied, FFT may use three times bigger than the trajectory size. So for each time correlation function, we record it in a text file, and read it again.
The thing we need to care about is $C_{\triangledown U x}$, We need to record an additional trajectory size. So in total, it will $\mathbf{8}$ times bigger memory. For the optimization, we just use one core, and didn't set the threads. For FFTW, it could use multi threads.
