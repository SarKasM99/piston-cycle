# Purpose 
This python script was written as part of a university assignment for the LMECA1210 course. The aim was to study various physical properties of a 4 stroke engine cycle, as well as estimate the theoretical thickness of a piston rod to avoid flexion. 

### Assumptions:
* The air in the engine is a perfect gas
* The combustion begins at a fix rod angle $\theta_c$
* The cycle begins with an air admission (from $\theta = -2\pi$ to $\theta = -\pi$)
> Followed by a compression (from $\theta = -\pi$ to $\theta = 0$) and a decompression (from $\theta = 0$ to $\theta = \pi$)
* The cycle ends with an exhaust phase (from $\theta = \pi$ to $\theta = 2\pi$)
* The ratio between the length of the connecting rod and the crank throw ($\beta$) is really large
* The gas is diatomic with a heat capacity ratio $\gamma=1.3$
* The atmospheric pressure is equal to 1 bar
# Requirements
The script uses the following python packages:
* numpy
* scipy
* matplotlib
# Usage
First, specify the following engine parameters at the beginning of the script:
```python
# An example of parameters
tau = 9.3        # Pressure ratio [-]
D = 0.082        # Diameter of the piston [m] 
C = 0.08         # Piston stroke [m] 
L = 0.136        # Length of the connecting rod [m] 
mpiston = 0.347  # Mass of the piston [kg]
mbielle = 0.6113 # Mass of the connecting rod [kg]
Q = 2800000      # Heat per mass of gas [J/kg_inlet gas]   
```
After specifying the parameters, use `myfunc` function to calculate the physical properties. You can find the following example at the bottom of the script:
```py
theta = np.linspace(-180,180,1000)
V_out, Q_out, Fp, Ft, p_out, t, dQ, dV = myfunc(2000,1,theta,24,40)
```
