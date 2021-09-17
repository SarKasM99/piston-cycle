import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

#Parameters of the engine
tau = 9.3 #[-]                  #Compression ration
D = 0.082 #[m]                  #Diameter of the piston
C = 0.08 #[m]                   #Piston stroke
L = 0.136 #[m]                  #Length of the connecting rod
mpiston = 0.347 #[kg]           #Mass of the piston
mbielle = 0.6113 #[kg]          #Mass of the connecting rod
Q = 2800000 #[J/kg_inlet gas]   #Heat per mass of gas

#Constants
R = C/2     #[m]
m_mol = 0.02897  #[kg/mol]
Rgaz = 8.3145  #[J/K*mol]
Rmas = Rgaz/m_mol #[J/kg*K]
T = 303.15     #[K]

def Fcrit(Fp,Ft):
    """
    Fcrit is a function that calculates the maximum force of compression experienced by a connecting rod

    :param Fp:
    An array that represents forces experienced by the lower part of the rod (pied in french) during the engine cycle
    :param Ft:
    An array that represents forces experienced by the upper part of the rod (tete in french) during the engine cycle

    :return: F_crit:
    The maximum force of compression experienced by the rod
    """
    F_pair = np.stack((Fp, Ft), axis=1)
    F_pair = np.abs(F_pair[F_pair[:, 0] * F_pair[:,1] < 0])  # Keeps only opposite forces and converts all forces to their absolute value
    res = np.minimum(F_pair[:, 0], F_pair[:, 1])
    F_crit = np.max(res)
    return F_crit

def t_crit(Fc,plan):
    """
    t_crit is a function that calculates the theoretical thickness of the connecting rod

    :param Fc:
    The maximum force of compression
    :param plan:
    The direction of the buckling (either "x" or "y")

    :return:
    The theoretical thickness of the rod
    """
    E = 2e+11  # Module of elasticity [Pa]
    sigma_c = 4.5e+8  # Resistance to compression [Pa]
    if(plan == "x"):
        alpha = (419 * np.pi*np.pi * E) / (12 * L*L)
    elif(plan == "y"):
        alpha = (524 * np.pi * np.pi * E) / (12 * L*L)

    coef = np.array([11*alpha*sigma_c,(-alpha)*Fc,(-11)*Fc*sigma_c])

    roots = np.roots(coef)
    roots = roots[roots > 0]
    return np.max(np.sqrt(roots))

def myfunc(rpm, s, theta, thetaC, deltaThetaC):
    """
    myfunc calculates different physical properties of the engine cycle.

    :param rpm:
    Rotations per minute of the engine
    :param s:
    Ration between the initial pressure and atmospheric pressure
    :param theta:
    An array that represents the angle of rotation of the connecting rod in degrees
    :param thetaC:
    The angle at which the combustion begins in degrees
    :param deltaThetaC:
    The duration of the combustion in degrees

    :return:
    This function returns a tuple that contains the following elements respectively:
    -V_output is an array that represents the volume of the cylinder during the cycle (same size as theta)
    -Q_output is an array that represents the heat of inside the cylinder during the cycle (same size as theta)
    -F_pied_output is an array that represents the force experienced by the upper part of the connecting rod (same size as theta)
    -F_tete_output is an array that represents the force experienced by the lower part of the connecting rod (same size as theta)
    -p_output is an array that represents the pressure experienced by the piston during the cycle (same size as theta)
    -t is the theoretical minimal thickness of the connecting rod
    -dQ and dV represents the derivatives of Q_output and V_output
    """
    Vc = (np.pi * C * D * D) / 4
    beta = L / R
    p_in = s * 100000     #Initial pressure
    m = p_in*Vc /(Rmas*T)             #The mass of the gas

    thetaC = -(thetaC*(np.pi)/180)
    deltaThetaC = (deltaThetaC*(np.pi)/180)
    theta = (theta*(np.pi)/180)
    w = (2*np.pi*rpm)/60                     #Angular frequency [rad/s]

    V_output = (Vc / 2) * (1 - np.cos(theta) + beta - np.sqrt(beta * beta - np.sin(theta) * np.sin(theta))) + ( Vc / (tau - 1))

    Q_output = ( (Q / 2) * (1 - np.cos(np.pi * ((theta - thetaC) / deltaThetaC))) ) * m

    Q_output[theta < thetaC] = 0
    Q_output[theta > thetaC + deltaThetaC] = 0

    dQ = ((Q * np.pi) / (2 * deltaThetaC) * (np.sin(np.pi * ((theta - thetaC) / (deltaThetaC))))) * m
    dQ[theta < thetaC] = 0
    dQ[theta > thetaC + deltaThetaC] = 0

    dV = (Vc / 2) * np.sin(theta) * (1 + np.cos(theta) / np.sqrt(beta * beta - np.sin(theta) * np.sin(theta)))

    def dP(t, P):
        V_output = (Vc / 2) * (1 - np.cos(t) + beta - np.sqrt(beta * beta - np.sin(t) * np.sin(t))) + (Vc / (tau - 1))
        dV = (Vc / 2) * np.sin(t) * (1 + np.cos(t) / np.sqrt(beta * beta - np.sin(t) * np.sin(t)))
        dQ = 0
        if (t > thetaC and t < thetaC + deltaThetaC):
            dQ = ( (Q * np.pi) / (2 * deltaThetaC) * (np.sin(np.pi * ((t - thetaC) / (deltaThetaC) ))) ) * m
        return -1.3 * (P / V_output) * dV + (0.3) * dQ / V_output

    sol = integrate.solve_ivp(dP, (theta[0], theta[-1]), [p_in], t_eval=theta)
    p_output = sol.y[0]
    F_pied_output = (np.pi * D*D * p_output) / 4 - mpiston*R*w*w*np.cos(theta)
    F_tete_output = -(np.pi * D*D * p_output) / 4 + (mpiston + mbielle)*R*w*w*np.cos(theta)

    F_crit = Fcrit(F_pied_output,F_tete_output)
    t_x = t_crit(F_crit,"x")
    t_y = t_crit(F_crit,"y")
    t = np.maximum(t_x,t_y)   #We pick the maximum thickness
    return (V_output, Q_output, F_pied_output, F_tete_output, p_output, t,dQ,dV);

#Here is an example on how you can use our function
theta = np.linspace(-180,180,1000)
V_out,Q_out,Fp,Ft,p_out,t,dQ,dV = myfunc(2000,1,theta,24,40)

#plt.plot(theta,dV,'g-')

#plt.plot(theta,V_out,'g-')
#plt.plot(theta,Q_out)
#plt.plot(theta,dQ,'g')
#plt.plot(theta,p_out,'r')
#plt.plot(theta,p2,'g')

plt.plot(theta,Fp,'r')
plt.plot(theta,Ft,'g')
print(t)
plt.grid('axes=both')
plt.show()
