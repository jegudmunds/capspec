import numpy as np
import pdata as pd
from matplotlib import pyplot as plt

# Remember that an often used unit for visosity is the 'Poise' = g/cm/s
# However, the SI unit for viscosity is (Pa s) but one (Pa s) = 10 Poise

def eta_he(T):
    # See Amdur1947
    # Units of Poise = g/cm/s = 0.1 Pa s
    eta0 = np.array([0.549, 1.277, 2.846, 3.147, 3.523, 8.205, 8.845, 9.470, 10.095,
                    10.680, 11.255, 11.815, 12.365, 12.90, 13.43, 13.95, 14.46, 14.96])
    T0 = np.array([1.64, 4.23, 14.21, 16.55, 20.38, 80, 90, 100, 110, 120, 130, 140,
                   150, 160, 170, 180, 190, 200])
    eta0 = np.append(eta0,19)
    eta0 = eta0*1e-6
    T0 = np.append(T0,300)
    eta_he = pd.extrap(T,T0,eta0)
    return eta_he

def eta_n2(T):
    # See Cole1985
    # Units of g/cm/s = 0.1 Pa s
    eta0 = np.array([8.24, 8.85, 9.45, 10.05, 12.88, 15.48, 17.90])
    T0 = np.array([120, 130, 140, 150, 200, 250, 300])
    eta0 = eta0*1e-6
    eta_n2 = pd.extrap(T,T0,eta0)
    return eta_n2

if __name__ == "__main__":
    T_He = np.linspace(1,300,100)
    T_N2 = np.linspace(120,300,100)
    plt.plot(T_He,eta_he(T_He))
    plt.plot(T_N2,eta_n2(T_N2))
    plt.show()
