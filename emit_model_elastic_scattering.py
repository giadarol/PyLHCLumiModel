import cmath
from math import *
import numpy as np
import sys
import matplotlib.pyplot as plt



#gamma = 7460.54

nIP = 2
Np = 1.1e11
nb = 2232
L = 5e33*1e-4
beta_star_m = 80/100 

sigma_el = 2.97e-26

def dedtElastic(nIP, Np, nb, L, beta_star_m, sigma_el_m2):

  L=L*1e-4
  beta_star = beta_star_m*100
  sigma_el = sigma_el_m2*1e4
  En = 13	# energy [TeV]
  B = 19.9
  P = 6500.

  t = 1./B

  gamma = 6927.64
  theta = np.sqrt(t/P**2)
  theta_x = theta/sqrt(2.)

  den_x_dt = 0.5*nIP*gamma*beta_star*theta_x**2*L*sigma_el/nb/Np
  den_x_dt_mmmradpersec = den_x_dt*10*1000  

  den_dt = 0.5*nIP*gamma*beta_star*theta**2*L*sigma_el/nb/Np
  den_dt_mmmradpersec = den_dt*10*1000 

  return den_x_dt_mmmradpersec*1e-6, den_dt_mmmradpersec*1e-6
