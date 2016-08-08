#!/usr/bin/python
# -*- coding: utf-8 -*-

import cmath
from math import *
import numpy as np
from IBSmodel_GY import IBSmodel
import sys
import matplotlib.pyplot as plt

from scipy.constants import c

#~ def Fgeom(bl,ex,betastar,phi):
	#~ FF=1/(sqrt(1+(bl*phi/sqrt(ex*betastar)/2)**2));
	#~ return FF
    
fr = 11.2455e3

def Lumi_inst(Number_bunches, bunch_intensity_p, ex_norm_m, ey_norm_m, bl_4sigma_s, betastar_m, phi_full_rad, gamma):
    sigmaz_m = bl_4sigma_s/4.*c
    sigmax_m = np.sqrt(ex_norm_m/gamma*betastar_m)
    sigmay_m = np.sqrt(ey_norm_m/gamma*betastar_m)
    FF=1./(sqrt(1.+(sigmaz_m/sigmax_m*phi_full_rad/2.)**2));
    return FF*fr*Number_bunches*bunch_intensity_p**2/(sqrt(2*sigmax_m**2)*sqrt(2*sigmay_m**2))/2/pi



def LumiModel(gamma=None, betastar_m=None, phi_full_rad=None, bunch_intensityin_p=None, exin_norm_m=None, 
        eyin_norm_m=None, blin_4sigma_s=None, tFill_s=None,tauSRxy_s=None, tauSRl_s=None, sigmaBOff_m2=None, 
        VRF_V=None, IBSON=1, emitBU=None, BoffON=1, nIPs = 2., dt_s=15*60.):
     
    ex_norm_m = [exin_norm_m]
    ey_norm_m = [eyin_norm_m]
    bl_4sigma_s = [blin_4sigma_s]
    bunch_intensity_p = [bunch_intensityin_p]
    
    
    Lumi_init = Lumi_inst(Number_bunches=1., bunch_intensity_p=bunch_intensityin_p, 
                        ex_norm_m=exin_norm_m, ey_norm_m=eyin_norm_m, bl_4sigma_s=blin_4sigma_s, 
                        betastar_m=betastar_m, phi_full_rad=phi_full_rad, gamma=gamma)
    
    
    Luminosity_invm2s = [Lumi_init]
    tt_s = np.arange(0., tFill_s, dt_s)
    N_steps = len(tt_s)
    
    tauh = 0
    tauv = 0

    for i_step in xrange(1, N_steps):
        
        
        ex_IBS_norm_m, IBSx, bl_IBS_4sigma_s, IBSl, ey_IBS_norm_m = IBSmodel(IBSON, gamma, 
                                bunch_intensity_p[i_step-1],ex_norm_m[i_step-1],ey_norm_m[i_step-1],
                                bl_4sigma_s[i_step-1],VRF_V,dt_s)
        if emitBU=='Model':	
            ex_norm_m.append(ex_IBS_norm_m)
            ey_norm_m.append(ey_norm_m[i_step-1]*exp(-2*dt_s/tauSRxy_s));   
        elif emitBU=='EmpiricalBlowup':
                #raise ValueError('Not checked yet for units (Gianni&Yannis)')
                #ex0.append(ex_norm_m[i-1]+tauh*dt)
                #ey0.append(ey_norm_m[i-1]+tauv*dt);
            ex_norm_m.append(ex_norm_m[i_step-1]+tauh*dt_s)
            ey_norm_m.append(ey_norm_m[i_step-1]+tauv*dt_s);

        bl_4sigma_s.append(bl_IBS_4sigma_s)


        if BoffON==1:
            #~ tauBOff_a = nb*Nb0[i-1]/(nIPs*sigmaBOff*L0[i-1]));#   % in minutes
            tauBOff_s = bunch_intensity_p[i_step-1]/(sigmaBOff_m2*nIPs*Luminosity_invm2s[i_step-1])
            bunch_intensity_p.append(bunch_intensity_p[i_step-1]/(1.+dt_s/tauBOff_s));        
        else:
            bunch_intensity_p.append(bunch_intensity_p[i_step-1]);

        Luminosity_invm2s.append(Lumi_inst(Number_bunches=1., bunch_intensity_p=bunch_intensity_p[i_step], 
                        ex_norm_m=ex_norm_m[i_step], ey_norm_m=ey_norm_m[i_step], bl_4sigma_s=bl_4sigma_s[i_step], 
                        betastar_m=betastar_m, phi_full_rad=phi_full_rad, gamma=gamma))


    return tt_s, Luminosity_invm2s, bunch_intensity_p, ex_norm_m, ey_norm_m, bl_4sigma_s


