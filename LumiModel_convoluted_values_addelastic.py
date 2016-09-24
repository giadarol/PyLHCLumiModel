#!/usr/bin/python


import numpy as np
from IBSmodel_GY import IBSmodel
import sys
import matplotlib.pyplot as plt
from emit_model_elastic_scattering import dedtElastic

from scipy.constants import c

#~ def Fgeom(bl,ex,betastar,phi):
	#~ FF=1/(sqrt(1+(bl*phi/sqrt(ex*betastar)/2)**2));
	#~ return FF
    
fr = 11.2455e3

def Lumi_inst(Number_bunches, bunch_intensity_p1, bunch_intensity_p2, ex_norm_m1, ey_norm_m1, ex_norm_m2, ey_norm_m2, bl_4sigma_s1, bl_4sigma_s2, betastar_m, phi_full_rad, gamma, experiment):

    sigmaz_m = (bl_4sigma_s1+bl_4sigma_s2)/2./4.*c
    sigmax_m = np.sqrt(ex_norm_m1/gamma*betastar_m+ex_norm_m2/gamma*betastar_m)/np.sqrt(2.)
    sigmay_m = np.sqrt(ey_norm_m1/gamma*betastar_m+ey_norm_m2/gamma*betastar_m)/np.sqrt(2.)

    if not experiment: 
        experiment = 'CMS'

    if experiment=='ATLAS':
        FF=1./((1.+(sigmaz_m/sigmay_m*phi_full_rad/2.)**2)**0.5);
    elif experiment=='CMS':
        FF=1./((1.+(sigmaz_m/sigmax_m*phi_full_rad/2.)**2)**.5);

    return FF*fr*Number_bunches*bunch_intensity_p1*bunch_intensity_p2/(sigmax_m*sigmay_m)/4./np.pi


def LumiModel(gamma=None, betastar_m=None, phi_full_rad_ATLAS=None, phi_full_rad_CMS=None, experiment=None, bunch_intensityin_p1=None, bunch_intensityin_p2=None, exin_norm_m1=None,	exin_norm_m2=None, eyin_norm_m1=None, eyin_norm_m2=None, blin_4sigma_s1=None, blin_4sigma_s2=None, tFill_s=None,tauSRxy_s=None, tauSRl_s=None, sigmaBOff_m2=None, sigmaElastic_m2=None, VRF_V=None, IBSON=1, emitBU="Model", BoffON=1, nIPs = 2., dt_s=15*60., tau_empirical_h1=None, tau_empirical_v1=None, tau_empirical_h2=None, tau_empirical_v2=None):
     
    ex_norm_m1 = [exin_norm_m1]
    ey_norm_m1 = [eyin_norm_m1]
    ex_norm_m2 = [exin_norm_m2]
    ey_norm_m2 = [eyin_norm_m2]
    bl_4sigma_s1 = [blin_4sigma_s1]
    bl_4sigma_s2 = [blin_4sigma_s2]
    bunch_intensity_p1 = [bunch_intensityin_p1]
    bunch_intensity_p2 = [bunch_intensityin_p2]
    
    
    Lumi_init_ATLAS = Lumi_inst(Number_bunches=1., bunch_intensity_p1=bunch_intensityin_p1, bunch_intensity_p2=bunch_intensityin_p2, 
                        ex_norm_m1=exin_norm_m1, ey_norm_m1=eyin_norm_m1, ex_norm_m2=exin_norm_m2, ey_norm_m2=eyin_norm_m2, 
                        bl_4sigma_s1=blin_4sigma_s1, bl_4sigma_s2=blin_4sigma_s2, 
                        betastar_m=betastar_m, phi_full_rad=phi_full_rad_ATLAS, gamma=gamma, experiment='ATLAS')

    Lumi_init_CMS = Lumi_inst(Number_bunches=1., bunch_intensity_p1=bunch_intensityin_p1, bunch_intensity_p2=bunch_intensityin_p2, 
                        ex_norm_m1=exin_norm_m1, ey_norm_m1=eyin_norm_m1, ex_norm_m2=exin_norm_m2, ey_norm_m2=eyin_norm_m2, 
                        bl_4sigma_s1=blin_4sigma_s1, bl_4sigma_s2=blin_4sigma_s2, 
                        betastar_m=betastar_m, phi_full_rad=phi_full_rad_CMS, gamma=gamma, experiment='CMS')

    
    Luminosity_ATLAS_invm2s = [Lumi_init_ATLAS]
    Luminosity_CMS_invm2s = [Lumi_init_CMS]    
    Luminosity_invm2s = [(Lumi_init_ATLAS + Lumi_init_CMS)/2.]
    tt_s = np.arange(0., tFill_s, dt_s)
    N_steps = len(tt_s)
    
    for i_step in xrange(1, N_steps):
        #print i_step, N_steps        
       
        ex_IBS_norm_m1, IBSx1, bl_IBS_4sigma_s1, IBSl1, ey_IBS_norm_m1 = IBSmodel(IBSON, gamma, 
                                bunch_intensity_p1[i_step-1],ex_norm_m1[i_step-1],ey_norm_m1[i_step-1],
                                bl_4sigma_s1[i_step-1],VRF_V,dt_s)
        ex_IBS_norm_m2, IBSx2, bl_IBS_4sigma_s2, IBSl2, ey_IBS_norm_m2 = IBSmodel(IBSON, gamma, 
                                bunch_intensity_p2[i_step-1],ex_norm_m2[i_step-1],ey_norm_m2[i_step-1],
                                bl_4sigma_s2[i_step-1],VRF_V,dt_s)

        if emitBU=='Model':
	    dexdt1, deconvdt1 = dedtElastic(2, np.array(bunch_intensity_p1[i_step-1]), 1, np.array(Luminosity_invm2s[i_step-1]), betastar_m, sigmaElastic_m2)
	    dexdt2, deconvdt2 = dedtElastic(2, np.array(bunch_intensity_p2[i_step-1]), 1, np.array(Luminosity_invm2s[i_step-1]), betastar_m, sigmaElastic_m2)
            ex_norm_m1.append(ex_IBS_norm_m1+dexdt1*dt_s)
            ey_norm_m1.append(ey_norm_m1[i_step-1]*np.exp(-2*dt_s/tauSRxy_s)+dexdt1*dt_s);                                                                                                                                                                                
            ex_norm_m2.append(ex_IBS_norm_m2+dexdt2*dt_s)
            ey_norm_m2.append(ey_norm_m2[i_step-1]*np.exp(-2*dt_s/tauSRxy_s)+dexdt2*dt_s);     
	    
        elif emitBU=='EmpiricalBlowup':
            ex_norm_m1.append(ex_norm_m1[i_step-1]*np.exp(dt_s/tau_empirical_h1))
            ey_norm_m1.append(ey_norm_m1[i_step-1]*np.exp(dt_s/tau_empirical_v1))
            ex_norm_m2.append(ex_norm_m2[i_step-1]*np.exp(dt_s/tau_empirical_h2))
            ey_norm_m2.append(ey_norm_m2[i_step-1]*np.exp(dt_s/tau_empirical_v2))
        else:
            raise ValueError('Not understood!')


        bl_4sigma_s1.append(bl_IBS_4sigma_s1)
        bl_4sigma_s2.append(bl_IBS_4sigma_s2)


        if BoffON==1:
            #~ tauBOff_a = nb*Nb0[i-1]/(nIPs*sigmaBOff*L0[i-1]));#   % in minutes
            tauBOff_s1 = bunch_intensity_p1[i_step-1]/(sigmaBOff_m2*nIPs*(Luminosity_invm2s[i_step-1]))
            bunch_intensity_p1.append(bunch_intensity_p1[i_step-1]/(1.+dt_s/tauBOff_s1));
            tauBOff_s2 = bunch_intensity_p2[i_step-1]/(sigmaBOff_m2*nIPs*Luminosity_invm2s[i_step-1])
            bunch_intensity_p2.append(bunch_intensity_p2[i_step-1]/(1.+dt_s/tauBOff_s2));                
        else:
            bunch_intensity_p1.append(bunch_intensity_p1[i_step-1]);
            bunch_intensity_p2.append(bunch_intensity_p2[i_step-1]);

        Lumi_ATLAS = Lumi_inst(Number_bunches=1., bunch_intensity_p1=bunch_intensity_p1[i_step], bunch_intensity_p2=bunch_intensity_p2[i_step], ex_norm_m1=ex_norm_m1[i_step], ey_norm_m1=ey_norm_m1[i_step], ex_norm_m2=ex_norm_m2[i_step], ey_norm_m2=ey_norm_m2[i_step], bl_4sigma_s1=bl_4sigma_s1[i_step], bl_4sigma_s2=bl_4sigma_s2[i_step], betastar_m=betastar_m, phi_full_rad=phi_full_rad_ATLAS, gamma=gamma, experiment='ATLAS') 

        Lumi_CMS = Lumi_inst(Number_bunches=1., bunch_intensity_p1=bunch_intensity_p1[i_step], bunch_intensity_p2=bunch_intensity_p2[i_step], ex_norm_m1=ex_norm_m1[i_step], ey_norm_m1=ey_norm_m1[i_step], ex_norm_m2=ex_norm_m2[i_step], ey_norm_m2=ey_norm_m2[i_step], bl_4sigma_s1=bl_4sigma_s1[i_step], bl_4sigma_s2=bl_4sigma_s2[i_step], betastar_m=betastar_m, phi_full_rad=phi_full_rad_CMS, gamma=gamma, experiment='CMS')  

        Luminosity_invm2s.append((Lumi_ATLAS+Lumi_CMS)/2.)
        Luminosity_ATLAS_invm2s.append(Lumi_ATLAS)
        Luminosity_CMS_invm2s.append(Lumi_CMS)


    return tt_s, Luminosity_ATLAS_invm2s, Luminosity_CMS_invm2s, bunch_intensity_p1, bunch_intensity_p2, ex_norm_m1, ex_norm_m2, ey_norm_m1, ey_norm_m2, bl_4sigma_s1, bl_4sigma_s2




