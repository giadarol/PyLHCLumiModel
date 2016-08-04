import numpy as np
import pylab as pl
from LumiModel_GY import LumiModel

gamma = 6927.64 #6.5TeV
betastar_m = .4
full_Xing_angle=2*185e-6
bunch_intensity = 1.15e11
n_emittx_init_m = 2.5e-6
n_emitty_init_m = 2.e-6
bunch_length_init_ns = 1e-9#4 sigma
tFill_s  = 20*3600.

VRF_V = 12e6

tauSRxy_s = 64.7*3600; 
tauSRl_s = 32.35*3600; 

sigmaBOff_m2 = 110*1e-31;

N_bunches = 2064

#~ LumiModel(gamma=gamma, betastar_m=betastar_m, phi_full_rad=full_Xing_angle, bunch_intensityin_p=bunch_intensity, 
                #~ exin_norm_m=n_emittx_init_m, eyin_norm_m=n_emitty_init_m, 
                #~ blin_4sigma_s=bunch_length_init_ns, tFill_s=tFill_s,
                #~ tauh_s=None, tauv_s=None, tauSRxy_s=tauSRxy_s, tauSRl_s=tauSRl_s, 
                #~ sigmaBOff_m2=sigmaBOff_m2, VRF_V=VRF_V, 
                #~ IBSON=1, emitBU="Model", BoffON=1, nIPs = 2., dt_s=15*60.)
                    
tt_s, Luminosity_invm2s, bunch_intensity_p, ex_norm_m, ey_norm_m, bl_4sigma_s = LumiModel(gamma=gamma, betastar_m=betastar_m, phi_full_rad=full_Xing_angle, 
            bunch_intensityin_p=bunch_intensity, exin_norm_m=n_emittx_init_m, eyin_norm_m=n_emitty_init_m, 
            blin_4sigma_s=bunch_length_init_ns, tFill_s=tFill_s,
            tauSRxy_s=tauSRxy_s, tauSRl_s=tauSRl_s, 
            sigmaBOff_m2=sigmaBOff_m2, VRF_V=VRF_V,
            IBSON=1, emitBU="Model", BoffON=1, nIPs = 2., dt_s=15*60.)

na = np.array

Luminosity_invcm2s = na(Luminosity_invm2s)*1e-4

pl.close('all')
pl.figure(2, figsize=(10,8))
sp = pl.subplot(2,2,1)
pl.plot(tt_s/3600., na(Luminosity_invcm2s)*N_bunches)
pl.ylabel('Luminosity [cm^-2 s^-1]')
pl.gca().ticklabel_format(style='sci', scilimits=(0,0),axis='y') 
pl.xlabel('Time [h]')
pl.grid()

pl.subplot(2,2,2, sharex=sp)
pl.plot(tt_s/3600., na(bunch_intensity_p))
pl.gca().ticklabel_format(style='sci', scilimits=(0,0),axis='y') 
pl.ylabel('Bunch intensity [p]')
pl.xlabel('Time [h]')
pl.grid()

pl.subplot(2,2,3, sharex=sp)
pl.plot(tt_s/3600., na(bl_4sigma_s))
pl.ylabel('Bunch length (4 sigmas) [s]')
pl.xlabel('Time [h]')
pl.gca().ticklabel_format(style='sci', scilimits=(0,0),axis='y') 
pl.grid()

pl.subplot(2,2,4, sharex=sp)
pl.plot(tt_s/3600., na(ex_norm_m), label='H')
pl.plot(tt_s/3600., na(ey_norm_m), label='V')
pl.ylabel('Norm emittance [m]')
pl.xlabel('Time [h]')
pl.gca().ticklabel_format(style='sci', scilimits=(0,0),axis='y') 
pl.grid()
pl.legend(loc='best')

pl.show()
