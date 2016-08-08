import numpy as np
import pylab as pl
from LumiModel_BU_GY import LumiModel

#~ gamma = 6927.64 #6.5TeV
#~ betastar_m = .4
#~ full_Xing_angle=2*185e-6
#~ bunch_intensity = 1.15e8
#~ n_emittx_init_m = 2.5e-6
#~ n_emitty_init_m = 2.e-6
#~ bunch_length_init_ns = 1e-9*100#4 sigma
#~ tFill_s  = 20*3600.
#~ VRF_V = 12e6
#~ tauSRxy_s = np.inf #64.7*3600; 
#~ tauSRl_s = np.inf  #32.35*3600; 
#~ sigmaBOff_m2 = 110*1e-31;


gamma = 6927.64 #6.5TeV
betastar_m = .4
full_Xing_angle= 2*145e-6#2*185e-6
bunch_intensity = 1.1e11
n_emittx_init_m = 3e-6
n_emitty_init_m = 3e-6
bunch_length_init_ns = 1.1e-9#4 sigma
tFill_s  = 20*3600.
VRF_V = 12e6
tauSRxy_s = 64.7*3600; 
tauSRl_s = 32.35*3600; 
sigmaBOff_m2 = 110*1e-31;

N_bunches = 2064
                    
tt_s, Luminosity_invm2s, bunch_intensity_p, ex_norm_m, ey_norm_m, bl_4sigma_s = LumiModel(gamma=gamma, betastar_m=betastar_m, phi_full_rad=full_Xing_angle, 
            bunch_intensityin_p=bunch_intensity, exin_norm_m=n_emittx_init_m, eyin_norm_m=n_emitty_init_m, 
            blin_4sigma_s=bunch_length_init_ns, tFill_s=tFill_s,
            tauSRxy_s=tauSRxy_s, tauSRl_s=tauSRl_s, 
            sigmaBOff_m2=sigmaBOff_m2, VRF_V=VRF_V,
            IBSON=1, emitBU="EmpiricalBlowup", BoffON=1, nIPs = 2., dt_s=15*60.)

na = np.array

Luminosity_invcm2s = na(Luminosity_invm2s)*1e-4

Luminositytot_fb = sum(na(Luminosity_invm2s)*15*60*N_bunches/10**43)
#print Luminositytot_fb

pl.close('all')
fig_h = pl.figure(2, figsize=(10,8))
fig_h.patch.set_facecolor('w')
sp = pl.subplot(2,2,1)
pl.plot(tt_s/3600., na(Luminosity_invcm2s)*N_bunches, linewidth=2.,label='Lint =  %.2f fb^-1'%(Luminositytot_fb))
pl.ylabel('Luminosity [cm^-2 s^-1]')
pl.gca().ticklabel_format(style='sci', scilimits=(0,0),axis='y') 
pl.xlabel('Time [h]')
pl.legend(loc='best')
pl.grid()

pl.subplot(2,2,2, sharex=sp)
pl.plot(tt_s/3600., na(bunch_intensity_p), linewidth=2.)
pl.gca().ticklabel_format(style='sci', scilimits=(0,0),axis='y') 
pl.ylabel('Bunch intensity [p]')
pl.xlabel('Time [h]')
pl.grid()

pl.subplot(2,2,3, sharex=sp)
pl.plot(tt_s/3600., na(bl_4sigma_s), linewidth=2.)
pl.ylabel('Bunch length (4 sigmas) [s]')
pl.xlabel('Time [h]')
pl.gca().ticklabel_format(style='sci', scilimits=(0,0),axis='y') 
pl.grid()

pl.subplot(2,2,4, sharex=sp)
pl.plot(tt_s/3600., na(ex_norm_m), label='H', linewidth=2.)
pl.plot(tt_s/3600., na(ey_norm_m), label='V', linewidth=2.)
pl.ylabel('Norm emittance [m]')
pl.xlabel('Time [h]')
pl.gca().ticklabel_format(style='sci', scilimits=(0,0),axis='y') 
pl.grid()
pl.legend(loc='best')
pl.ylim(1.e-6, 3.e-6)

pl.show()
