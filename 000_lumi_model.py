import numpy as np
import pylab as pl
from LumiModel_GY import LumiModel
from scipy.integrate import cumtrapz

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
#~ emitBU = "Model"


gamma = 6927.64 #6.5TeV
betastar_m = .4
full_Xing_angle= 2*185e-6
bunch_intensity = 1.1e11
n_emittx_init_m = 2e-6
n_emitty_init_m = 2e-6
bunch_length_init_ns = 1.2e-9#4 sigma
tFill_s  = 40*3600.
VRF_V = 12e6
tauSRxy_s = 64.7*3600; 
tauSRl_s = 32.35*3600; 
sigmaBOff_m2 = 110*1e-31;
emitBU="EmpiricalBlowup"
tau_empirical_h = np.inf
tau_empirical_v = np.inf
time_t_around_h = 5.7

N_bunches = 2064
                    
tt_s, Luminosity_invm2s, bunch_intensity_p, ex_norm_m, ey_norm_m, bl_4sigma_s = LumiModel(gamma=gamma, betastar_m=betastar_m, phi_full_rad=full_Xing_angle, 
            bunch_intensityin_p=bunch_intensity, exin_norm_m=n_emittx_init_m, eyin_norm_m=n_emitty_init_m, 
            blin_4sigma_s=bunch_length_init_ns, tFill_s=tFill_s,
            tauSRxy_s=tauSRxy_s, tauSRl_s=tauSRl_s, 
            sigmaBOff_m2=sigmaBOff_m2, VRF_V=VRF_V,
            IBSON=1, emitBU=emitBU, BoffON=1, nIPs = 2., dt_s=15*60., tau_empirical_h=tau_empirical_h, tau_empirical_v=tau_empirical_v)

na = np.array

integrated_luminosity_inv_fb = cumtrapz(y=Luminosity_invm2s, x=tt_s)*1e-43

Luminosity_invcm2s = na(Luminosity_invm2s)*1e-4
sp=None
pl.close('all')
fig_h = pl.figure(2, figsize=(16,8))
fig_h.patch.set_facecolor('w')
sp = pl.subplot(2,3,1, sharex=sp)
pl.plot(tt_s/3600., na(bunch_intensity_p), linewidth=2.)
pl.gca().ticklabel_format(style='sci', scilimits=(0,0),axis='y') 
pl.ylabel('Bunch intensity [p]')
pl.xlabel('Time [h]')
pl.grid()

sp = pl.subplot(2,3,2, sharex=sp)
pl.plot(tt_s/3600., na(bl_4sigma_s), linewidth=2.)
pl.ylabel('Bunch length (4 sigmas) [s]')
pl.xlabel('Time [h]')
pl.gca().ticklabel_format(style='sci', scilimits=(0,0),axis='y') 
pl.grid()

sp = pl.subplot(2,3,3, sharex=sp)
pl.plot(tt_s/3600., na(ex_norm_m), label='H', linewidth=2.)
pl.plot(tt_s/3600., na(ey_norm_m), label='V', linewidth=2.)
pl.ylabel('Norm emittance [m]')
pl.xlabel('Time [h]')
pl.gca().ticklabel_format(style='sci', scilimits=(0,0),axis='y') 
pl.grid()
pl.legend(loc='best')
pl.ylim(1.e-6, 3.5e-6)


sp = pl.subplot(2,3,4, sharex=sp)
pl.plot(tt_s/3600., na(Luminosity_invcm2s)*N_bunches, linewidth=2.)
pl.ylabel('Luminosity [cm^-2 s^-1]')
pl.gca().ticklabel_format(style='sci', scilimits=(0,0),axis='y') 
pl.xlabel('Time [h]')
pl.grid()

sp = pl.subplot(2,3,5, sharex=sp)
pl.plot(tt_s[:-1]/3600., na(integrated_luminosity_inv_fb)*N_bunches, linewidth=2.)
pl.ylabel('Integr. luminosity per week [fb^-1]')
pl.xlabel('Time [h]')
pl.grid()

sp = pl.subplot(2,3,6, sharex=sp)
pl.plot(tt_s[:-1]/3600., na(integrated_luminosity_inv_fb)*N_bunches/(tt_s[:-1]+time_t_around_h*3600)*24*3600*7, linewidth=2.)
pl.ylabel('Integr. luminosity [fb^-1]')
pl.xlabel('Time [h]')
pl.grid()


pl.subplots_adjust(left=.09, right=.98, wspace=.3)
pl.show()
