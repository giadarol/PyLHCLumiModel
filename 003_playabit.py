import numpy as np
import pylab as pl
from LumiModel_GY import LumiModel
from scipy.integrate import cumtrapz
import mystyle as ms


gamma = 6927.64 #6.5TeV
betastar_m = .4

bunch_intensity = 1.11e11
#n_emittx_init_m = 2.35e-6
#n_emitty_init_m = 1.8e-6
n_emittx_init_m = 1e-6*(2.1+1.4)/2.
n_emitty_init_m = 1e-6*(2.1+1.4)/2.
bunch_length_init_ns = 1.07e-9#4 sigma
tFill_s  = 20*3600.
VRF_V = 12e6
tauSRxy_s = 64.7*3600; 
tauSRl_s = 32.35*3600; 
sigmaBOff_m2 = 110*1e-31;
emitBU="EmpiricalBlowup"
tau_empirical_h = np.inf
tau_empirical_v = np.inf
time_t_around_h = 6 #3.5 #5.7

N_bunches = 2208

fontsz =14
                    
                    
#full_Xing_angle_list = [2*185e-6, 2*142.5e-6]
full_Xing_angle_list = [2*185e-6]
color_list = ['b', 'g']

na = np.array

sp=None
pl.close('all')
ms.mystyle(fontsz =fontsz)

fig_h = pl.figure(2, figsize=(16,8))
fig_h.patch.set_facecolor('w')

spint = pl.subplot(2,3,1)
spbl = pl.subplot(2,3,2, sharex=spint)
spemi = pl.subplot(2,3,3, sharex=spint)
splumi = pl.subplot(2,3,4, sharex=spint)
spinteg = pl.subplot(2,3,5, sharex=spint)
spweek = pl.subplot(2,3,6, sharex=spint)


for ii in xrange(len(full_Xing_angle_list)):
    
    full_Xing_angle = full_Xing_angle_list[ii]
    col = color_list[ii]
                    
    tt_s, Luminosity_invm2s, bunch_intensity_p, ex_norm_m, ey_norm_m, bl_4sigma_s = LumiModel(gamma=gamma, betastar_m=betastar_m, phi_full_rad=full_Xing_angle, 
                bunch_intensityin_p=bunch_intensity, exin_norm_m=n_emittx_init_m, eyin_norm_m=n_emitty_init_m, 
                blin_4sigma_s=bunch_length_init_ns, tFill_s=tFill_s,
                tauSRxy_s=tauSRxy_s, tauSRl_s=tauSRl_s, 
                sigmaBOff_m2=sigmaBOff_m2, VRF_V=VRF_V,
                IBSON=1, emitBU=emitBU, BoffON=1, nIPs = 2., dt_s=15*60., tau_empirical_h=tau_empirical_h, tau_empirical_v=tau_empirical_v)



    integrated_luminosity_inv_fb = cumtrapz(y=Luminosity_invm2s, x=tt_s)*1e-43
    
    L_week_full_beam = na(integrated_luminosity_inv_fb)*N_bunches/(tt_s[:-1]+time_t_around_h*3600)*24*3600*7
    
    i_max = np.argmax(L_week_full_beam)
    
    
    
    print np.sum(integrated_luminosity_inv_fb)

    Luminosity_invcm2s = na(Luminosity_invm2s)*1e-4


    
    spint.plot(tt_s/3600., na(bunch_intensity_p), linewidth=2.)
    spint.ticklabel_format(style='sci', scilimits=(0,0),axis='y') 
    spint.set_ylabel('Bunch intensity [p]')
    spint.set_xlabel('Time [h]')
    spint.grid('on')


    spbl.plot(tt_s/3600., na(bl_4sigma_s), linewidth=2.)
    spbl.set_ylabel('Bunch length (4 sigmas) [s]')
    spbl.set_xlabel('Time [h]')
    spbl.ticklabel_format(style='sci', scilimits=(0,0),axis='y') 
    spbl.grid('on')


    spemi.plot(tt_s/3600., na(ex_norm_m), label='H', linewidth=2.)
    spemi.plot(tt_s/3600., na(ey_norm_m), label='V', linewidth=2.)
    spemi.set_ylabel('Norm emittance [m]')
    spemi.set_xlabel('Time [h]')
    spemi.ticklabel_format(style='sci', scilimits=(0,0),axis='y') 
    spemi.grid('on')
    #~ spemi.legend(loc='best')
    spemi.set_ylim(1.e-6, 3.5e-6)



    splumi.plot(tt_s/3600., na(Luminosity_invcm2s)*N_bunches, linewidth=2.)
    splumi.set_ylabel('Luminosity [cm^-2 s^-1]')
    splumi.ticklabel_format(style='sci', scilimits=(0,0),axis='y') 
    splumi.set_xlabel('Time [h]')
    splumi.grid('on')

    spinteg.plot(tt_s[:-1]/3600., na(integrated_luminosity_inv_fb)*N_bunches, linewidth=2.)
    spinteg.set_ylabel('Integr. luminosity [fb^-1]')
    spinteg.set_xlabel('Time [h]')
    spinteg.grid('on')

    spweek.plot(tt_s[:-1]/3600., L_week_full_beam, linewidth=2., 
    label='t_opt=%.1f, L_week_opt=%.1f'%(tt_s[i_max]/3600., L_week_full_beam[i_max]))
    spweek.set_ylabel('Integr. luminosity per week [fb^-1]')
    spweek.set_xlabel('Time [h]')
    spweek.grid('on')

spweek.legend(loc='best', prop={'size':fontsz})
pl.subplots_adjust(left=.09, right=.98, wspace=.3)
pl.show()
