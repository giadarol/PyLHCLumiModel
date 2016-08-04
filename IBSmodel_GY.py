#!/usr/bin/python

import cmath
from math import *
import numpy 
from scipy.constants import c

ene_tol_GeV = 20.

def IBSmodel(IBSON, gamma, bunch_intensity_p,ex_norm_m,ey_norm_m,bl_4sigma_s,VRF_V,dt_s):
    
    En = gamma*0.93828
    Nb = bunch_intensity_p*1e-11
    ex0 = ex_norm_m*1e6
    ey0 = ey_norm_m*1e6
    bl0 = bl_4sigma_s*1e9
    V0 = VRF_V/1e6
    t = dt_s/60.
    
    if IBSON==True or IBSON==1:
    
        #print 'Energy', En
        if fabs(En-6500)<ene_tol_GeV:
            
            if V0==12:
                a00=5.8128e-06*t**1.9244;
                a01=0.033455*t**0.71595-4.3735;
                a10=0.0033935*t**0.92978;
                a11=0.0069376*t**0.8735-1.6269;
                b00=1.6958e-07*t**1.9284+1.1511e-05;
                b01=0.069374*t**0.69181-4.6381;
                b10=1.4673e-07*t**2-5.5269e-05*t+5.0618e-05;
                b11=-0.00028513*t**1.7199+2.6338;
                b20=-0.0005309*t**0.99039+1.0001;
                b21=1.6569e-07*t**2-5.2329e-06*t+5.223e-05;
            
                a00l=2.6454e-06*t**2+1.925e-05*t-8.4302e-05;
                a01l=11.6706*t**-0.049427-10.9141;
                a02l=-7.5066e-08*t**2+1.2167e-05*t+1.9744e-05;
                a10l=0.0012145*t**0.92003;
                a11l=0.0062127*t**1.1036-1.3455;
                a12l=-4.8327e-07*t**2.3482;
                b00l=7.7846e-07*t**2+3.5015e-05*t-4.6149e-05;
                b01l=0.00014708*t**2-0.019588*t-0.70599;
                b02l=-9.6701e-09*t**2+1.8945e-06*t-1.6834e-05;
                b10l=-4.9054e-07*t**2+0.00086599*t+0.55636;
                b11l=2.5997e-07*t**2+8.1857e-06*t-7.235e-05;
                b12l=4.1863e-07*t**2-0.0014262*t+0.44368;
            elif V0==10:
                a00=8.1034e-06*t**1.8783;
                a01=0.28811*t**0.36532-4.6617;
                a10=0.0037095*t**0.90923;
                a11=0.020951*t**0.72569-1.5621;
                b00=2.6472e-07*t**1.8488+1.4178e-05;
                b01=0.042724*t**0.85404-4.1085;
                b10=2.2905e-07*t**2+-6.3663e-05*t+6.6759e-05;
                b11=-0.0012795*t**1.5417+2.2332;
                b20=-0.00051568*t**0.99814+0.99999;
                b21=4.2472e-07*t**2-9.6913e-06*t+6.8151e-05;
            
                a00l=3.8677e-06*t**2+0.00014648*t-0.0013837;
                a01l=14.9302*t**-1.3934-1.9835;
                a02l=-4.9863e-08*t**2+3.6098e-06*t+1.3622e-05;
                a10l=0.001452*t**0.8945;
                a11l=0.090579*t**0.61505-1.6815;
                a12l=-8.2446e-07*t**2.3908;
                b00l=1.851e-06*t**2+0.00013284*t-0.00057515;
                b01l=0.00021686*t**2+-0.022711*t-0.86553;
                b02l=-5.0828e-08*t**2+4.1981e-06*t-3.0944e-05;
                b10l=-4.4717e-07*t**2+0.0018648*t+0.54738;
                b11l=2.4981e-07*t**2+4.9872e-05*t-0.00040738;
                b12l=1.6328e-07*t**2+-0.0024397*t+0.45289;   
            elif V0==8:
                a00=1.0093e-05*t**1.8717;
                a01=-6.876*t**-0.46105+-2.3664;
                a10=0.0040721*t**0.89564;
                a11=0.019936*t**0.74094+-1.5277;
                b00=4.5268e-07*t**1.7642+1.8187e-05;
                b01=0.093895*t**0.67778+-3.972;
                b10=3.0534e-07*t**2+-7.5303e-05*t+7.9589e-05;
                b11=-0.0013357*t**1.5252+2.0487;
                b20=-1.2699e-08*t**2+-0.00050988*t+0.99997;
                b21=4.6533e-07*t**2+-1.3315e-05*t+7.2448e-05;
            else:
                raise ValueError('Voltage not available!')
                
            a0l=a00l*ex0**a01l+a02l;
            a1l=a10l*ex0**a11l+a12l;
            b0l=b00l*ex0**b01l+b02l;
            b1l=b10l*ex0**b11l+b12l;

            C0l=a1l*Nb+a0l;
            ccSRl=b0l*Nb+b1l;

            ey=ey0

            IBSl=C0l*bl0**-3.3+ccSRl;
            bl=IBSl*bl0;

            a0=a00*bl0**a01;
            a1=a10*bl0**a11;
            b0=b00*bl0**b01;
            b1=b10*bl0**b11;
            b2=b20*bl0**b21;

            C0=a1*Nb+a0;
            ccSR=b0*Nb**2+b1*Nb+b2;

            IBSx=ccSR+C0/ex0**2;
            ex=IBSx*ex0;
        
        elif  fabs(En-450)<ene_tol_GeV:
            raise ValueError('Injection not implemented yet')
            
        else:
            raise ValueError('Energy not available!')



        
        ex_out_norm_m = ex*1e-6
        ey_out_norm_m = ey*1e-6
        bl_out_4sigma_s =  bl*1e-9
    
    else:
        ex_out_norm_m = ex_norm_m
        IBSx = 0.
        bl_out_4sigma_s = bl_4sigma_s
        IBSl = 0.
        ey_out_norm_m = ey_norm_m

    return ex_out_norm_m, IBSx, bl_out_4sigma_s, IBSl, ey_out_norm_m 

#~ test=IBSmodel(6500,1,3.3,3.3,1.25,12,10);
#~ print 'test=',test
	
	
