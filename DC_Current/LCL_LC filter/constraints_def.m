function [att,Ltot_min,L1_min,Cmax,Ltot_max,wres]=constraints_def(factor,qfact,irpp)

global vg fs Ig fc vdcmin vdcmax Pn Inom

Ltot_max=round(sqrt(vdcmin^2/8-vg^2)/(2*pi*fc*Ig),6);      % max total inductance / votage drop


Cmax=qfact*Pn/3/(vg^2*2*pi*fc);                         % max Cap / max Q absorption

L1_min=round(vg/(2*sqrt(6)*fs*Ig*irpp),6);              % min converter inductance / ripple

wres=factor*fs*2*pi;                                    % LCL resonance frequency

mh=(fs-2*fc);               % 1st V harmonic to mitigate (fs-2fc) 
wh=2*pi*mh;                 % 1st V harmonic to mitigate (fs-2fc) 

vh=2*vdcmin/(pi*(1-2*fc/fs))*besselj(2,pi*sqrt(2)*(1-2*fc/fs)*vg/vdcmin); % harmonic of fs-2fc / 
%first significant sideband
% vh calculated with vmin because it gives the maximum value

att=0.003*Ig/(vh/sqrt(2)); % required attenuation / 0.3% of Ig

Ltot_min1=wres^2/(wh*(wh^2-wres^2)*att);
Ltot_min2=(vdcmax/2-vg*0.05)/(fs*4*Inom); % ride through constraints

Ltot_min=round(max(Ltot_min1,Ltot_min2),6);


end