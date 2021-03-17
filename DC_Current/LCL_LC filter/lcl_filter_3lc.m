ccc
%%%%  LCL FILTER DESIGN  %%%%
%%%  Three level converter  %%%
%        Marco Stecca       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set parameters
vg=400/sqrt(3);             % phase voltage
fc=50;                      % grid frequency
Pn=100000;                  % Phase nominal current
Ig=Pn/sqrt(3)/400;          % Phase nominal current
vdcmin=774;                 % DC voltage - minimum 
vdcmax=1004;                % DC voltage - minimum maximum
Inom=300;                   % seminconductor nominal current

Ltot_max=round(sqrt(vdcmin^2/8-vg^2)/(2*pi*fc*Ig),6);      % max total inductance / votage drop

Lmin=1e-6;                  % min inductance

Lb=(vg/(2*pi*fc*Ig));       % base inductance
Cb=(Ig/(2*pi*fc*vg));       % base capacitance

%% Variables for the filter design

fs=15000;                                               % switching frequency
wfactor=0.1;                                            % position of fres respect to fs
qfact=0.05;                                             % max Q absorption
irpp=0.5;                                               % max Peak to peak ripple / p.u.

wres=wfactor*fs*2*pi;                                   % LCL resonance frequency

Cmax=qfact*Pn/3/(vg^2*2*pi*fc);                         % max Cap / max Q absorption

L1_min=round(vdcmin/(6*fs*Ig*irpp),6);              	% min converter inductance

%% harmonic calculation according to fs
mh=(fs-2*fc);               % 1st V harmonic to mitigate (fs-2fc) 
wh=2*pi*mh;                 % 1st V harmonic to mitigate (fs-2fc) 

for k=1:1:1000000
vh2(k)=4*vdcmin/(pi^2*(2*(1-2*fc/fs)-1))*(2*k-1)*besselj(2*k-1,1*pi*sqrt(2)*(1-...
    2*fc/fs)*vg/vdcmin)/(2*k-1-4)/(2*k-1+4); % harmonic of fs-2fc / 
end

vh3=sum(vh2);
vh=vdcmin/(pi*2*(1-2*fc/fs))*besselj(-3,4*pi*sqrt(2)*(1-2*fc/fs)*vg/vdcmin); % harmonic of fs-2fc /
vhtot=abs(vh3+vh);

%first significant sideband
% vh calculated with vmin because it gives the maximum value

att=0.003*Ig/(vhtot/sqrt(2)); % required attenuation / 0.3% of Ig

Ltot_min1=wres^2/(wh*(wh^2-wres^2)*att);
Ltot_min2=(vdcmax/2-vg*0.05)/(fs*4*Inom);

Ltot_min=max(Ltot_min1,Ltot_min2);
%%  check if LC is enough

if Ltot_min<=L1_min % LC filter is enough
    
    i=0;
    
    for t=L1_min:0.0000001:Ltot_max
        i=i+1;

        L1(i)=t;

        C(i)=1/(wres^2*L1(i));

        irppact(i)=vg/(2*sqrt(6)*fs*L1(i))/Ig;
        
        L2(i)=0;

        if L1(i)<=Lmin || C(i)>=Cmax
            break
        end
    end
    
else % i need LCL filter  

    i=0;
    for t=L1_min:0.0000001:Ltot_max
        i=i+1;

        L1(i)=t;

        L2(i)=Ltot_min1-L1(i);

        C(i)=(L1(i)+L2(i))/(L1(i)*L2(i)*wres^2);

        irppact(i)=vg/(2*sqrt(6)*fs*L1(i))/Ig;

        if L1(i)<=Lmin || L2(i)<=Lmin || C(i)>=Cmax
            break
        end
    end
    
end

L1pu=L1/Lb;
L2pu=L2/Lb;
Ltot_max_pu=Ltot_max/Lb;

Cpu=C/Cb;


figure;
hold on
plot(L1pu)
plot(L2pu)
plot(Cpu)
legend('L_{c}','L_{g}','C')
title('3 Level converter possible filter design')

