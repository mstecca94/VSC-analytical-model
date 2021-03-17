ccc
%%%%  LCL FILTER DESIGN  %%%%
%%%  Two level converter  %%%
%        Marco Stecca       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global vg fs Ig fc  vdcmin vdcmax Pn Inom Lmin

%% Set parameters
vg=400/sqrt(3);             % phase voltage
fc=50;                      % grid frequency
Pn=100000;                  % Phase nominal current
Ig=Pn/sqrt(3)/400;          % Phase nominal current
vdcmin=774;                 % DC voltage - minimum 
vdcmax=1004;                % DC voltage - minimum maximum
Inom=300;                   % seminconductor nominal current

Lmin=1e-6;                  % min inductance

Lb=(vg/(2*pi*fc*Ig));       % base inductance
Cb=(Ig/(2*pi*fc*vg));       % base capacitance

%% Variables for the filter design
fs=15000;                                               % switching frequency
factor=0.3;                                            % position of fres respect to fs
qfact=0.05;                                             % max Q absorption
irpp=0.5;                                               % max Peak to peak ripple / p.u.
%% define L/C constraints 

[att,Ltot_min,L1_min,Cmax,Ltot_max,wres]=constraints_def(factor,qfact,irpp);
Ltot_max_pu=Ltot_max/Lb;

if Ltot_min>Ltot_max
    
    disp('Design not possible')
    
end
%% LC filter 
    
% [L1lc,Clc,irppactlc]=filter_lc(L1_min,wres,Cmax,Lmin,Ltot_min);
% L1pulc=L1lc/Lb;
% Cpulc=Clc/Cb;

% figure;
% sgtitle('2 Level converter possible LC filter design')
% subplot(2,1,1)
% hold on
% plot(L1pulc)
% plot(Cpulc)
% legend('L','C')
% subplot(2,1,2)
% plot(irppactlc)

%% LCL filter  

% deisgn at minimum inductance

[L1lcl,L2lcl,Clcl,irppactlcl]=lcl_min_L(L1_min,wres,Ltot_min,Ltot_max,Cmax);


L1pulcl=L1lcl/Lb;
L2pulcl=L2lcl/Lb;
Cpulcl=Clcl/Cb;

figure;
sgtitle('2 Level converter possible LCL filter design')
subplot(2,1,1)
hold on
plot(L1pulcl)
plot(L2pulcl)
plot(Cpulcl)
legend('L_{c}','L_{g}','C')
subplot(2,1,2)
plot(irppactlcl)

