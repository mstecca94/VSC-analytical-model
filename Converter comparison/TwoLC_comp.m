function [Pc_comp,Ps_comp,Ptot_comp,Ptot,I_DC]=TwoLC_comp(U__DC,Vac,Sac,ddeg,f__s)

d=ddeg/180.*pi;

Iac_rms=abs(Sac/(sqrt(3)*Vac));
I__ac=sqrt(2)*Iac_rms; %%% peak AC current

m=Vac*2*sqrt(2)/sqrt(3)/U__DC; %% modulation index

%% IGBT parameters

U__Bt=600;
b__0t=0; 
b__1t=1.7328e-04;
b__2t=4.8252e-08; 

rt=0.0050;
Vt=0.56;

%% Diode parameters

U__Bd=600;
b__0d=0;
b__1d=0.00011336152231240664;
b__2d=-1.1406623528819162e-7;

rd=0.0038;
Vd=0.50;

%% RMS currents

irms_t12=sqrt(2)*sqrt(3)*sqrt(8*m*cos(d)+3*pi)*abs(I__ac)/(12*sqrt(pi));

% sqrt(3)*sqrt((8*m-3*sin(abs(d))).*cos(abs(d))+3*abs(d))*abs(I__ac)/(6*sqrt(pi));

irms_d12=sqrt(2)*sqrt(3)*sqrt(-8*m*cos(d)+3*pi)*abs(I__ac)/(12*sqrt(pi));

% sqrt(3)*sqrt((-8*m+3*sin(abs(d))).*cos(abs(d))-3*abs(d)+3*pi)*abs(I__ac)/(6*sqrt(pi));


%% AVG currents

iavg_t12=I__ac*(pi*cos(d)*m+4)/(8*pi);

% I__ac.*(m.*cos(d).*pi-2.*cos(d)+2)/(4.*pi);

iavg_d12=-I__ac*(pi*cos(d)*m-4)/(8*pi);

% -I__ac.*(m.*cos(d).*pi-2.*cos(d)-2)/(4.*pi);


%% Conduction losses

Pc_t1=irms_t12.^2*rt+iavg_t12*Vt;

Pc_t2=Pc_t1;

Pc_d1=irms_d12.^2*rd+iavg_d12*Vd;

Pc_d2=Pc_d1;

Pc_comp=[Pc_t1 Pc_d1 Pc_t2 Pc_d2];

Pc_tot=sum(Pc_comp);

%% Switching losses

Ps_t12=2*f__s*U__DC*(((1/2)*b__2t*I__ac^2+b__0t)*pi+2*b__1t*I__ac)/(4*pi*U__Bt);

Ps_d12=2*f__s*U__DC*(((1/2)*b__2d*I__ac^2+b__0d)*pi+2*b__1d*I__ac)/(4*pi*U__Bd);

Ps_t1=Ps_t12;

Ps_t2=Ps_t1;

Ps_d1=Ps_d12;

Ps_d2=Ps_d1;

Ps_comp=[Ps_t1 Ps_d1 Ps_t2 Ps_d2];

Ps_tot=sum(Ps_comp);

%% total losses per component

Pt1=Pc_t1+Ps_t1;

Pt2=Pc_t2+Ps_t2;

Pd1=Pc_d1+Ps_d1;

Pd2=Pc_d2+Ps_d2;

Ptot_comp=[Pt1 Pd1 Pt2 Pd2];

Ptot=(Pc_tot+Ps_tot)*3;

%% DC Currents

I_dc_avg=3/4*I__ac*cos(d)*m;

I_dc_rms=sqrt(4*sqrt(3)/pi*m(1/4+(cos(d))^2));

I_hf_rms=sqrt(m*(sqrt(3)/(4*pi))*I__ac^2+(sqrt(3)/pi-9*m/16));

I_DC=[I_dc_avg I_dc_rms I_hf_rms];

