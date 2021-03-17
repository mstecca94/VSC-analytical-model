function [Pc_comp,Ps_comp,Ptot_comp,Ptot]=TwoLC_freqSiC(U__DC,Vac,Sac,~,f__s)


Iac_rms=abs(Sac/(sqrt(3)*Vac));
I__ac=sqrt(2)*Iac_rms; %%% peak AC current

%% mosfet

U__Bt=600;
b__0t=0.0014953978290906684; 
b__1t=0.000023084336683723477;
b__2t=3.763236812059803e-8; 

rt=0.0042;

%% RMS currents

irms_m=abs(I__ac)/(2);

%% Conduction losses

Pc_t1=irms_m.^2*rt;

Pc_t2=Pc_t1;

Pc_d1=0;

Pc_d2=Pc_d1;

Pc_comp=[Pc_t1 Pc_d1 Pc_t2 Pc_d2];

Pc_tot=sum(Pc_comp);

%% Switching losses

Ps_t12=2*f__s*U__DC*(((1/2)*b__2t*I__ac^2+b__0t)*pi+2*b__1t*I__ac)/(4*pi*U__Bt);

Ps_d12=0;

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
