function [Pc_comp,Ps_comp,Ptot_comp,Ptot]=npcc_freq(U__DC,Vac,Sac,ddeg,f__s)

d=ddeg/180.*pi;

Iac_rms=abs(Sac/(sqrt(3)*Vac));

I__ac=sqrt(2)*Iac_rms;

m=Vac*2*sqrt(2)/sqrt(3)/U__DC;

%% IGBT 1 parameters

U__Bt1=300;
b__0t1=0; 
b__1t1=6.8375e-05 ;
b__2t1=-1.5335e-08; 

rt=0.0031;
Vt=0.82;

%% IGBT 2 parameters

U__Bt2=300;
b__0t2=0; 
b__1t2=6.1203e-05;
b__2t2=-1.3600e-08; 

rt=0.0031;
Vt=0.82;


%% Diode 1-2 parameters

U__Bd12=300;
b__0d12=0;
b__1d12=0.000040843026433002205;
b__2d12=-3.930898280888438e-8;

rd12=0.0027;
Vd12=0.85;

%% Diode 5 parameters

U__Bd5=300;
b__0d5=0;
b__1d5=0.000023833567117242575;
b__2d5=-2.610058019068224e-8;

rd5=0.0027;
Vd5=0.85;

%% RMS currents

irms_t23=sqrt(3).*sqrt(-2.*m.*cos(d).^2+4.*m.*cos(d)+3.*pi-2.*m).*abs(I__ac)/(6.*sqrt(pi));

irms_t14=sqrt(2).*sqrt(3).*(cos(d)+1).*sqrt(m).*abs(I__ac)/(6.*sqrt(pi));

irms_d1234=-sqrt(2).*sqrt(3).*(cos(d)-1).*sqrt(m).*abs(I__ac)/(6.*sqrt(pi));

irms_d56=sqrt(3).*sqrt(-(2.*(m.*cos(2.*d)+3.*m-3.*pi.*(1/2))).*I__ac^2)/(6.*sqrt(pi));

%% AVG currents

iavg_t23=-I__ac.*(-m.*cos(abs(d)).*abs(d)+m.*sin(abs(d))-4)/(4.*pi);

iavg_t14=-I__ac.*((-pi+abs(d)).*cos(abs(d))-sin(abs(d))).*m/(4.*pi);

iavg_d1234=m.*I__ac.*(-cos(abs(d)).*abs(d)+sin(abs(d)))/(4.*pi);

iavg_d56=I__ac.*(m.*(abs(d)-(1/2).*pi).*cos(abs(d))-m.*sin(abs(d))+2)/(2.*pi);

%% Conduction losses

Pc_t1=irms_t14.^2*rt+iavg_t14*Vt;

Pc_t2=irms_t23.^2*rt+iavg_t23*Vt;

Pc_t3=Pc_t2;

Pc_t4=Pc_t1;

Pc_d1=irms_d1234.^2*rd12+iavg_d1234*Vd12;

Pc_d2=irms_d1234.^2*rd12+iavg_d1234*Vd12;

Pc_d3=Pc_d2;

Pc_d4=Pc_d1;

Pc_d5=irms_d56.^2*rd5+iavg_d56*Vd5;

Pc_d6=Pc_d5;

Pc_comp=[Pc_t1 Pc_d1 Pc_t2 Pc_d2 Pc_t3 Pc_d3  Pc_t4 Pc_d4 Pc_d5  Pc_d6];

Pc_tot=sum(Pc_comp);

%% Switching losses

Ps_t14=-(-(1/2).*I__ac.*(b__2t1.*I__ac.*sin(abs(d))+2.*b__1t1).*cos(abs(d))+...
    (1/2).*b__2t1.*(-pi+abs(d)).*I__ac^2-b__1t1.*I__ac+b__0t1.*(-pi+abs(d))).*U__DC.*f__s/(4.*pi.*U__Bt1);

Ps_t23=-f__s.*U__DC.*(b__2t2.*I__ac^2.*cos(abs(d)).*sin(abs(d))-b__2t2.*I__ac^2.*abs(d)+...
    2.*b__1t2.*I__ac.*cos(abs(d))-2.*b__1t2.*I__ac-2.*b__0t2.*abs(d))/(8.*pi.*U__Bt2);

Ps_d14=-f__s.*U__DC.*(b__2d12.*I__ac^2.*cos(abs(d)).*sin(abs(d))-b__2d12.*I__ac^2.*abs(d)+...
    2.*b__1d12.*I__ac.*cos(abs(d))-2.*b__1d12.*I__ac-2.*b__0d12.*abs(d))/(8.*pi.*U__Bd12);

Ps_d56=-(-(1/2).*I__ac.*(b__2d5.*I__ac.*sin(abs(d))+2.*b__1d5).*cos(abs(d))+...
    (1/2).*b__2d5.*(-pi+abs(d)).*I__ac^2-b__1d5.*I__ac+b__0d5.*(-pi+abs(d))).*U__DC.*f__s/(4.*pi.*U__Bd5);

Ps_t1=Ps_t14;

Ps_t2=Ps_t23;

Ps_t3=Ps_t2;

Ps_t4=Ps_t1;

Ps_d1=Ps_d14;

Ps_d2=0;

Ps_d3=0;

Ps_d4=Ps_d1;

Ps_d5=Ps_d56;

Ps_d6=Ps_d5;

Ps_comp=[Ps_t1 Ps_d1 Ps_t2 Ps_d2 Ps_t3 Ps_d3  Ps_t4 Ps_d4 Ps_d5  Ps_d6];

Ps_tot=sum(Ps_comp);

%% total losses per component

Pt1=Pc_t1+Ps_t1;

Pt2=Pc_t2+Ps_t2;

Pt3=Pc_t3+Ps_t3;

Pt4=Pc_t4+Ps_t4;

Pd1=Pc_d1+Ps_d1;

Pd2=Pc_d2+Ps_d2;

Pd3=Pc_d3+Ps_d3;

Pd4=Pc_d4+Ps_d4;

Pd5=Pc_d5+Ps_d5;

Pd6=Pc_d6+Ps_d6;

Ptot_comp=[Pt1 Pd1 Pt2 Pd2 Pt3 Pd3  Pt4 Pd4 Pd5  Pd6];

Ptot=(Pc_tot+Ps_tot)*3;

end