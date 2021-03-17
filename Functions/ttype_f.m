function [Pc_comp,Ps_comp,Ptot_comp,Ptot]=ttype_f(U__DC,Vac,Sac,ddeg,f__s)

d=ddeg/180.*pi;

Iac_rms=abs(Sac/(sqrt(3)*Vac));
I__ac=sqrt(2)*Iac_rms;

m=Vac*2*sqrt(2)/sqrt(3)/U__DC;

%% IGBT 1 parameters

U__Bt1=300;
b__0t1=0; 
b__1t1=9.1730e-05;
b__2t1=-1.5072e-08; 

rt1=0.005;
Vt1=0.70;

%% IGBT 2 parameters

U__Bt2=300;
b__0t2=0; 
b__1t2=7.9360e-05;
b__2t2=-3.7275e-08;

rt2=0.0031;
Vt2=0.82;

%% Diode 1 parameters

U__Bd1=300;
b__0d1=0;
b__1d1=0.000079362517886555;
b__2d1=-1.670978594570541e-7;

rd1=0.0042;
Vd1=0.9;


%% Diode 2 parameters


U__Bd2=300;
b__0d2=0;
b__1d2=0.00006452882219560093;
b__2d2=-1.1808213513239772e-7;

rd2=0.00179;
Vd2=0.85;


%% RMS currents

irms_t23d23=sqrt(3).*sqrt(-4.*m.*cos(d).^2+3.*pi-4.*m).*abs(I__ac)/(6.*sqrt(pi));

irms_t14=sqrt(2).*sqrt(3).*(cos(d)+1).*sqrt(m).*abs(I__ac)/(6.*sqrt(pi));

irms_d14=-sqrt(2).*sqrt(3).*(cos(d)-1).*sqrt(m).*abs(I__ac)/(6.*sqrt(pi));

%% AVG currents

iavg_t23d23=(m.*(abs(d)-(1/2).*pi).*cos(abs(d))-m.*sin(abs(d))+2).*I__ac/(2.*pi);

iavg_t14=-m.*I__ac.*((-pi+abs(d)).*cos(abs(d))-sin(abs(d)))/(4.*pi);

iavg_d14=m.*I__ac.*(-cos(abs(d)).*abs(d)+sin(abs(d)))/(4.*pi);


%% Conduction losses

Pc_t1=irms_t14.^2*rt1+iavg_t14*Vt1;

Pc_t2=irms_t23d23.^2*rt2+iavg_t23d23*Vt2;

Pc_t3=Pc_t2;

Pc_t4=Pc_t1;

Pc_d1=irms_d14.^2*rd1+iavg_d14*Vd1;

Pc_d2=irms_t23d23.^2*rd2+iavg_t23d23*Vd2;

Pc_d3=Pc_d2;

Pc_d4=Pc_d1;

Pc_comp=[Pc_t1 Pc_d1 Pc_t2 Pc_d2 Pc_t3 Pc_d3  Pc_t4 Pc_d4];

Pc_tot=sum(Pc_comp);

%% Switching losses


Ps_t14=-(-(1/2).*I__ac.*(b__2t1.*I__ac.*sin(abs(d))+2.*b__1t1).*cos(abs(d))+...
    (1/2).*b__2t1.*(-pi+abs(d)).*I__ac^2-b__1t1.*I__ac+b__0t1.*(-pi+abs(d))).*U__DC.*f__s/(4.*pi.*U__Bt1);

Ps_d23=-(-(1/2).*I__ac.*(b__2d2.*I__ac.*sin(abs(d))+2.*b__1d2).*cos(abs(d))+...
    (1/2).*b__2d2.*(-pi+abs(d)).*I__ac^2-b__1d2.*I__ac+b__0d2.*(-pi+abs(d))).*U__DC.*f__s/(4.*pi.*U__Bd2);

Ps_t23=-f__s.*U__DC.*(b__2t2.*I__ac^2.*cos(abs(d)).*sin(abs(d))-b__2t2.*I__ac^2.*abs(d)+...
    2.*b__1t2.*I__ac.*cos(abs(d))-2.*b__1t2.*I__ac-2.*b__0t2.*abs(d))/(8.*pi.*U__Bt2);

Ps_d14=-f__s.*U__DC.*(b__2d1.*I__ac^2.*cos(abs(d)).*sin(abs(d))-b__2d1.*I__ac^2.*abs(d)+...
    2.*b__1d1.*I__ac.*cos(abs(d))-2.*b__1d1.*I__ac-2.*b__0d1.*abs(d))/(8.*pi.*U__Bd1);


Ps_t1=Ps_t14;

Ps_t2=Ps_t23;

Ps_t3=Ps_t2;

Ps_t4=Ps_t1;

Ps_d1=Ps_d14;

Ps_d2=Ps_d23;

Ps_d3=Ps_d2;

Ps_d4=Ps_d1;

Ps_comp=[Ps_t1 Ps_d1 Ps_t2 Ps_d2 Ps_t3 Ps_d3  Ps_t4 Ps_d4];

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


Ptot_comp=[Pt1 Pd1 Pt2 Pd2 Pt3 Pd3  Pt4 Pd4 ];


Ptot=(Pc_tot+Ps_tot)*3;

end