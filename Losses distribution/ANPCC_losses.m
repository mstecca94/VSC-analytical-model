%% Design parameters
f__s=20000;
U__DC=900;
Vac=400;
Pac=100000; %% 200 kW
Iac_rms=Pac/(sqrt(3)*Vac);

ddeg=0;
% ddeg=-180:1:180;
d=ddeg/180.*pi;
I__ac=sqrt(2)*Iac_rms; %%% peak AC current

m=Vac*2*sqrt(2)/sqrt(3)/U__DC; % modulation index

t__a=0.5; %%% I divide the zero path in the two cases 
% t__a indicates the zero sequences P <--> 0U1/0U2 and N <--> 0L1/0L2
% t__a is the time when the zero voltage current flows in the same
% direction of P or N
t__b=1-t__a;  %%%  t__b indicates the zero sequences opposite to P or N / t__b=1-t__a
U__s=1; %%% time p.u. on sequence N <--> 0U1 / N <--> 0U2 is (1-U__s)
L__s=1; %%% time p.u. on sequence N <--> 0L1 / N <--> 0L2 is (1-L__s)

%%% for optimal losses sharing t__a=0.5 / U__s=1 / L__s=1

%% IGBT 1 parameters

U__Bt=300;
b__0t=0; 
b__1t=6.1203e-05;
b__2t=-1.3600e-08; 

rt=0.0038;
Vt=0.735;

% % % % % %%% IGBT 2 parameters
% % % % % 
% % % % % U__Bt=600;
% % % % % b__0t=0; 
% % % % % b__1t=1.8682e-04;
% % % % % b__2t=3.4499e-08; 
% % % % % 
% % % % % rt=0.0055;
% % % % % Vt=0.56;

%% Diode 1-2 parameters

U__Bd=300;
b__0d=0;
b__1d=0.000040843026433002205;
b__2d=-3.930898280888438e-8;

rd=0.0024;
Vd=0.97;

%% RMS currents

irms_t14=sqrt(2).*sqrt(3).*(cos(d)+1).*sqrt(m).*abs(I__ac)/(6.*sqrt(pi));

irms_t23=sqrt(3).*sqrt(((-6.*t__b+3).*sin(abs(d))+8.*m.*t__b).*cos(abs(d))+...
    (6.*abs(d)-3.*pi).*t__b-3.*abs(d)+3.*pi).*abs(I__ac)/(6.*sqrt(pi));

% sqrt(3).*sqrt(-2.*m.*(t__a-t__b).*cos(d).^2+((3.*t__a+3.*t__b-3).*sin(d)+...
%     4.*m.*(t__a+t__b)).*cos(d)+(-3.*d-2.*m).*t__a+(-3.*d+2.*m-3.*pi).*t__b+3.*d+3.*pi).*abs(I__ac)/(6.*sqrt(pi));

irms_t56=sqrt(3).*sqrt(-2.*m.*(t__a+t__b).*cos(abs(d))^2+(4.*(m-3.*sin(abs(d)).*(1/4))).*(t__a-t__b).*cos(abs(d))+...
    (-3.*abs(d)-2.*m+3.*pi).*t__b+3.*t__a.*(abs(d)-2.*m.*(1/3))).*abs(I__ac)/(6.*sqrt(pi));


% sqrt(3).*sqrt(-4.*t__b.*(m.*cos(d).^2+m-3.*pi.*(1/4))).*m.*abs(I__ac)/(6.*sqrt(pi));

irms_d14=-sqrt(2).*sqrt(3).*(cos(d)-1).*sqrt(m).*abs(I__ac)/(6.*sqrt(pi));

irms_d23=sqrt(3).*sqrt(((6.*t__b-3).*sin(abs(d))-8.*m.*t__b).*cos(abs(d))+...
    (-6.*abs(d)+3.*pi).*t__b+3.*abs(d)).*abs(I__ac)/(6.*sqrt(pi));

% sqrt(3).*sqrt(2.*m.*(t__a-t__b).*cos(d).^2+((-3.*t__a-3.*t__b+3).*sin(d)-4.*m.*(t__a+...
%     t__b)).*cos(d)+(3.*d+2.*m).*t__a+(3.*d-2.*m+3.*pi).*t__b-3.*d).*abs(I__ac)/(6.*sqrt(pi));

irms_d56=sqrt(3).*sqrt(-2.*cos(abs(d))^2+(4.*(1-3.*sin(abs(d)).*(1/4))).*cos(abs(d))-2+...
    3.*abs(d))/(6.*sqrt(pi));


% sqrt(3).*sqrt(-(4.*(m.*cos(d).^2+m-3.*pi.*(1/4))).*t__a).*abs(I__ac)/(6.*sqrt(pi));

%% AVG currents

iavg_t14=-m.*((-pi+abs(d)).*cos(abs(d))-sin(abs(d))).*I__ac/(4.*pi);

iavg_t23=I__ac*(2+(pi*m*t__b-4*t__b+2)*cos(abs(d)))/(4*pi);

% I__ac.*(((-2+(pi-abs(d)).*m).*t__b+m.*t__a.*abs(d)-2.*t__a+2).*cos(abs(d))-m.*(t__a-...
%     t__b).*sin(abs(d))+2.*t__a-2.*t__b+2)/(4.*pi);

iavg_t56=I__ac*(((2+(-pi+abs(d))*m)*t__b+t__a*(abs(d)*m-2))*cos(abs(d))-(t__a+...
    t__b)*(m*sin(abs(d))-2))/(4*pi);

% I__ac.*(m.*(abs(d)-(1/2).*pi).*cos(abs(d))-m.*sin(abs(d))+2).*t__b/(2.*pi);

iavg_d14=-I__ac.*m.*(cos(abs(d)).*abs(d)-sin(abs(d)))/(4.*pi);

iavg_d23=-(-2+(pi*m*t__b-4*t__b+2)*cos(abs(d)))*I__ac/(4*pi);

% -I__ac.*(((-2+(pi-abs(d)).*m).*t__b+m.*t__a.*abs(d)-2.*t__a+2).*cos(abs(d))-m.*(t__a-...
%     t__b).*sin(abs(d))+2.*t__a-2.*t__b-2)/(4.*pi);

iavg_d56=I__ac*(((2+(-pi+abs(d))*m)*t__a+t__b*(abs(d)*m-2))*cos(abs(d))-...
    (t__a+t__b)*(m*sin(abs(d))-2))/(4*pi);

% I__ac.*(m.*(abs(d)-(1/2).*pi).*cos(abs(d))-m.*sin(abs(d))+2).*t__a/(2.*pi);

%% Conduction losses

Pc_t1=irms_t14.^2*rt+iavg_t14*Vt;

Pc_t2=irms_t23.^2*rt+iavg_t23*Vt;

Pc_t3=Pc_t2;

Pc_t4=Pc_t1;

Pc_t5=irms_t56.^2*rt+iavg_t56*Vt;

Pc_t6=Pc_t5;

Pc_d1=irms_d14.^2*rd+iavg_d14*Vd;

Pc_d2=irms_d23.^2*rd+iavg_d23*Vd;

Pc_d3=Pc_d2;

Pc_d4=Pc_d1;

Pc_d5=irms_d56.^2*rd+iavg_d56*Vd;

Pc_d6=Pc_d5;

Pc_comp=[Pc_t1 Pc_d1 Pc_t2 Pc_d2 Pc_t3 Pc_d3 Pc_t4 Pc_d4 Pc_t5 Pc_d5 Pc_t6 Pc_d6];

Pc_tot=sum(Pc_comp);

%% Switching losses

Ps_t14=-(t__a+(1-t__a)*(1-L__s)).*f__s.*U__DC.*(-(1/2).*I__ac.*(b__2t.*I__ac.*sin(abs(d))+2.*b__1t).*cos(abs(d))+...
    (1/2).*b__2t.*(-pi+abs(d)).*I__ac^2-b__1t.*I__ac+b__0t.*(-pi+abs(d)))/(4.*U__Bt.*pi);

Ps_t23=U__DC.*(-(b__2t.*I__ac.*sin(abs(d))+2.*b__1t).*I__ac.*((-1+t__a).*L__s+t__b).*cos(abs(d))+...
    ((-1+t__a).*(-pi+abs(d)).*L__s+abs(d).*t__b).*b__2t.*I__ac^2+(2.*((1-t__a).*L__s+t__b)).*b__1t.*I__ac+...
    (2.*((-1+t__a).*(-pi+abs(d)).*L__s+abs(d).*t__b)).*b__0t).*f__s/(8.*pi.*U__Bt);

% (-I__ac.*((-1+t__a).*L__s+t__a).*(b__2t.*I__ac.*sin(abs(d))+2.*b__1t).*cos(abs(d))+...
%     b__2t.*((-1+t__a).*(-pi+abs(d)).*L__s+abs(d).*t__a).*I__ac^2-(2.*((-1+...
%     t__a).*L__s-t__a)).*b__1t.*I__ac+(2.*((-1+t__a).*(-pi+abs(d)).*L__s+...
%     abs(d).*t__a)).*b__0t).*U__DC.*f__s/(8.*pi.*U__Bt);

Ps_t56=-(1-t__b).*f__s.*U__DC.*(b__2t.*I__ac^2.*cos(abs(d)).*sin(abs(d))-b__2t.*I__ac^2.*abs(d)+...
    2.*b__1t.*I__ac.*cos(abs(d))-2.*b__1t.*I__ac-2.*b__0t.*abs(d))/(8.*U__Bt.*pi);

Ps_d14=-(1-U__s*t__b).*f__s.*U__DC.*(b__2d.*I__ac^2.*cos(abs(d)).*sin(abs(d))-b__2d.*I__ac^2.*abs(d)+...
    2.*b__1d.*I__ac.*cos(abs(d))-2.*b__1d.*I__ac-2.*b__0d.*abs(d))/(8.*U__Bd.*pi);

Ps_d23=(-I__ac.*(b__2d.*I__ac.*sin(abs(d))+2.*b__1d).*(U__s.*t__b+...
    t__a-1).*cos(abs(d))+b__2d.*((-pi+abs(d)).*t__a+abs(d).*t__b.*U__s-...
    abs(d)+pi).*I__ac^2+2.*b__1d.*(U__s.*t__b-t__a+1).*I__ac+(2.*((-pi+...
    abs(d)).*t__a+abs(d).*t__b.*U__s-abs(d)+pi)).*b__0d).*U__DC.*f__s/(8.*pi.*U__Bd);

% (-(-1+(U__s+1).*t__a).*I__ac.*(b__2d.*I__ac.*sin(abs(d))+2.*b__1d).*cos(abs(d))...
%     +b__2d.*(((U__s+1).*abs(d)-pi).*t__a-abs(d)+pi).*I__ac^2+(2.*(1+...
%     (-1+U__s).*t__a)).*b__1d.*I__ac+(2.*(((U__s+1).*abs(d)-pi).*t__a-abs(d)+...
%     pi)).*b__0d).*U__DC.*f__s/(8.*pi.*U__Bd);

Ps_d56=-t__a.*f__s.*U__DC.*(-(1/2).*I__ac.*(b__2d.*I__ac.*sin(abs(d))+2.*b__1d).*cos(abs(d))+...
    (1/2).*b__2d.*(-pi+abs(d)).*I__ac^2-b__1d.*I__ac+b__0d.*(-pi+abs(d)))/(4.*U__Bd.*pi);


Ps_t1=Ps_t14;

Ps_t2=Ps_t23;

Ps_t3=Ps_t2;

Ps_t4=Ps_t1;

Ps_t5=Ps_t56;

Ps_t6=Ps_t5;

Ps_d1=Ps_d14;

Ps_d2=Ps_d23;

Ps_d3=Ps_d2;

Ps_d4=Ps_d1;

Ps_d5=Ps_d56;

Ps_d6=Ps_d5;

Ps_comp=[Ps_t1 Ps_d1 Ps_t2 Ps_d2 Ps_t3 Ps_d3  Ps_t4 Ps_d4 Ps_t5 Ps_d5 Ps_t6  Ps_d6];

Ps_tot=sum(Ps_comp);

%% total losses per component

Pt1=Pc_t1+Ps_t1;

Pt2=Pc_t2+Ps_t2;

Pt3=Pc_t3+Ps_t3;

Pt4=Pc_t4+Ps_t4;

Pt5=Pc_t5+Ps_t5;

Pt6=Pc_t6+Ps_t6;

Pd1=Pc_d1+Ps_d1;

Pd2=Pc_d2+Ps_d2;

Pd3=Pc_d3+Ps_d3;

Pd4=Pc_d4+Ps_d4;

Pd5=Pc_d5+Ps_d5;

Pd6=Pc_d6+Ps_d6;

Ptot_comp=[Pt1 Pd1 Pt2 Pd2 Pt3 Pd3 Pt4 Pd4 Pt5 Pd5 Pt6 Pd6];

Ptot=(Pc_tot+Ps_tot)*3;

%% Plotting 

Pc_comp=round(Pc_comp);
figure;
bar(transpose([Pc_comp; Ps_comp]),'stacked')
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','T5','D5','T6','D6'})
title('ANPCC phi=0 / ')