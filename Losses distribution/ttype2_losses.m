% figure;
% bar(transpose([Pc_comp; Ps_comp]),'stacked')
% xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4'})
% title('T-type2 phi=0')
% clear all

%%

f__s=20000;
U__DC=900;
Vac=400;
Pac=100000; %% 200 kW
Iac_rms=Pac/(sqrt(3)*Vac);

%%% IGBT 1 parameters

U__Bt1=300;
b__0t1=0; 
b__1t1=9.1730e-05;
b__2t1=-1.5072e-08; 

rt1=0.0056;
Vt1=0.80;

%%% IGBT 2 parameters

U__Bt2=300;
b__0t2=0; 
b__1t2=7.9360e-05;
b__2t2=-3.7275e-08;

rt2=0.0041;
Vt2=0.70;

%%% Diode 1 parameters

U__Bd1=300;
b__0d1=0;
b__1d1=0.00004174076629555432;
b__2d1=-4.5913405069770814e-8;

rd1=0.0014;
Vd1=2.05;

%%% Diode 2 parameters

U__Bd2=300;
b__0d2=0;
b__1d2=0.00004104053136013355;
b__2d2=-4.458356879059264e-8;

rd2=0.0016;
Vd2=1.15;

ddeg=0;
d=ddeg/180.*pi;
I__ac=sqrt(2)*Iac_rms

m=Vac*2*sqrt(2)/sqrt(3)/U__DC%1;%0.05:0.05:0.95;

% d=-pi:1/180.*pi:pi;

i=0;
for f__s=1000:1000:50000
i=i+1;
%%% RMS currents

irms_t23d23=sqrt(3).*sqrt(-4.*m.*cos(d).^2+3.*pi-4.*m).*abs(I__ac)/(6.*sqrt(pi));

irms_t14=sqrt(2).*sqrt(3).*(cos(d)+1).*sqrt(m).*abs(I__ac)/(6.*sqrt(pi));

irms_d14=-sqrt(2).*sqrt(3).*(cos(d)-1).*sqrt(m).*abs(I__ac)/(6.*sqrt(pi));

%%% AVG currents

iavg_t23d23=(m.*(abs(d)-(1/2).*pi).*cos(abs(d))-m.*sin(abs(d))+2).*I__ac/(2.*pi);

iavg_t14=-m.*I__ac.*((-pi+abs(d)).*cos(abs(d))-sin(abs(d)))/(4.*pi);

iavg_d14=m.*I__ac.*(-cos(abs(d)).*abs(d)+sin(abs(d)))/(4.*pi);


%%% Conduction losses

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

%%% Switching losses


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

%%% total losses per component

Pt1=Pc_t1+Ps_t1;

Pt2=Pc_t2+Ps_t2;

Pt3=Pc_t3+Ps_t3;

Pt4=Pc_t4+Ps_t4;

Pd1=Pc_d1+Ps_d1;

Pd2=Pc_d2+Ps_d2;

Pd3=Pc_d3+Ps_d3;

Pd4=Pc_d4+Ps_d4;


Ptot_comp=[Pt1 Pd1 Pt2 Pd2 Pt3 Pd3  Pt4 Pd4 ];


Ptot(i)=(Pc_tot+Ps_tot)*3;
% eff1=(Pac-Ptot)/Pac*100
eff(i)=(Pac-Ptot(i))/Pac*100;

eff_cond(i)=(Pac-Pc_tot)/Pac*100;
eff_sw(i)=(Pac-Ps_tot)/Pac*100;
end

figure(98)
hold on
plot(eff_cond)

figure(99)
hold on
plot(eff_sw,'--g')

figure(100)
hold on
plot(eff)
%% Sweeping delta
% close all
% figure;
% plot(d,irms_t14,d,irms_t23d23,d,irms_d14)
% title('I_{rms}')
% 
% figure;
% plot(d,iavg_t14,d,iavg_t23d23,d,iavg_d14)
% title('I_{avg}')
% 
% 
% figure;
% plot(d,Pc_t1,d,Pc_d1,d,Pc_t2,d,Pc_d2)
% title('Conduction losses')
% 
% figure;
% plot(d,Ps_t14,d,Ps_d14,d,Ps_t23,d,Ps_d23)
% title('Switching losses')

% figure;
% plot(d,Pc_tot,d,Ps_tot,d,Ptot)
% title('Cond and sw losses')

% figure;
% 
% subplot(1,2,2)
% plot(d,irms_t14,d,irms_t23d23,d,irms_d14,'Linewidth',1.5)
% title('I_{rms}')
% xlim([-pi pi])
% xticks([-pi -pi/2 0 pi/2 pi ])
% xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
% xlabel('\delta')
% ylabel('[A]')
% 
% subplot(1,2,1)
% plot(d,iavg_t14,d,iavg_t23d23,d,iavg_d14,'Linewidth',1.5)
% title('I_{avg}')
% xlim([-pi pi])
% xticks([-pi -pi/2 0 pi/2 pi ])
% xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
% legend('T1-T4','T2-T3-D2-D3','D1-D4')
% ylabel('[A]')
% xlabel('\delta')
% suptitle('T-type 2 currents')
%% Sweeping m
% close all
% figure;
% plot(m,irms_t23d56,m,irms_t14,m,irms_d1234)
% title('I_{rms}')
% 
% figure;
% plot(m,iavg_t23d56,m,iavg_t14,m,iavg_d1234)
% title('I_{avg}')
% 
% figure;
% hold on
% plot(m,Pc_tot,m,Ps_tot,m,Ptot)
% scatter(m,Ptot)
% title('Cond and sw losses')

figure;
bar(transpose([Pc_comp; Ps_comp]),'stacked')
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','D5','D6'})