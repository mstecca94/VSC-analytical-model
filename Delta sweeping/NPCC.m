% figure;
% bar(transpose([Pc_comp; Ps_comp]),'stacked')
% xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','D5','D6'})
% title('NPCC phi=180')
% clear all

%%

f__s=20000;
U__DC=900;
Vac=400;
Pac=100000; %% 100 kW
Iac_rms=Pac/(sqrt(3)*Vac);

%%% IGBT 1 parameters

U__Bt=600;
b__0t=0; 
b__1t=1.7410e-04 ;
b__2t=5.4937e-08; 

rt=0.0055;
Vt=0.56;

% % % % % %%% IGBT 2 parameters
% % % % % 
% % % % % U__Bt=600;
% % % % % b__0t=0; 
% % % % % b__1t=1.8682e-04;
% % % % % b__2t=3.4499e-08; 
% % % % % 
% % % % % rt=0.0055;
% % % % % Vt=0.56;

%%%% IGBT 1 and 2 Eon and Eoff are same, slightly different due to the
%%%% curve fitting - I can assume them equal

%%% Diode 1-2 parameters

U__Bd12=600;
b__0d12=0;
b__1d12=0.00010777181728433188;
b__2d12=-9.544835941123982e-8;

rd12=0.0046;
Vd12=0.98;

%%% Diode 5 parameters

U__Bd5=600;
b__0d5=0;
b__1d5=0.00010895583901692339;
b__2d5=-8.525717376201372e-8;

rd5=0.0044;
Vd5=0.95;

ddeg=0;
d=ddeg/180.*pi;
I__ac=sqrt(2)*Iac_rms

m=Vac*2*sqrt(2)/sqrt(3)/U__DC%1;%0.05:0.05:0.95;

% d=-pi:1/180.*pi:pi;
i=0;
for f__s=1000:1000:50000
i=i+1;
%%% RMS currents

irms_t23=sqrt(3).*sqrt(-2.*m.*cos(d).^2+4.*m.*cos(d)+3.*pi-2.*m).*abs(I__ac)/(6.*sqrt(pi));

irms_t14=sqrt(2).*sqrt(3).*(cos(d)+1).*sqrt(m).*abs(I__ac)/(6.*sqrt(pi));

irms_d1234=-sqrt(2).*sqrt(3).*(cos(d)-1).*sqrt(m).*abs(I__ac)/(6.*sqrt(pi));

irms_d56=sqrt(3).*sqrt(-(2.*(m.*cos(2.*d)+3.*m-3.*pi.*(1/2))).*I__ac^2)/(6.*sqrt(pi));

%%% AVG currents

iavg_t23=-I__ac.*(-m.*cos(abs(d)).*abs(d)+m.*sin(abs(d))-4)/(4.*pi);

iavg_t14=-I__ac.*((-pi+abs(d)).*cos(abs(d))-sin(abs(d))).*m/(4.*pi);

iavg_d1234=m.*I__ac.*(-cos(abs(d)).*abs(d)+sin(abs(d)))/(4.*pi);

iavg_d56=I__ac.*(m.*(abs(d)-(1/2).*pi).*cos(abs(d))-m.*sin(abs(d))+2)/(2.*pi);

%%% Conduction losses

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

%%% Switching losses

Ps_t14=-(-(1/2).*I__ac.*(b__2t.*I__ac.*sin(abs(d))+2.*b__1t).*cos(abs(d))+...
    (1/2).*b__2t.*(-pi+abs(d)).*I__ac^2-b__1t.*I__ac+b__0t.*(-pi+abs(d))).*U__DC.*f__s/(4.*pi.*U__Bt);

Ps_t23=-f__s.*U__DC.*(b__2t.*I__ac^2.*cos(abs(d)).*sin(abs(d))-b__2t.*I__ac^2.*abs(d)+...
    2.*b__1t.*I__ac.*cos(abs(d))-2.*b__1t.*I__ac-2.*b__0t.*abs(d))/(8.*pi.*U__Bt);

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

%%% total losses per component

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


Ptot(i)=(Pc_tot+Ps_tot)*3;
% Ptot=(Pc_tot+Ps_tot)*3;
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
plot(eff_sw,'ok')

figure(100)
hold on
plot(eff)


%% Sweeping delta
% close all
% figure
% 
% subplot(1,2,2)
% plot(d,irms_t14,d,irms_t23,d,irms_d1234,d,irms_d56,'Linewidth',1.5)
% title('I_{rms}')
% xlim([-pi pi])
% xticks([-pi -pi/2 0 pi/2 pi ])
% xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
% xlabel('\delta')
% ylabel('[A]')
% 
% subplot(1,2,1)
% plot(d,iavg_t14,d,iavg_t23,d,iavg_d1234,d,iavg_d56,'Linewidth',1.5)
% title('I_{avg}')
% xlim([-pi pi])
% xticks([-pi -pi/2 0 pi/2 pi ])
% xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
% legend('T1-T4','T2-T3','D1-D2-D3-D4','D5-D6')
% ylabel('[A]')
% xlabel('\delta')
% suptitle('NPCC currents')

% figure;
% plot(d,Pc_t1,d,Pc_t2,d,Pc_d1,d,Pc_d2,d,Pc_d5)
% title('Conduction losses')
% 
% figure;
% plot(d,Ps_t14,d,Ps_t23,d,Ps_d14,d,Ps_d56)
% title('Switching losses')
% 
% figure;
% plot(d,Pc_tot,d,Ps_tot,d,Ptot)
% title('Cond and sw losses')

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
% 
figure;
bar(transpose([Pc_comp; Ps_comp]),'stacked')
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','D5','D6'})