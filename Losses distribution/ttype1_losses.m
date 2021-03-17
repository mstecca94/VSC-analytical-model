% figure;
% bar(transpose([Pc_comp; Ps_comp]),'stacked')
% xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','D5','D6'})
% title('T-type1 phi=180')
clear all

%%
f__s=20000;
U__DC=1000;
Vac=690;
Pac=250000; %% 200 kW
Iac_rms=Pac/(sqrt(3)*Vac);
%%% IGBT parameters

U__Bt=600;
b__0t=0; 
b__1t=1.8738e-04 ;
b__2t=2.4378e-08; 

rt=0.0039;
Vt=0.61;

%%% Diode parameters

U__Bd=600;
b__0d=0;
b__1d=0.00011521308191658408;
b__2d=-1.1367801754859876e-7;

rd=0.0024;
Vd=0.63;

ddeg=0;
d=ddeg/180.*pi;
I__ac=sqrt(2)*Iac_rms

m=1;%0.05:0.05:0.95;

d=-pi:1/180.*pi:pi;
% i=0;
% for f__s=1000:1000:50000
% i=i+1;
%%% RMS currents

irms_t23d56=sqrt(3).*sqrt(-4.*m.*cos(d).^2+3.*pi-4.*m).*abs(I__ac)/(6.*sqrt(pi));

irms_t14=sqrt(2).*sqrt(3).*(cos(d)+1).*sqrt(m).*abs(I__ac)/(6.*sqrt(pi));

irms_d1234=-sqrt(2).*sqrt(3).*(cos(d)-1).*sqrt(m).*abs(I__ac)/(6.*sqrt(pi));

%%% AVG currents

iavg_t23d56=I__ac.*(m.*(abs(d)-(1/2).*pi).*cos(abs(d))-sin(abs(d)).*m+2)/(2.*pi);

iavg_t14=-((-pi+abs(d)).*cos(abs(d))-sin(abs(d))).*m.*I__ac/(4.*pi);

iavg_d1234=I__ac.*m.*(-cos(abs(d)).*abs(d)+sin(abs(d)))/(4.*pi);


%%% Conduction losses

Pc_t1=irms_t14.^2*rt+iavg_t14*Vt;

Pc_t2=irms_t23d56.^2*rt+iavg_t23d56*Vt;

Pc_t3=Pc_t2;

Pc_t4=Pc_t1;

Pc_d1=irms_d1234.^2*rd+iavg_d1234*Vd;

Pc_d2=irms_d1234.^2*rd+iavg_d1234*Vd;

Pc_d3=Pc_d2;

Pc_d4=Pc_d1;

Pc_d5=irms_t23d56.^2*rd+iavg_t23d56*Vd;

Pc_d6=Pc_d5;

Pc_comp=[Pc_t1 Pc_d1 Pc_t2 Pc_d2 Pc_t3 Pc_d3  Pc_t4 Pc_d4 Pc_d5  Pc_d6];

Pc_tot=sum(Pc_comp);

%%% Switching losses

Ps_t14=-(-(1/2).*I__ac.*(b__2t.*I__ac.*sin(abs(d))+2.*b__1t).*cos(abs(d))+...
    (1/2).*b__2t.*(-pi+abs(d)).*I__ac^2-b__1t.*I__ac+b__0t.*(-pi+abs(d))).*U__DC.*f__s/(4.*pi.*U__Bt);

Ps_d56=-(-(1/2).*I__ac.*(b__2d.*I__ac.*sin(abs(d))+2.*b__1d).*cos(abs(d))+...
    (1/2).*b__2d.*(-pi+abs(d)).*I__ac^2-b__1d.*I__ac+b__0d.*(-pi+abs(d))).*U__DC.*f__s/(4.*pi.*U__Bd);

Ps_t23=-f__s.*U__DC.*(b__2t.*I__ac^2.*cos(abs(d)).*sin(abs(d))-b__2t.*I__ac^2.*abs(d)+...
    2.*b__1t.*I__ac.*cos(abs(d))-2.*b__1t.*I__ac-2.*b__0t.*abs(d))/(8.*pi.*U__Bt);

Ps_d14=-f__s.*U__DC.*(b__2d.*I__ac^2.*cos(abs(d)).*sin(abs(d))-b__2d.*I__ac^2.*abs(d)+...
    2.*b__1d.*I__ac.*cos(abs(d))-2.*b__1d.*I__ac-2.*b__0d.*abs(d))/(8.*pi.*U__Bd);

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

% Ptot(i)=(Pc_tot+Ps_tot)*3;
% % eff1=(Pac-Ptot)/Pac*100
% eff(i)=(Pac-Ptot(i))/Pac*100;
% 
% eff_cond(i)=(Pac-Pc_tot)/Pac*100;
% eff_sw(i)=(Pac-Ps_tot)/Pac*100;
% end
% 
% figure(98)
% hold on
% plot(eff_cond)
% 
% figure(99)
% hold on
% plot(eff_sw,'ob')
% 
% figure(100)
% hold on
% plot(eff)
%% Sweeping delta
% close all
% figure;
% plot(d,irms_t23d56,d,irms_t14,d,irms_d1234)
% title('I_{rms}')
% 
% figure;
% plot(d,iavg_t23d56,d,iavg_t14,d,iavg_d1234)
% title('I_{avg}')
% 
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

figure;

subplot(1,2,2)
plot(d,irms_t14,d,irms_t23d56,d,irms_d1234,'Linewidth',1.5)
title('I_{rms}')
xlim([-pi pi])
xticks([-pi -pi/2 0 pi/2 pi ])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
xlabel('\delta')
ylabel('[A]')

subplot(1,2,1)
plot(d,iavg_t14,d,iavg_t23d56,d,iavg_d1234,'Linewidth',1.5)
title('I_{avg}')
xlim([-pi pi])
xticks([-pi -pi/2 0 pi/2 pi ])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
legend('T1-T4','T2-T3-D5-D6','D1-D2-D3-D4')
ylabel('[A]')
xlabel('\delta')
suptitle('T-type 1 currents')

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
% hold on% figure;
% subplot(1,4,1)
% bar(transpose([Pc_comp; Ps_comp]),'stacked')
% xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','D5','D6'})
% 
% subplot(1,4,2)
% bar(transpose([Pc_comp; Ps_comp]),'stacked')
% xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','D5','D6'})
% 
% subplot(1,4,3)
% bar(transpose([Pc_comp; Ps_comp]),'stacked')
% xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','D5','D6'})
% 
% subplot(1,4,4)
% bar(transpose([Pc_comp; Ps_comp]),'stacked')
% xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','D5','D6'})

% plot(m,Pc_tot,m,Ps_tot,m,Ptot)
% scatter(m,Ptot)
% title('Cond and sw losses')

