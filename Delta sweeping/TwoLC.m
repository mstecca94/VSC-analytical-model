% figure;
% bar(transpose([Pc_comp; Ps_comp]),'stacked')
% xticklabels({'T1','D1','T2','D2'})
% title('TwoLC phi=0')

%%


f__s=8000;
U__DC=900;
Vac=400;
Pac=100000; %% 200 kW
Iac_rms=Pac/(sqrt(3)*Vac)
%%% IGBT parameters

U__Bt=600;
b__0t=0; 
b__1t=1.7328e-04;
b__2t=4.8252e-08; 

rt=0.0050;
Vt=0.56;

%%% Diode parameters

U__Bd=600;
b__0d=0;
b__1d=0.00011336152231240664;
b__2d=-1.1406623528819162e-7;


rd=0.0038;
Vd=0.50;

ddeg=0;
d=ddeg/180.*pi;
I__ac=sqrt(2)*Iac_rms %%% peak AC current

m=1;%0.05:0.05:0.95;

% d=-pi:1/180.*pi:pi;

i=0;
for f__s=1000:1000:50000
i=i+1;
%%% RMS currents

irms_t12=sqrt(3)*sqrt((8*m-3*sin(abs(d))).*cos(abs(d))+3*abs(d))*abs(I__ac)/(6*sqrt(pi));

irms_d12=sqrt(3)*sqrt((-8*m+3*sin(abs(d))).*cos(abs(d))-3*abs(d)+3*pi)*abs(I__ac)/(6*sqrt(pi));


%%% AVG currents

iavg_t12=I__ac.*(m.*cos(d).*pi-2.*cos(d)+2)/(4.*pi);

iavg_d12=-I__ac.*(m.*cos(d).*pi-2.*cos(d)-2)/(4.*pi);


%%% Conduction losses

Pc_t1=irms_t12.^2*rt+iavg_t12*Vt;

Pc_t2=Pc_t1;

Pc_d1=irms_d12.^2*rt+iavg_d12*Vt;

Pc_d2=Pc_d1;

Pc_comp=[Pc_t1 Pc_d1 Pc_t2 Pc_d2];

Pc_tot=sum(Pc_comp);

%%% Switching losses

Ps_t12=2*f__s*U__DC*(((1/2)*b__2t*I__ac^2+b__0t)*pi+2*b__1t*I__ac)/(4*pi*U__Bt);

Ps_d12=2*f__s*U__DC*(((1/2)*b__2d*I__ac^2+b__0d)*pi+2*b__1d*I__ac)/(4*pi*U__Bd);

Ps_t1=Ps_t12;

Ps_t2=Ps_t1;

Ps_d1=Ps_d12;

Ps_d2=Ps_d1;

Ps_comp=[Ps_t1 Ps_d1 Ps_t2 Ps_d2];

Ps_tot=sum(Ps_comp);

%%% total losses per component

Pt1=Pc_t1+Ps_t1;

Pt2=Pc_t2+Ps_t2;

Pd1=Pc_d1+Ps_d1;

Pd2=Pc_d2+Ps_d2;

Ptot_comp=[Pt1 Pd1 Pt2 Pd2];

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
plot(eff_sw,'--c')

figure(100)
hold on
plot(eff)
%% Sweeping delta
% close all
% figure;
% plot(d,irms_t12,d,irms_d12)
% title('I_{rms}')
% 
% figure;
% plot(d,iavg_t12,d,iavg_d12)
% title('I_{avg}')
% 
% figure;
% plot(d,Pc_t1,d,Pc_d1)
% title('Conduction losses')
% 
% figure;
% plot(d,Ps_t1,d,Ps_d1)
% title('Switching losses')
% 
% figure;
% plot(d,Pc_tot,d,Ps_tot,d,Ptot)
% title('Cond and sw losses')
% figure;
% 
% subplot(1,2,2)
% plot(d,irms_t12,d,irms_d12,'Linewidth',1.5)
% title('I_{rms}')
% xlim([-pi pi])
% xticks([-pi -pi/2 0 pi/2 pi ])
% xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
% xlabel('\delta')
% ylabel('[A]')
% 
% subplot(1,2,1)
% plot(d,iavg_t12,d,iavg_d12,'Linewidth',1.5)
% title('I_{avg}')
% xlim([-pi pi])
% xticks([-pi -pi/2 0 pi/2 pi ])
% xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
% legend('T1-T2','D1-D2')
% ylabel('[A]')
% xlabel('\delta')
% suptitle('2LC currents')
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
% figure;
% bar(transpose([Pc_comp; Ps_comp]),'stacked')
% xticklabels({'T1','D1','T2','D2'})