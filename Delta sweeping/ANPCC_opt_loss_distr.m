% figure;
% bar(transpose([Pc_comp; Ps_comp]),'stacked')
% xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','T5','D5','T6','D6'})
% title('ANPCC phi=0 / ta=1')
clear all
%%
f__s=20000;
U__DC=900;
Vac=400;
Pac=100000; %% 200 kW
Iac_rms=Pac/(sqrt(3)*Vac)

t__a=0.5; %%% follows path P<-->0U1
t__b=1-t__a;

ldf=t__a; %%% it means i am using the sequence of thiago paper to optimally distribute losses
%%% ldf = 1 i only use sequence P<-->0U1 of thiago paper

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

U__Bd=600;
b__0d=0;
b__1d=0.00010777181728433188;
b__2d=-9.544835941123982e-8;

rd=0.0046;
Vd=0.98;

ddeg=0;
d=ddeg/180.*pi;
I__ac=sqrt(2)*Iac_rms %%% peak AC current

m=Vac*2*sqrt(2)/sqrt(3)/U__DC%1;%0.05:0.05:0.95;

% d=-pi:1/180.*pi:pi;
i=0;
for f__s=1000:1000:50000
i=i+1;
%%% RMS currents

irms_t14=sqrt(2).*sqrt(3).*(cos(d)+1).*sqrt(m).*abs(I__ac)/(6.*sqrt(pi));

irms_t23=sqrt(3).*sqrt(-2.*m.*(t__a-t__b).*cos(d).^2+((3.*t__a+3.*t__b-3).*sin(d)+...
    4.*m.*(t__a+t__b)).*cos(d)+(-3.*d-2.*m).*t__a+(-3.*d+2.*m-3.*pi).*t__b+3.*d+3.*pi).*abs(I__ac)/(6.*sqrt(pi));

irms_t56=sqrt(3).*sqrt(-4.*t__b.*(m.*cos(d).^2+m-3.*pi.*(1/4))).*m.*abs(I__ac)/(6.*sqrt(pi));

irms_d14=-sqrt(2).*sqrt(3).*(cos(d)-1).*sqrt(m).*abs(I__ac)/(6.*sqrt(pi));

irms_d23=sqrt(3).*sqrt(2.*m.*(t__a-t__b).*cos(d).^2+((-3.*t__a-3.*t__b+3).*sin(d)-4.*m.*(t__a+...
    t__b)).*cos(d)+(3.*d+2.*m).*t__a+(3.*d-2.*m+3.*pi).*t__b-3.*d).*abs(I__ac)/(6.*sqrt(pi));

irms_d56=sqrt(3).*sqrt(-(4.*(m.*cos(d).^2+m-3.*pi.*(1/4))).*t__a).*abs(I__ac)/(6.*sqrt(pi));

%%% AVG currents

iavg_t14=-m.*((-pi+abs(d)).*cos(abs(d))-sin(abs(d))).*I__ac/(4.*pi);

iavg_t23=I__ac.*(((-2+(pi-abs(d)).*m).*t__b+m.*t__a.*abs(d)-2.*t__a+2).*cos(abs(d))-m.*(t__a-...
    t__b).*sin(abs(d))+2.*t__a-2.*t__b+2)/(4.*pi);

iavg_t56=I__ac.*(m.*(abs(d)-(1/2).*pi).*cos(abs(d))-m.*sin(abs(d))+2).*t__b/(2.*pi);

iavg_d14=-I__ac.*m.*(cos(abs(d)).*abs(d)-sin(abs(d)))/(4.*pi);

iavg_d23=-I__ac.*(((-2+(pi-abs(d)).*m).*t__b+m.*t__a.*abs(d)-2.*t__a+2).*cos(abs(d))-m.*(t__a-...
    t__b).*sin(abs(d))+2.*t__a-2.*t__b-2)/(4.*pi);

iavg_d56=I__ac.*(m.*(abs(d)-(1/2).*pi).*cos(abs(d))-m.*sin(abs(d))+2).*t__a/(2.*pi);

%%% Conduction losses

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

%%% Switching losses

Ps_t14=-f__s.*U__DC.*ldf.*(-(1/2).*I__ac.*(b__2t.*I__ac.*sin(abs(d))+2.*b__1t).*cos(abs(d))+...
    (1/2).*b__2t.*(-pi+abs(d)).*I__ac^2-b__1t.*I__ac+b__0t.*(-pi+abs(d)))/(4.*U__Bt.*pi);

Ps_t23=-U__DC.*(((1/2).*b__2t.*I__ac^2+b__0t).*pi+2.*b__1t.*I__ac).*f__s.*(-1+ldf)/(4.*U__Bt.*pi);

Ps_t56=-ldf.*f__s.*U__DC.*(b__2t.*I__ac^2.*cos(abs(d)).*sin(abs(d))-b__2t.*I__ac^2.*abs(d)+...
    2.*b__1t.*I__ac.*cos(abs(d))-2.*b__1t.*I__ac-2.*b__0t.*abs(d))/(8.*U__Bt.*pi);

Ps_d14=-ldf.*f__s.*U__DC.*(b__2d.*I__ac^2.*cos(abs(d)).*sin(abs(d))-b__2d.*I__ac^2.*abs(d)+...
    2.*b__1d.*I__ac.*cos(abs(d))-2.*b__1d.*I__ac-2.*b__0d.*abs(d))/(8.*U__Bd.*pi);

Ps_d23=-U__DC.*(((1/2).*b__2d.*I__ac^2+b__0d).*pi+2.*b__1d.*I__ac).*f__s.*(-1+ldf)/(4.*U__Bd.*pi);

Ps_d56=-f__s.*U__DC.*ldf.*(-(1/2).*I__ac.*(b__2d.*I__ac.*sin(abs(d))+2.*b__1d).*cos(abs(d))+...
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

%%% total losses per component

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
plot(eff_sw)%,'--r')

figure(100)
hold on
plot(eff)
%% Sweeping delta


% close all
% figure;
% plot(d,irms_t14,d,irms_t23,d,irms_t56,d,irms_d14,d,irms_d23,d,irms_d56)
% 
% title('I_{rms}')
% 
% figure;
% plot(d,iavg_t14,d,iavg_t23,d,iavg_t56,d,iavg_d14,d,iavg_d23,d,iavg_d56)
% 
% title('I_{avg}')
% 
% figure;
% plot(d,Pc_t1,d,Pc_t2,d,Pc_d1,d,Pc_d2,d,Pc_d5)
% legend
% title('Conduction losses')
% 
% figure;
% hold on
% plot(d,Ps_t14,d,Ps_t56,d,Ps_d14,d,Ps_d56)
% plot(d,Ps_t23*ones(1,361))
% hold off
% legend
% title('Switching losses')
% 
% figure;
% plot(d,Pc_tot,d,Ps_tot,d,Ptot)
% title('Cond and sw losses')
% close all
% figure;
% 
% subplot(1,2,2)
% plot(d,irms_t14,d,irms_t23,d,irms_t56,d,irms_d14,d,irms_d23,d,irms_d56,'Linewidth',1.5)
% title('I_{rms}')
% xlim([-pi pi])
% xticks([-pi -pi/2 0 pi/2 pi ])
% xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
% xlabel('\delta')
% ylabel('[A]')
% 
% subplot(1,2,1)
% plot(d,iavg_t14,d,iavg_t23,d,iavg_t56,d,iavg_d14,d,iavg_d23,d,iavg_d56,'Linewidth',1.5)
% title('I_{avg}')
% xlim([-pi pi])
% xticks([-pi -pi/2 0 pi/2 pi ])
% xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
% legend('T1-T4','T2-T3','T5-T6','D1-D4','D2-D3','D5-D6')
% ylabel('[A]')
% xlabel('\delta')
% % suptitle('ANPCC currents - m_{a} = 0.5')

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
% xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','T5','D5','T6','D6'})
% title('ANPCC phi=0 / ')