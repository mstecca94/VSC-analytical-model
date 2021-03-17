ccc
%% Design parameters
U__DC=900;
Vac=400;
Pac=100000; %% 100 kW
ddeg=0;
figindx=1;
%% Sweeping frequency - Inverter operation delta=0
f__s1=4000;
f__s2=16000;
i=0;

for f__s=f__s1:1000:f__s2
i=i+1;
[Pc_comp_ANPCC,Ps_comp_ANPCC,Ptot_comp_ANPCC,Ptot_ANPCC]=anpcc_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_ANPCC(i)=(Pac-Ptot_ANPCC)/Pac*100;

[Pc_comp_NPCC,Ps_comp_NPCC,Ptot_comp_NPCC,Ptot_NPCC]=npcc_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_NPCC(i)=(Pac-Ptot_NPCC)/Pac*100;

[Pc_comp_TTYPE,Ps_comp_TTYPE,Ptot_comp_TTYPE,Ptot_TTYPE]=ttype2_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_TTYPE(i)=(Pac-Ptot_TTYPE)/Pac*100;

[Pc_comp_2LC,Ps_comp_2LC,Ptot_comp_2LC,Ptot_2LC]=TwoLC_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_2LC(i)=(Pac-Ptot_2LC)/Pac*100;

end

eff_ANPCC_inv=eff_ANPCC;
eff_NPCC_inv=eff_NPCC;
eff_TTYPE_inv=eff_TTYPE;
eff_2LC_inv=eff_2LC;

f__s_range=f__s1/1000:1:f__s2/1000;
figure(figindx)
plot(f__s_range,eff_ANPCC,f__s_range,eff_NPCC,f__s_range,eff_TTYPE,f__s_range,eff_2LC,'LineWidth',2)
legend('ANPCC','NPCC','T-Type','2 LC')
xlabel('Switching Frequency [kHz]')
ylabel('Efficiency [%]')
xlim([f__s1/1000 f__s2/1000])
title('Inverter operation')
% title('Semiconductor efficiency varying f_{s} - \delta=0^{circ}')

figindx=figindx+1;

%% Sweeping frequency - Inverter operation delta=0

f__s1=4000;
f__s2=16000;
i=0;
ddeg=180;
for f__s=f__s1:1000:f__s2
i=i+1;
[Pc_comp_ANPCC,Ps_comp_ANPCC,Ptot_comp_ANPCC,Ptot_ANPCC]=anpcc_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_ANPCC(i)=(Pac-Ptot_ANPCC)/Pac*100;

[Pc_comp_NPCC,Ps_comp_NPCC,Ptot_comp_NPCC,Ptot_NPCC]=npcc_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_NPCC(i)=(Pac-Ptot_NPCC)/Pac*100;

[Pc_comp_TTYPE,Ps_comp_TTYPE,Ptot_comp_TTYPE,Ptot_TTYPE]=ttype2_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_TTYPE(i)=(Pac-Ptot_TTYPE)/Pac*100;

[Pc_comp_2LC,Ps_comp_2LC,Ptot_comp_2LC,Ptot_2LC]=TwoLC_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_2LC(i)=(Pac-Ptot_2LC)/Pac*100;

end

f__s_range=f__s1/1000:1:f__s2/1000;
figure(figindx)
plot(f__s_range,eff_ANPCC_inv,f__s_range,eff_NPCC_inv,f__s_range,eff_TTYPE_inv,f__s_range,eff_2LC_inv,'LineWidth',2)
xlabel('Switching Frequency [kHz]')
ylabel('Efficiency [%]')
% title('Semiconductor efficiency varying f_{s} - \delta=180^{circ}')
title('Rectifier operation')
xlim([f__s1/1000 f__s2/1000])
figindx=figindx+1;

%%
figure(figindx)
subplot(1,2,1)
hold on
plot(f__s_range,eff_ANPCC_inv,'k','LineWidth',2)
plot(f__s_range,eff_NPCC_inv,'g','LineWidth',2)
plot(f__s_range,eff_TTYPE_inv,'b','LineWidth',2)
plot(f__s_range,eff_2LC_inv,'r','LineWidth',2)
grid on
xlabel('Switching Frequency [kHz]')
ylabel('Efficiency [%]')
% title('Semiconductor efficiency varying f_{s} - \delta=180^{circ}')
title('Rectifier operation')
xlim([f__s1/1000 f__s2/1000])
title('Inverter operation')

subplot(1,2,2)
hold on
plot(f__s_range,eff_ANPCC,'k','LineWidth',2)
plot(f__s_range,eff_NPCC,'g','LineWidth',2)
plot(f__s_range,eff_TTYPE,'b','LineWidth',2)
plot(f__s_range,eff_2LC,'r','LineWidth',2)
grid on
legend('ANPCC','NPCC','T-Type','2 LC')
xlabel('Switching Frequency [kHz]')
ylabel('Efficiency [%]')
title('Rectifier operation')
xlim([f__s1/1000 f__s2/1000])
figindx=figindx+1;

figure(figindx)
subplot(1,2,1)
hold on
plot(f__s_range,eff_ANPCC_inv,'k','LineWidth',2)
plot(f__s_range,eff_NPCC_inv,'r','LineWidth',2)
plot(f__s_range,eff_TTYPE_inv,'b','LineWidth',2)
plot(f__s_range,eff_2LC_inv,'g','LineWidth',2)
grid on
xlabel('Switching Frequency [kHz]')
ylabel('Efficiency [%]')
% title('Semiconductor efficiency varying f_{s} - \delta=180^{circ}')
title('Rectifier operation')
xlim([f__s1/1000 f__s2/1000])
title('Inverter operation')

subplot(1,2,2)
hold on
plot(f__s_range,eff_ANPCC,'k','LineWidth',2)
plot(f__s_range,eff_NPCC,'r','LineWidth',2)
plot(f__s_range,eff_TTYPE,'b','LineWidth',2)
plot(f__s_range,eff_2LC,'g','LineWidth',2)
grid on
legend('ANPCC','NPCC','T-Type','2 LC')
xlabel('Switching Frequency [kHz]')
ylabel('Efficiency [%]')
title('Rectifier operation')
xlim([f__s1/1000 f__s2/1000])
figindx=figindx+1;

figure(figindx)
subplot(1,2,1)
hold on
plot(f__s_range,eff_NPCC_inv,'r','LineWidth',2)
plot(f__s_range,eff_TTYPE_inv,'b','LineWidth',2)
plot(f__s_range,eff_2LC_inv,'g','LineWidth',2)
grid on
xlabel('Switching Frequency [kHz]')
ylabel('Efficiency [%]')
% title('Semiconductor efficiency varying f_{s} - \delta=180^{circ}')
title('Rectifier operation')
xlim([f__s1/1000 f__s2/1000])
title('Inverter operation')
legend('NPCC','T-Type','2 LC')

subplot(1,2,2)
hold on
plot(f__s_range,eff_NPCC,'r','LineWidth',2)
plot(f__s_range,eff_TTYPE,'b','LineWidth',2)
plot(f__s_range,eff_2LC,'g','LineWidth',2)
grid on
legend('NPCC','T-Type','2 LC')
xlabel('Switching Frequency [kHz]')
ylabel('Efficiency [%]')
title('Rectifier operation')
xlim([f__s1/1000 f__s2/1000])
figindx=figindx+1;
%%

figure(figindx)
subplot(2,1,1)
hold on
plot(f__s_range,eff_ANPCC_inv,'r','LineWidth',2)
plot(f__s_range,eff_NPCC_inv,'b','LineWidth',2)
plot(f__s_range,eff_TTYPE_inv,'g','LineWidth',2)
plot(f__s_range,eff_2LC_inv,'k','LineWidth',2)
grid on
xlabel('Switching Frequency [kHz]')
ylabel('Efficiency [%]')
% title('Semiconductor efficiency varying f_{s} - \delta=180^{circ}')
title('Rectifier operation')
xlim([f__s1/1000 f__s2/1000])
title('Inverter operation')

subplot(2,1,2)
hold on
plot(f__s_range,eff_ANPCC,'r','LineWidth',2)
plot(f__s_range,eff_NPCC,'b','LineWidth',2)
plot(f__s_range,eff_TTYPE,'g','LineWidth',2)
plot(f__s_range,eff_2LC,'k','LineWidth',2)
grid on
legend('ANPCC','NPCC','T-Type','2 LC')
xlabel('Switching Frequency [kHz]')
ylabel('Efficiency [%]')
% title('Semiconductor efficiency varying f_{s} - \delta=180^{circ}')
title('Rectifier operation')
xlim([f__s1/1000 f__s2/1000])
figindx=figindx+1;

%% Sweeping delta - Inverter operation delta=0

f__s=12000; % f at 10 kHz
delta_1=-180;
delta_2=180;
i=0;

for d=delta_1:1:delta_2
i=i+1;
ddeg=d;
[Pc_comp_ANPCC,Ps_comp_ANPCC,Ptot_comp_ANPCC,Ptot_ANPCC]=anpcc_freq(U__DC,Vac,Pac,ddeg,f__s);

Pc_comp_ANPCC_delta(i,:)=Pc_comp_ANPCC;
Ps_comp_ANPCC_delta(i,:)=Ps_comp_ANPCC;
Ptot_comp_ANPCC_delta(i,:)=Ptot_comp_ANPCC;

eff_ANPCC_delta(i)=(Pac-Ptot_ANPCC)/Pac*100;

[Pc_comp_NPCC,Ps_comp_NPCC,Ptot_comp_NPCC,Ptot_NPCC]=npcc_freq(U__DC,Vac,Pac,ddeg,f__s);

Pc_comp_NPCC_delta(i,:)=Pc_comp_NPCC;
Ps_comp_NPCC_delta(i,:)=Ps_comp_NPCC;
Ptot_comp_NPCC_delta(i,:)=Ptot_comp_NPCC;

eff_NPCC_delta(i)=(Pac-Ptot_NPCC)/Pac*100;

[Pc_comp_TTYPE,Ps_comp_TTYPE,Ptot_comp_TTYPE,Ptot_TTYPE]=ttype2_freq(U__DC,Vac,Pac,ddeg,f__s);

Pc_comp_TTYPE_delta(i,:)=Pc_comp_TTYPE;
Ps_comp_TTYPE_delta(i,:)=Ps_comp_TTYPE;
Ptot_comp_TTYPE_delta(i,:)=Ptot_comp_TTYPE;

eff_TTYPE_delta(i)=(Pac-Ptot_TTYPE)/Pac*100;

[Pc_comp_2LC,Ps_comp_2LC,Ptot_comp_2LC,Ptot_2LC]=TwoLC_freq(U__DC,Vac,Pac,ddeg,f__s);

Pc_comp_2LC_delta(i,:)=Pc_comp_2LC;
Ps_comp_2LC_delta(i,:)=Ps_comp_2LC;
Ptot_comp_2LC_delta(i,:)=Ptot_comp_2LC;

eff_2LC_delta(i)=(Pac-Ptot_2LC)/Pac*100;

end

delta_range=delta_1:1:delta_2;
figure(figindx)
plot(delta_range,eff_ANPCC_delta,delta_range,eff_NPCC_delta...
    ,delta_range,eff_TTYPE_delta,delta_range,eff_2LC_delta,'LineWidth',2)
title('Semiconductor efficiency')
xlabel('\delta [\circ]')
xlim([-180 180])
ylabel('Efficiency [%]')
legend('ANPCC','NPCC','T-Type','2 LC')
figindx=figindx+1;


figure(figindx)
sgtitle('Conduction losses varying \delta')
subplot(2,2,1)
hold on
plot(delta_range,Pc_comp_ANPCC_delta(:,1:4),'LineWidth',2)
plot(delta_range,Pc_comp_ANPCC_delta(:,9:10),'LineWidth',2)
hold off
legend('T1-T4','D1-D4','T2-T3','D2-D3','T5-T6','D5-D6')
title('ANPCC')

subplot(2,2,2)
hold on
plot(delta_range,Pc_comp_NPCC_delta(:,1:4),'LineWidth',2)
plot(delta_range,Pc_comp_NPCC_delta(:,9))
hold off
legend('T1-T4','D1-D4','T2-T3','D2-D3','D5-D6')
title('NPCC')

subplot(2,2,3)
plot(delta_range,Pc_comp_TTYPE_delta(:,1:4),'LineWidth',2)
legend('T1-T4','D1-D4','T2-T3','D2-D3')
title('T-Type')

subplot(2,2,4)
plot(delta_range,Pc_comp_2LC_delta(:,1:2),'LineWidth',2)
legend('T1-T2','D1-D2')
title('Two Level converter')

figindx=figindx+1;


figure(figindx)
sgtitle('Switching losses varying \delta')

subplot(2,2,1)
hold on
plot(delta_range,Ps_comp_ANPCC_delta(:,1:4),'LineWidth',2)
plot(delta_range,Ps_comp_ANPCC_delta(:,9:10),'LineWidth',2)
hold off
legend('T1-T4','D1-D4','T2-T3','D2-D3','T5-T6','D5-D6')
title('ANPCC')

subplot(2,2,2)
hold on
plot(delta_range,Ps_comp_NPCC_delta(:,1:4),'LineWidth',2)
plot(delta_range,Ps_comp_NPCC_delta(:,9),'LineWidth',2)
hold off
legend('T1-T4','D1-D4','T2-T3','D2-D3','D5-D6')
title('NPCC')

subplot(2,2,3)
plot(delta_range,Ps_comp_TTYPE_delta(:,1:4),'LineWidth',2)
legend('T1-T4','D1-D4','T2-T3','D2-D3')
title('T-Type')

subplot(2,2,4)
plot(delta_range,Ps_comp_2LC_delta(:,1:2),'LineWidth',2)
legend('T1-T2','D1-D2')
title('Two Level converter')


figindx=figindx+1;

figure(figindx)
sgtitle('Semiconductor losses efficiency varying \delta')
subplot(2,2,1)
hold on
plot(delta_range,Ptot_comp_ANPCC_delta(:,1:4),'LineWidth',2)
plot(delta_range,Ptot_comp_ANPCC_delta(:,9:10),'LineWidth',2)
hold off
title('ANPCC')
legend('T1-T4','D1-D4','T2-T3','D2-D3','T5-T6','D5-D6')

subplot(2,2,2)
hold on
plot(delta_range,Ptot_comp_NPCC_delta(:,1:4),'LineWidth',2)
plot(delta_range,Ptot_comp_NPCC_delta(:,9),'LineWidth',2)
hold off
title('NPCC')
legend('T1-T4','D1-D4','T2-T3','D2-D3','D5-D6')

subplot(2,2,3)
plot(delta_range,Ptot_comp_TTYPE_delta(:,1:4),'LineWidth',2)
title('T-Type')
legend('T1-T4','D1-D4','T2-T3','D2-D3')

subplot(2,2,4)
plot(delta_range,Ptot_comp_2LC_delta(:,1:2),'LineWidth',2)
title('Two Level converter')
legend('T1-T2','D1-D2')
figindx=figindx+1;


%% Sweeping m Modulation index

ddeg=0;

f__s=12000; % f at 10 kHz

m_1=0.5;
m_2=0.95;

i=0;

for msweeping=m_1:0.01:m_2

i=i+1;

U__DC=Vac*2*sqrt(2)/sqrt(3)/msweeping;

[Pc_comp_ANPCC,Ps_comp_ANPCC,Ptot_comp_ANPCC,Ptot_ANPCC]=anpcc_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_ANPCC_m(i)=(Pac-Ptot_ANPCC)/Pac*100;

[Pc_comp_NPCC,Ps_comp_NPCC,Ptot_comp_NPCC,Ptot_NPCC]=npcc_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_NPCC_m(i)=(Pac-Ptot_NPCC)/Pac*100;

[Pc_comp_TTYPE,Ps_comp_TTYPE,Ptot_comp_TTYPE,Ptot_TTYPE]=ttype2_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_TTYPE_m(i)=(Pac-Ptot_TTYPE)/Pac*100;

[Pc_comp_2LC,Ps_comp_2LC,Ptot_comp_2LC,Ptot_2LC]=TwoLC_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_2LC_m(i)=(Pac-Ptot_2LC)/Pac*100;

end

m_range=m_1:0.01:m_2;

figure(figindx)
plot(m_range,eff_ANPCC_m,m_range,eff_NPCC_m,m_range,eff_TTYPE_m,m_range,eff_2LC_m,'LineWidth',2)
title('Semiconductor efficiency')
xlabel('Modulation index [m]')
xlim([0.05 0.95])
ylabel('Efficiency [%]')
legend('ANPCC','NPCC','T-Type','2 LC')
figindx=figindx+1;


%% Sweeping DC voltage <--> m Modulation index \ Battery pack Samsung M3-R089

ddeg=0;

f__s=12000; % f at 10 kHz

U_DC_1=774;
U_DC_2=1004;
m1_M3R089=Vac*2*sqrt(2)/sqrt(3)/U_DC_1;
m2_M3R089=Vac*2*sqrt(2)/sqrt(3)/U_DC_2;

i=0;

for U__DC=U_DC_1:1:U_DC_2

i=i+1;

[Pc_comp_ANPCC,Ps_comp_ANPCC,Ptot_comp_ANPCC,Ptot_ANPCC]=anpcc_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_ANPCC_M3R089(i)=(Pac-Ptot_ANPCC)/Pac*100;

[Pc_comp_NPCC,Ps_comp_NPCC,Ptot_comp_NPCC,Ptot_NPCC]=npcc_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_NPCC_M3R089(i)=(Pac-Ptot_NPCC)/Pac*100;

[Pc_comp_TTYPE,Ps_comp_TTYPE,Ptot_comp_TTYPE,Ptot_TTYPE]=ttype2_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_TTYPE_M3R089(i)=(Pac-Ptot_TTYPE)/Pac*100;

[Pc_comp_2LC,Ps_comp_2LC,Ptot_comp_2LC,Ptot_2LC]=TwoLC_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_2LC_M3R089(i)=(Pac-Ptot_2LC)/Pac*100;

end

U__DC_range=U_DC_1:1:U_DC_2;

delta_range=delta_1:1:delta_2;

figure(figindx)
plot(U__DC_range,eff_ANPCC_M3R089,U__DC_range,eff_NPCC_M3R089,U__DC_range,...
    eff_TTYPE_M3R089,U__DC_range,eff_2LC_M3R089,'LineWidth',2)
title('Semiconductor efficiency - U_{DC} range of Battery Rack Samsung M3-R089')
xlabel('U_{DC} [V]')
xlim([774 1004])
ylabel('Efficiency [%]')
legend('ANPCC','NPCC','T-Type','2 LC')
figindx=figindx+1;

figure(figindx)
% sgtitle('Efficiency')
subplot(1,3,1)
plot(U__DC_range,eff_ANPCC_M3R089,U__DC_range,eff_NPCC_M3R089,U__DC_range,...
    eff_TTYPE_M3R089,U__DC_range,eff_2LC_M3R089,'LineWidth',2)
xlabel('U_{DC} [V] - range of Battery Rack Samsung M3-R089')
xlim([774 1004])
ylabel('Efficiency [%]')

subplot(1,3,2)

plot(m_range,eff_ANPCC_m,m_range,eff_NPCC_m,m_range,eff_TTYPE_m,m_range,eff_2LC_m,'LineWidth',2)
xlabel('Modulation index [m]')
xlim([0.5 0.95])
ylabel('Efficiency [%]')

subplot(1,3,3)
plot(delta_range,eff_ANPCC_delta,delta_range,eff_NPCC_delta...
    ,delta_range,eff_TTYPE_delta,delta_range,eff_2LC_delta,'LineWidth',2)
xlabel('\delta [\circ]')
xlim([-180 180])
xticks([-180  -90 0 90 180])
xticklabels({'-180','-90','0','90','180'})
ylabel('Efficiency [%]')
legend('ANPCC','NPCC','T-Type','2 LC')
figindx=figindx+1;

%% Sweeping P for partial load

ddeg=0;

f__s=12000; % f at 10 kHz
U__DC=900;
Vac=400;
P_1=10000;
P_2=100000;

i=0;
for P=P_1:1:P_2
i=i+1;

Pac=P;

[Pc_comp_ANPCC,Ps_comp_ANPCC,Ptot_comp_ANPCC,Ptot_ANPCC]=anpcc_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_ANPCC_P(i)=(Pac-750-Ptot_ANPCC)/Pac*100;

[Pc_comp_NPCC,Ps_comp_NPCC,Ptot_comp_NPCC,Ptot_NPCC]=npcc_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_NPCC_P(i)=(Pac-750-Ptot_NPCC)/Pac*100;

[Pc_comp_TTYPE,Ps_comp_TTYPE,Ptot_comp_TTYPE,Ptot_TTYPE]=ttype2_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_TTYPE_P(i)=(Pac-750-Ptot_TTYPE)/Pac*100;

[Pc_comp_2LC,Ps_comp_2LC,Ptot_comp_2LC,Ptot_2LC]=TwoLC_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_2LC_P(i)=(Pac-750-Ptot_2LC)/Pac*100;
% Losses_2LC(i,:) = [Ptot_2LC  Pac   eff_2LC_P(i) ] ;
end
inverter_eff = [P_1:1:P_2  eff_2LC_P] ; 
P_range=P_1:1:P_2;
figure(figindx)
plot(P_range,eff_ANPCC_P,P_range,eff_NPCC_P,P_range,eff_TTYPE_P,P_range,eff_2LC_P)
title('Semiconductor efficiency varying P')
legend('ANPCC','NPCC','T-Type','2 LC')
ylim([90 99])
figindx=figindx+1;

figure(figindx)
hold on
plot(P_range/1000,eff_2LC_P,'b','LineWidth',1.5)
title('Semiconductor efficiency varying P')
ylim([90 98])
xlabel('AC Power [kW]')
ylabel('Efficiency [%]')
grid on
figindx=figindx+1;

figure(figindx)
hold on
plot(P_range/1000,eff_ANPCC_P,'r','LineWidth',1.5)
plot(P_range/1000,eff_NPCC_P,'g','LineWidth',1.5)
plot(P_range/1000,eff_TTYPE_P,'b','LineWidth',1.5)
plot(P_range/1000,eff_2LC_P,'k','LineWidth',1.5)
title('Semiconductor efficiency varying P')
legend('ANPCC','NPCC','T-Type','2 LC')
ylim([90 99])
xlabel('AC Power [kW]')
ylabel('Efficiency [%]')
grid on
figindx=figindx+1;

%% Sweeping P for frequency regulation

load('Pfreqreg.mat')

ddeg=0;

f__s=12000; % f at 10 kHz

Psp=Psp2(1:86400); %%% this is 1 day --> 24h*60min*60s=86400
Psp=Psp*1000000;

for i=1:1:size(Psp,1)
Pac=Psp(i);
if Pac<0
    Pac=-Pac;
    ddeg=180;
end

[Pc_comp_ANPCC,Ps_comp_ANPCC,Ptot_comp_ANPCC,Ptot_ANPCC]=anpcc_freq(U__DC,Vac,Pac,ddeg,f__s);

Ptot_ANPCC_Pfreg(i)=Ptot_ANPCC;

[Pc_comp_NPCC,Ps_comp_NPCC,Ptot_comp_NPCC,Ptot_NPCC]=npcc_freq(U__DC,Vac,Pac,ddeg,f__s);

Ptot_NPCC_Pfreg(i)=Ptot_NPCC;

[Pc_comp_TTYPE,Ps_comp_TTYPE,Ptot_comp_TTYPE,Ptot_TTYPE]=ttype2_freq(U__DC,Vac,Pac,ddeg,f__s);

Ptot_TTYPE_Pfreg(i)=Ptot_TTYPE;

[Pc_comp_2LC,Ps_comp_2LC,Ptot_comp_2LC,Ptot_2LC]=TwoLC_freq(U__DC,Vac,Pac,ddeg,f__s);

Ptot_2LC_Pfreg(i)=Ptot_2LC;

end

Energy_cycled=sum(abs(Psp));

Ptot_Pfreg_final=[sum(Ptot_ANPCC_Pfreg) sum(Ptot_NPCC_Pfreg) sum(Ptot_TTYPE_Pfreg) sum(Ptot_2LC_Pfreg)];

Ptot_Pfreg_final=Ptot_Pfreg_final/size(Psp,1);

Ptot_Pfreg_eff=(Energy_cycled-Ptot_Pfreg_final)/Energy_cycled*100;

t_range=1:1:size(Psp,1);

figure(figindx)
plot(t_range,Ptot_ANPCC_Pfreg,t_range,Ptot_NPCC_Pfreg...
    ,t_range,Ptot_TTYPE_Pfreg,t_range,Ptot_2LC_Pfreg)
title('Semiconductor losses providing frequency regulation')
legend('ANPCC','NPCC','T-Type','2 LC')
figindx=figindx+1;

figure(figindx)
bar(Ptot_Pfreg_final)
xticklabels({'ANPCC','NPCC','T-type','2LC'})
title('Semiconductor losses providing frequency regulation')
figindx=figindx+1;

figure(figindx)
bar(Ptot_Pfreg_eff)
xticklabels({'ANPCC','NPCC','T-type','2LC'})
title('Semiconductor efficiency providing frequency regulation')
ylim([90 100])
figindx=figindx+1;

%% Losses distribution

i=0;
delta_testing=[0 90 180 -90];
Pac=100000;
f__s=4000; % f at 10 kHz

while size(delta_testing,2)>=1

ddeg=delta_testing(1);

[Pc_comp_ANPCC,Ps_comp_ANPCC,Ptot_comp_ANPCC,Ptot_ANPCC]=anpcc_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_ANPCC=round((Pac-Ptot_ANPCC)/Pac*100,2);

[Pc_comp_NPCC,Ps_comp_NPCC,Ptot_comp_NPCC,Ptot_NPCC]=npcc_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_NPCC=round((Pac-Ptot_NPCC)/Pac*100,2);

[Pc_comp_TTYPE,Ps_comp_TTYPE,Ptot_comp_TTYPE,Ptot_TTYPE]=ttype2_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_TTYPE=round((Pac-Ptot_TTYPE)/Pac*100,2);

[Pc_comp_2LC,Ps_comp_2LC,Ptot_comp_2LC,Ptot_2LC]=TwoLC_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_2LC=round((Pac-Ptot_2LC)/Pac*100,2);

ddeg=round(ddeg,2);
figure(figindx)
sgtitle('Semiconductor losses distribution 1')
hold on
subplot(4,4,i*4+1)
bar(transpose([Pc_comp_ANPCC; Ps_comp_ANPCC]),'stacked')
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','T5','D5','T6','D6'})
title(sprintf('ANPCC  \\delta = %g^{\\circ} \\eta = %g%%',ddeg,eff_ANPCC))

subplot(4,4,i*4+2)
bar(transpose([Pc_comp_NPCC; Ps_comp_NPCC]),'stacked')
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','D5','D6'})
title(sprintf('NPCC  \\delta = %g^{\\circ} \\eta = %g%%',ddeg,eff_NPCC))

subplot(4,4,i*4+3)
bar(transpose([Pc_comp_TTYPE; Ps_comp_TTYPE]),'stacked')
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4'})
title(sprintf('Ttype  \\delta = %g^{\\circ} \\eta = %g%%',ddeg,eff_TTYPE))

subplot(4,4,i*4+4)
bar(transpose([Pc_comp_2LC; Ps_comp_2LC]),'stacked')
xticklabels({'T1','D1','T2','D2'})
title(sprintf('2LC \\delta = %g^{\\circ}  \\eta = %g%%',ddeg,eff_2LC))



i=i+1;

delta_testing(1)=[];
end
figindx=figindx+1;

%% Losses distribution 2

ddeg=0;
Pac=100000;
f__s=12000; % f at 10 kHz

[Pc_comp_ANPCC,Ps_comp_ANPCC,Ptot_comp_ANPCC,Ptot_ANPCC]=anpcc_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_ANPCC=round((Pac-Ptot_ANPCC)/Pac*100,2);

[Pc_comp_NPCC,Ps_comp_NPCC,Ptot_comp_NPCC,Ptot_NPCC]=npcc_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_NPCC=round((Pac-Ptot_NPCC)/Pac*100,2);

[Pc_comp_TTYPE,Ps_comp_TTYPE,Ptot_comp_TTYPE,Ptot_TTYPE]=ttype2_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_TTYPE=round((Pac-Ptot_TTYPE)/Pac*100,2);

[Pc_comp_2LC,Ps_comp_2LC,Ptot_comp_2LC,Ptot_2LC]=TwoLC_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_2LC=round((Pac-Ptot_2LC)/Pac*100,2);

figure(figindx)
% sgtitle('Semiconductor losses distribution 1')
subplot(1,4,1)
bar(transpose([Pc_comp_2LC; Ps_comp_2LC/3 ; Ps_comp_2LC/3*2]),0.45,'stacked')
ylabel('Losses [W]')
xticklabels({'T1','D1','T2','D2'})
% title(sprintf('2LC - \\eta = %g%%',eff_2LC))
title('2LC')
ylim([0 300])
grid on

subplot(1,4,2)
bar(transpose([Pc_comp_TTYPE; Ps_comp_TTYPE/3 ; Ps_comp_TTYPE/3*2]),'stacked')
ylabel('Losses [W]')
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4'})
% title(sprintf('Ttype - \\eta = %g%%',eff_TTYPE))
title('Ttype')
ylim([0 300])
grid on

subplot(1,4,3)
bar(transpose([Pc_comp_NPCC; Ps_comp_NPCC/3 ; Ps_comp_NPCC/3*2]),'stacked')
ylabel('Losses [W]')
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','D5','D6'})
% title(sprintf('NPCC - \\eta = %g%%',eff_NPCC))
title('NPCC')
ylim([0 300])
grid on

subplot(1,4,4)
bar(transpose([Pc_comp_ANPCC; Ps_comp_ANPCC/3 ; Ps_comp_ANPCC/3*2]),'stacked')
ylabel('Losses [W]')
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','T5','D5','T6','D6'})
% title(sprintf('ANPCC - \\eta = %g%%',eff_ANPCC))
title('ANPCC')
legend('Conduction Losses','Switching Losses 4 kHz','Additional Switching Losses at 12 kHz')
figindx=figindx+1;
ylim([0 300])
grid on

figure(figindx)
% sgtitle('Semiconductor losses distribution 1')
subplot(1,4,1)
bar(transpose([Pc_comp_2LC; Ps_comp_2LC]),0.45,'stacked')
ylabel('Losses [W]')
xticklabels({'T1','D1','T2','D2'})
title(sprintf('2LC - \\eta = %g%%',eff_2LC))
ylim([0 300])
grid on

subplot(1,4,2)
bar(transpose([Pc_comp_TTYPE; Ps_comp_TTYPE]),'stacked')
ylabel('Losses [W]')
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4'})
title(sprintf('Ttype - \\eta = %g%%',eff_TTYPE))
ylim([0 300])
grid on

subplot(1,4,3)
bar(transpose([Pc_comp_NPCC; Ps_comp_NPCC]),'stacked')
ylabel('Losses [W]')
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','D5','D6'})
title(sprintf('NPCC - \\eta = %g%%',eff_NPCC))
ylim([0 300])
grid on

subplot(1,4,4)
bar(transpose([Pc_comp_ANPCC; Ps_comp_ANPCC]),'stacked')
ylabel('Losses [W]')
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','T5','D5','T6','D6'})
title(sprintf('ANPCC - \\eta = %g%%',eff_ANPCC))
legend('Conduction Losses','Switching Losses')
figindx=figindx+1;
ylim([0 300])
grid on

figure(figindx)
% sgtitle('Semiconductor losses distribution 1')
subplot(2,2,1)
bar(transpose([Pc_comp_2LC; Ps_comp_2LC]),'stacked')
ylabel('Losses [W]')
xticklabels({'T1','D1','T2','D2'})
title(sprintf('2LC - \\eta = %g%%',eff_2LC))
ylim([0 300])
grid on

subplot(2,2,2)
bar(transpose([Pc_comp_TTYPE; Ps_comp_TTYPE]),'stacked')
ylabel('Losses [W]')
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4'})
title(sprintf('Ttype - \\eta = %g%%',eff_TTYPE))
ylim([0 300])
grid on

subplot(2,2,3)
bar(transpose([Pc_comp_NPCC; Ps_comp_NPCC]),'stacked')
ylabel('Losses [W]')
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','D5','D6'})
title(sprintf('NPCC - \\eta = %g%%',eff_NPCC))
ylim([0 300])
grid on

subplot(2,2,4)
bar(transpose([Pc_comp_ANPCC; Ps_comp_ANPCC]),'stacked')
ylabel('Losses [W]')
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','T5','D5','T6','D6'})
title(sprintf('ANPCC - \\eta = %g%%',eff_ANPCC))
legend('Conduction Losses','Switching Losses')
figindx=figindx+1;
ylim([0 300])
grid on

%% Losses distribution top part of circl

i=0;
delta_testing=[-acosd(0.8) 0 acosd(0.8)];

f__s=12000; % f at 20 kHz

while size(delta_testing,2)>=1  
    
ddeg=delta_testing(1);

[Pc_comp_ANPCC,Ps_comp_ANPCC,Ptot_comp_ANPCC,Ptot_ANPCC]=anpcc_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_ANPCC=round((Pac-Ptot_ANPCC)/Pac*100,2);

[Pc_comp_NPCC,Ps_comp_NPCC,Ptot_comp_NPCC,Ptot_NPCC]=npcc_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_NPCC=round((Pac-Ptot_NPCC)/Pac*100,2);

[Pc_comp_TTYPE,Ps_comp_TTYPE,Ptot_comp_TTYPE,Ptot_TTYPE]=ttype2_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_TTYPE=round((Pac-Ptot_TTYPE)/Pac*100,2);

[Pc_comp_2LC,Ps_comp_2LC,Ptot_comp_2LC,Ptot_2LC]=TwoLC_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_2LC=round((Pac-Ptot_2LC)/Pac*100,2);

ddeg=round(ddeg,2);
figure(figindx)
sgtitle('Semiconductor losses distribution 2')
hold on
subplot(3,4,i*4+1)
bar(transpose([Pc_comp_ANPCC; Ps_comp_ANPCC]),'stacked')
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','T5','D5','T6','D6'})
title(sprintf('ANPCC  \\delta = %g^{\\circ} \\eta = %g%%',ddeg,eff_ANPCC))

subplot(3,4,i*4+2)
bar(transpose([Pc_comp_NPCC; Ps_comp_NPCC]),'stacked')
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','D5','D6'})
title(sprintf('NPCC  \\delta = %g^{\\circ} \\eta = %g%%',ddeg,eff_NPCC))

subplot(3,4,i*4+3)
bar(transpose([Pc_comp_TTYPE; Ps_comp_TTYPE]),'stacked')
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4'})
title(sprintf('Ttype  \\delta = %g^{\\circ} \\eta = %g%%',ddeg,eff_TTYPE))

subplot(3,4,i*4+4)
bar(transpose([Pc_comp_2LC; Ps_comp_2LC]),'stacked')
xticklabels({'T1','D1','T2','D2'})
title(sprintf('2LC \\delta = %g^{\\circ}  \\eta = %g%%',ddeg,eff_2LC))

i=i+1;

delta_testing(1)=[];
end
figindx=figindx+1;
