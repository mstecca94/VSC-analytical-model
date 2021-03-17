addpath(genpath('.\Functions'))
ccc
%% Design parameters
U__DC=900;
Vac=400;
Pac=100000; %% 100 kW
ddeg=0;

%% Sweeping frequency - Inverter operation delta=0
f__s1=4000;
f__s2=16000;
f__s_range = f__s1:1000:f__s2 ;
eff_2LC=zeros(size(f__s_range)) ;
eff_TTYPE=zeros(size(f__s_range)) ;
eff_NPCC=zeros(size(f__s_range)) ;
eff_ANPCC=zeros(size(f__s_range)) ;

i=0;
for f__s=f__s1:1000:f__s2
i=i+1;
[~,~,~,Ptot_ANPCC]=anpcc_f(U__DC,Vac,Pac,ddeg,f__s);

eff_ANPCC(i)=(Pac-Ptot_ANPCC)/Pac*100;

[~,~,~,Ptot_NPCC]=npcc_f(U__DC,Vac,Pac,ddeg,f__s);

eff_NPCC(i)=(Pac-Ptot_NPCC)/Pac*100;

[~,~,~,Ptot_TTYPE]=ttype_f(U__DC,Vac,Pac,ddeg,f__s);

eff_TTYPE(i)=(Pac-Ptot_TTYPE)/Pac*100;

[~,~,~,Ptot_2LC]=twoL_f(U__DC,Vac,Pac,ddeg,f__s);

eff_2LC(i)=(Pac-Ptot_2LC)/Pac*100;

end

eff_ANPCC_inv=eff_ANPCC;
eff_NPCC_inv=eff_NPCC;
eff_TTYPE_inv=eff_TTYPE;
eff_2LC_inv=eff_2LC;

%% Sweeping frequency - Inverter operation delta=0
i=0;
ddeg=180;

eff_2LC=zeros(size(f__s_range)) ;
eff_TTYPE=zeros(size(f__s_range)) ;
eff_NPCC=zeros(size(f__s_range)) ;
eff_ANPCC=zeros(size(f__s_range)) ;
for f__s=f__s_range
i=i+1;

[~,~,~,Ptot_ANPCC]=anpcc_f(U__DC,Vac,Pac,ddeg,f__s);
eff_ANPCC(i)=(Pac-Ptot_ANPCC)/Pac*100;

[~,~,~,Ptot_NPCC]=npcc_f(U__DC,Vac,Pac,ddeg,f__s);
eff_NPCC(i)=(Pac-Ptot_NPCC)/Pac*100;

[~,~,~,Ptot_TTYPE]=ttype_f(U__DC,Vac,Pac,ddeg,f__s);
eff_TTYPE(i)=(Pac-Ptot_TTYPE)/Pac*100;

[~,~,~,Ptot_2LC]=twoL_f(U__DC,Vac,Pac,ddeg,f__s);
eff_2LC(i)=(Pac-Ptot_2LC)/Pac*100;

end
f__s_range=f__s_range/1000;
%%

figure;
subplot(121)
hold on
plot(f__s_range,eff_ANPCC_inv,'r','LineWidth',2)
plot(f__s_range,eff_NPCC_inv,'b','LineWidth',2)
plot(f__s_range,eff_TTYPE_inv,'g','LineWidth',2)
plot(f__s_range,eff_2LC_inv,'k','LineWidth',2)
grid on
xlabel('Switching Frequency [kHz]')
ylabel('Efficiency [%]')
title('Rectifier operation')
xlim([f__s1/1000 f__s2/1000])
title('Inverter operation')

subplot(122)
hold on
plot(f__s_range,eff_ANPCC,'r','LineWidth',2)
plot(f__s_range,eff_NPCC,'b','LineWidth',2)
plot(f__s_range,eff_TTYPE,'g','LineWidth',2)
plot(f__s_range,eff_2LC,'k','LineWidth',2)
grid on
legend('ANPCC','NPCC','T-Type','2 LC')
xlabel('Switching Frequency [kHz]')
ylabel('Efficiency [%]')
title('Rectifier operation')
xlim([f__s1/1000 f__s2/1000])


%% Sweeping delta - Inverter operation delta=0

f__s=12000; % f at 10 kHz
delta_1=-180;
delta_2=180;
d_range = delta_1:1:delta_2 ;

eff_2LC_delta=zeros(size(d_range)) ;
eff_TTYPE_delta=zeros(size(d_range)) ;
eff_NPCC_delta=zeros(size(d_range)) ;
eff_ANPCC_delta=zeros(size(d_range)) ;

i=0;
for d=d_range
i=i+1;
ddeg=d;
[Pc_comp_ANPCC,Ps_comp_ANPCC,Ptot_comp_ANPCC,Ptot_ANPCC]=anpcc_f(U__DC,Vac,Pac,ddeg,f__s);

Pc_comp_ANPCC_delta(i,:)=Pc_comp_ANPCC;
Ps_comp_ANPCC_delta(i,:)=Ps_comp_ANPCC;
Ptot_comp_ANPCC_delta(i,:)=Ptot_comp_ANPCC;

eff_ANPCC_delta(i)=(Pac-Ptot_ANPCC)/Pac*100;

[Pc_comp_NPCC,Ps_comp_NPCC,Ptot_comp_NPCC,Ptot_NPCC]=npcc_f(U__DC,Vac,Pac,ddeg,f__s);

Pc_comp_NPCC_delta(i,:)=Pc_comp_NPCC;
Ps_comp_NPCC_delta(i,:)=Ps_comp_NPCC;
Ptot_comp_NPCC_delta(i,:)=Ptot_comp_NPCC;

eff_NPCC_delta(i)=(Pac-Ptot_NPCC)/Pac*100;

[Pc_comp_TTYPE,Ps_comp_TTYPE,Ptot_comp_TTYPE,Ptot_TTYPE]=ttype_f(U__DC,Vac,Pac,ddeg,f__s);

Pc_comp_TTYPE_delta(i,:)=Pc_comp_TTYPE;
Ps_comp_TTYPE_delta(i,:)=Ps_comp_TTYPE;
Ptot_comp_TTYPE_delta(i,:)=Ptot_comp_TTYPE;

eff_TTYPE_delta(i)=(Pac-Ptot_TTYPE)/Pac*100;

[Pc_comp_2LC,Ps_comp_2LC,Ptot_comp_2LC,Ptot_2LC]=twoL_f(U__DC,Vac,Pac,ddeg,f__s);

Pc_comp_2LC_delta(i,:)=Pc_comp_2LC;
Ps_comp_2LC_delta(i,:)=Ps_comp_2LC;
Ptot_comp_2LC_delta(i,:)=Ptot_comp_2LC;

eff_2LC_delta(i)=(Pac-Ptot_2LC)/Pac*100;

end

delta_range=delta_1:1:delta_2;

figure;
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

figure;
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


figure;
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

%% Sweeping m Modulation index

ddeg=0;

f__s=12000; % f at 10 kHz

m_1=0.5;
m_2=0.95;
m_range = m_1:0.01:m_2 ;

eff_2LC_m=zeros(size(m_range)) ;
eff_TTYPE_m=zeros(size(m_range)) ;
eff_NPCC_m=zeros(size(m_range)) ;
eff_ANPCC_m=zeros(size(m_range)) ;

i=0;
for msweeping=m_range

i=i+1;

U__DC=Vac*2*sqrt(2)/sqrt(3)/msweeping;

[~,~,~,Ptot_ANPCC]=anpcc_f(U__DC,Vac,Pac,ddeg,f__s);

eff_ANPCC_m(i)=(Pac-Ptot_ANPCC)/Pac*100;

[~,~,~,Ptot_NPCC]=npcc_f(U__DC,Vac,Pac,ddeg,f__s);

eff_NPCC_m(i)=(Pac-Ptot_NPCC)/Pac*100;

[~,~,~,Ptot_TTYPE]=ttype_f(U__DC,Vac,Pac,ddeg,f__s);

eff_TTYPE_m(i)=(Pac-Ptot_TTYPE)/Pac*100;

[~,~,~,Ptot_2LC]=twoL_f(U__DC,Vac,Pac,ddeg,f__s);

eff_2LC_m(i)=(Pac-Ptot_2LC)/Pac*100;

end


%% Sweeping DC voltage <--> m Modulation index \ Battery pack Samsung M3-R089

ddeg=0;

f__s=12000; % f at 10 kHz

U_DC_1=774;
U_DC_2=1004;
m1_M3R089=Vac*2*sqrt(2)/sqrt(3)/U_DC_1;
m2_M3R089=Vac*2*sqrt(2)/sqrt(3)/U_DC_2;

Udc_range = U_DC_1:1:U_DC_2 ;
eff_2LC_M3R089=zeros(size(Udc_range)) ;
eff_TTYPE_M3R089=zeros(size(Udc_range)) ;
eff_NPCC_M3R089=zeros(size(Udc_range)) ;
eff_ANPCC_M3R089=zeros(size(Udc_range)) ;

i=0;
for U__DC= Udc_range

i=i+1;

[~,~,~,Ptot_ANPCC]=anpcc_f(U__DC,Vac,Pac,ddeg,f__s);

eff_ANPCC_M3R089(i)=(Pac-Ptot_ANPCC)/Pac*100;

[~,~,~,Ptot_NPCC]=npcc_f(U__DC,Vac,Pac,ddeg,f__s);

eff_NPCC_M3R089(i)=(Pac-Ptot_NPCC)/Pac*100;

[~,~,~,Ptot_TTYPE]=ttype_f(U__DC,Vac,Pac,ddeg,f__s);

eff_TTYPE_M3R089(i)=(Pac-Ptot_TTYPE)/Pac*100;

[~,~,~,Ptot_2LC]=twoL_f(U__DC,Vac,Pac,ddeg,f__s);

eff_2LC_M3R089(i)=(Pac-Ptot_2LC)/Pac*100;

end


delta_range=delta_1:1:delta_2;

figure;
sgtitle('Efficiency varying different parameters')
subplot(1,3,1)
plot(Udc_range,eff_ANPCC_M3R089,Udc_range,eff_NPCC_M3R089,Udc_range,...
    eff_TTYPE_M3R089,Udc_range,eff_2LC_M3R089,'LineWidth',2)
xlabel('U_{DC} [V]')
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


%% Losses distribution 1

ddeg=180;
Pac=100000;
f__s=12000; % f at 10 kHz

[Pc_comp_ANPCC,Ps_comp_ANPCC,Ptot_comp_ANPCC,Ptot_ANPCC]=anpcc_f(U__DC,Vac,Pac,ddeg,f__s);

eff_ANPCC=round((Pac-Ptot_ANPCC)/Pac*100,2);

[Pc_comp_NPCC,Ps_comp_NPCC,Ptot_comp_NPCC,Ptot_NPCC]=npcc_f(U__DC,Vac,Pac,ddeg,f__s);

eff_NPCC=round((Pac-Ptot_NPCC)/Pac*100,2);

[Pc_comp_TTYPE,Ps_comp_TTYPE,Ptot_comp_TTYPE,Ptot_TTYPE]=ttype_f(U__DC,Vac,Pac,ddeg,f__s);

eff_TTYPE=round((Pac-Ptot_TTYPE)/Pac*100,2);

[Pc_comp_2LC,Ps_comp_2LC,Ptot_comp_2LC,Ptot_2LC]=twoL_f(U__DC,Vac,Pac,ddeg,f__s);

eff_2LC=round((Pac-Ptot_2LC)/Pac*100,2);

figure;
sgtitle('Loss distribtution at inverter operation')
subplot(221)
bar(transpose([Pc_comp_2LC; Ps_comp_2LC/3 ; Ps_comp_2LC/3*2]),0.45,'stacked')
ylabel('Losses [W]')
xticklabels({'T1','D1','T2','D2'})
title('2LC')
grid on

subplot(222)
bar(transpose([Pc_comp_TTYPE; Ps_comp_TTYPE/3 ; Ps_comp_TTYPE/3*2]),'stacked')
ylabel('Losses [W]')
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4'})
title('Ttype')
grid on

subplot(223)
bar(transpose([Pc_comp_NPCC; Ps_comp_NPCC/3 ; Ps_comp_NPCC/3*2]),'stacked')
ylabel('Losses [W]')
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','D5','D6'})
title('NPCC')
grid on

subplot(224)
bar(transpose([Pc_comp_ANPCC; Ps_comp_ANPCC/3 ; Ps_comp_ANPCC/3*2]),'stacked')
ylabel('Losses [W]')
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','T5','D5','T6','D6'})
title('ANPCC')
legend('Conduction Losses','Switching Losses 4 kHz','Additional Switching Losses at 12 kHz')
grid on

%% Losses distribution 2

ddeg=0;
Pac=100000;
f__s=12000; % f at 10 kHz

[Pc_comp_ANPCC,Ps_comp_ANPCC,Ptot_comp_ANPCC,Ptot_ANPCC]=anpcc_f(U__DC,Vac,Pac,ddeg,f__s);

eff_ANPCC=round((Pac-Ptot_ANPCC)/Pac*100,2);

[Pc_comp_NPCC,Ps_comp_NPCC,Ptot_comp_NPCC,Ptot_NPCC]=npcc_f(U__DC,Vac,Pac,ddeg,f__s);

eff_NPCC=round((Pac-Ptot_NPCC)/Pac*100,2);

[Pc_comp_TTYPE,Ps_comp_TTYPE,Ptot_comp_TTYPE,Ptot_TTYPE]=ttype_f(U__DC,Vac,Pac,ddeg,f__s);

eff_TTYPE=round((Pac-Ptot_TTYPE)/Pac*100,2);

[Pc_comp_2LC,Ps_comp_2LC,Ptot_comp_2LC,Ptot_2LC]=twoL_f(U__DC,Vac,Pac,ddeg,f__s);

eff_2LC=round((Pac-Ptot_2LC)/Pac*100,2);

figure;
sgtitle('Loss distribtution at rectifier operation')
subplot(221)
bar(transpose([Pc_comp_2LC; Ps_comp_2LC/3 ; Ps_comp_2LC/3*2]),0.45,'stacked')
ylabel('Losses [W]')
xticklabels({'T1','D1','T2','D2'})
title('2LC')
grid on

subplot(222)
bar(transpose([Pc_comp_TTYPE; Ps_comp_TTYPE/3 ; Ps_comp_TTYPE/3*2]),'stacked')
ylabel('Losses [W]')
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4'})
title('Ttype')
grid on

subplot(223)
bar(transpose([Pc_comp_NPCC; Ps_comp_NPCC/3 ; Ps_comp_NPCC/3*2]),'stacked')
ylabel('Losses [W]')
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','D5','D6'})
title('NPCC')
grid on

subplot(224)
bar(transpose([Pc_comp_ANPCC; Ps_comp_ANPCC/3 ; Ps_comp_ANPCC/3*2]),'stacked')
ylabel('Losses [W]')
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','T5','D5','T6','D6'})
title('ANPCC')
legend('Conduction Losses','Switching Losses 4 kHz','Additional Switching Losses at 12 kHz')

grid on

%% Heatmap

f__s1=4000;
f__s2=16000;
P_1=-100000;
P_2=100000;
P_range=P_1:10000:P_2;
f_range=f__s1:1000:f__s2;
P_range(find(P_range==0))=[];
Paux=0;

ii=0;
for f__s=f_range
ii=ii+1;
i=0;
for P=P_range
i=i+1;
ddeg=0;
if P<0
    P=abs(P);
    ddeg=180;
end
Pac=P;

if P==0
    eff_ANPCC_P(ii,i)=100;
    eff_NPCC_P(ii,i)=100;
    eff_TTYPE_P(ii,i)=100;
    eff_2LC_P(ii,i)=100;
    continue
end
[Pc_comp_ANPCC,Ps_comp_ANPCC,Ptot_comp_ANPCC,Ptot_ANPCC]=anpcc_f(U__DC,Vac,Pac,ddeg,f__s);

eff_ANPCC_P(ii,i)=(Pac-Paux-Ptot_ANPCC)/Pac*100;

[Pc_comp_NPCC,Ps_comp_NPCC,Ptot_comp_NPCC,Ptot_NPCC]=npcc_f(U__DC,Vac,Pac,ddeg,f__s);

eff_NPCC_P(ii,i)=(Pac-Paux-Ptot_NPCC)/Pac*100;

[Pc_comp_TTYPE,Ps_comp_TTYPE,Ptot_comp_TTYPE,Ptot_TTYPE]=ttype_f(U__DC,Vac,Pac,ddeg,f__s);

eff_TTYPE_P(ii,i)=(Pac-Paux-Ptot_TTYPE)/Pac*100;

[Pc_comp_2LC,Ps_comp_2LC,Ptot_comp_2LC,Ptot_2LC]=twoL_f(U__DC,Vac,Pac,ddeg,f__s);

eff_2LC_P(ii,i)=(Pac-Paux-Ptot_2LC)/Pac*100;


end

end
figure;
subplot(1,2,1)
h = heatmap(P_range/1000,f_range/1000,eff_TTYPE_P-eff_2LC_P);
xlabel('Output Power [%]')
ylabel('Switching Frequency [kHz]')
colormap(jet);
title('Efficiency difference [%] - 3L-TType vs 2LC')
caxis([0 1.5]);
h.GridVisible = 'off' ;

subplot(1,2,2)
h = heatmap(P_range/1000,f_range/1000,eff_NPCC_P-eff_2LC_P);
xlabel('Output Power [%] \leftarrow \rightarrow')
ylabel('Switching Frequency [kHz]')
colormap(jet)
title('Efficiency difference [%] - 3L-NPC vs 2LC')
caxis([0 1.5]);
h.GridVisible = 'off' ;

%% Heatsink

f__s1=4000;
f__s2=16000;
f_range=f__s1:1000:f__s2;
Paux = 0 ;

eff_2LC_hs=zeros(size(f_range)) ;
eff_TTYPE_hs=zeros(size(f_range)) ;
eff_NPCC_hs=zeros(size(f_range)) ;

Rs_2LC_P=zeros(size(f_range)) ;
Rs_TTYPE_P=zeros(size(f_range)) ;
Rs_NPCC_P=zeros(size(f_range)) ;

i=0;
for f__s=f_range
i=i+1;

%%%%%%%%%%%% NPCC
[Pc_comp_NPCC,Ps_comp_NPCC,Ptot_comp_NPCC,Ptot_NPCC]=npcc_f(U__DC,Vac,Pac,ddeg,f__s);
eff_NPCC_hs(i)=(Pac-Paux-Ptot_NPCC)/Pac*100;

Semicond.Type = 'SEMiX305MLI07E4_NPCC_650V' ;
Semicond.topology = '3L_NPC' ;
[ Rs ]= heat_sink_design ( Ptot_comp_NPCC , Semicond );
Rs_NPCC_P(i)=Rs;

%%%%%%%%%%%% TTYPE
[Pc_comp_TTYPE,Ps_comp_TTYPE,Ptot_comp_TTYPE,Ptot_TTYPE]=ttype_f(U__DC,Vac,Pac,ddeg,f__s);
eff_TTYPE_hs(i)=(Pac-Paux-Ptot_TTYPE)/Pac*100;

Semicond.Type = 'SEMiX305TMLI12E4B' ;
Semicond.topology = '3L_TTYPE' ;
[ Rs ]= heat_sink_design ( Ptot_comp_TTYPE , Semicond );
Rs_TTYPE_P(i)=Rs;

%%%%%%%%%% 2LC
[Pc_comp_2LC,Ps_comp_2LC,Ptot_comp_2LC,Ptot_2LC]=twoL_f(U__DC,Vac,Pac,ddeg,f__s);
eff_2LC_hs(i)=(Pac-Paux-Ptot_2LC)/Pac*100;

Semicond.Type = 'FF300R12KE4' ;
Semicond.topology = '2LC' ;
[ Rs ]= heat_sink_design ( Ptot_comp_2LC , Semicond );
Rs_2LC_P(i)=Rs;


end
width=2;  

figure;
hold on
plot(f_range/1000,Rs_NPCC_P,'r','Linewidth',width)
plot(f_range/1000,Rs_TTYPE_P,'k','Linewidth',width)
plot(f_range/1000,Rs_2LC_P,'b','Linewidth',width)
legend('NPCC','TTYPE','2LC-Si')
xlabel('Switching frequency [kHz]')
ylabel('Rs [K/W]')
grid on
title('HeatSink Thermal Resistance')