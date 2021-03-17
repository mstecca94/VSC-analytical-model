%% Design parameters
% ccc
clear all
U__DC=900;
Vac=400;
ddeg=0;
figindx=1;

%% Sweeping frequency - P

f__s1=4000;
f__s2=16000;
P_1=-100000;
P_2=100000;
P_range=P_1:10000:P_2;
f_range=f__s1:1000:f__s2;
P_range(find(P_range==0))=[];
Paux=0;

i=0;
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
[Pc_comp_ANPCC,Ps_comp_ANPCC,Ptot_comp_ANPCC,Ptot_ANPCC]=anpcc_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_ANPCC_P(ii,i)=(Pac-Paux-Ptot_ANPCC)/Pac*100;

[Pc_comp_NPCC,Ps_comp_NPCC,Ptot_comp_NPCC,Ptot_NPCC]=npcc_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_NPCC_P(ii,i)=(Pac-Paux-Ptot_NPCC)/Pac*100;

[Pc_comp_TTYPE,Ps_comp_TTYPE,Ptot_comp_TTYPE,Ptot_TTYPE]=ttype2_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_TTYPE_P(ii,i)=(Pac-Paux-Ptot_TTYPE)/Pac*100;

[Pc_comp_2LC,Ps_comp_2LC,Ptot_comp_2LC,Ptot_2LC]=TwoLC_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_2LC_P(ii,i)=(Pac-Paux-Ptot_2LC)/Pac*100;

[Pc_comp_2LC_SiC,Ps_comp_2LC_SiC,Ptot_comp_2LC_SiC,Ptot_2LC_SiC]...
    =TwoLC_freqSiC(U__DC,Vac,Pac,ddeg,f__s);

eff_2LC_SiC_P(ii,i)=(Pac-Paux-Ptot_2LC_SiC)/Pac*100;

[Pc_comp_2LC_SiC_HP,Ps_comp_2LC_SiC_HP,Ptot_comp_2LC_SiC_HP,Ptot_2LC_SiC_HP]...
    =TwoLC_freqSiC_HP(U__DC,Vac,Pac,ddeg,f__s);

% [ HeatSink ]= heat_sink_design ( Ptot_comp_2LC_SiC_HP  )

eff_2LC_SiC_HP_P(ii,i)=(Pac-Paux-Ptot_2LC_SiC_HP)/Pac*100;

end

end

%%

figure;
hold on
plot(f_range/1000,eff_ANPCC_P(:,1))
plot(f_range/1000,eff_NPCC_P(:,1))
plot(f_range/1000,eff_TTYPE_P(:,1))
plot(f_range/1000,eff_2LC_P(:,1))
plot(f_range/1000,eff_2LC_SiC_P(:,1))
plot(f_range/1000,eff_2LC_SiC_HP_P(:,1))
% plot(f__s/1000,n*100)
legend('ANPCC','NPCC','TTYPE','2LC-Si','2LC-SiC','2LC-SiC-HP','2LC-2xSi')

widht=2;
figure;
hold on
plot(f_range/1000,eff_NPCC_P(:,1),'r','Linewidth',widht)
plot(f_range/1000,eff_TTYPE_P(:,1),'r--','Linewidth',widht)
plot(f_range/1000,eff_2LC_P(:,1),'b','Linewidth',widht)
plot(f_range/1000,eff_2LC_SiC_P(:,1),'g','Linewidth',widht)
% plot(fs/1000,n(1,:)*100,'b--','Linewidth',widht)
% plot(fs/1000,n(2,:)*100,'m--x','Linewidth',widht)
% plot(fs/1000,n(3,:)*100,'m--','Linewidth',widht)
% plot(fs/1000,n(4,:)*100,'g--','Linewidth',widht)
legend('NPCC','TTYPE','2LC-Si','2LC-SiC','2LC-2xSi','2LC-Si+4xSiC','2LC-SiC+4xSi','2LC-2xSiC')
xlabel('Switching frequency [kHz]')
ylabel('Efficienccy [%]')
grid on
%%
figure;
subplot(1,2,1)
h = heatmap(P_range/1000,f_range/1000,eff_TTYPE_P-eff_2LC_P)
xlabel('Output Power [%]')
ylabel('Switching Frequency [kHz]')
colormap(jet);
title('Efficiency difference [%] - 3L-TType vs 2LC')
caxis([0 1.5]);
h.GridVisible = 'off' ;

subplot(1,2,2)
h = heatmap(P_range/1000,f_range/1000,eff_NPCC_P-eff_2LC_P)
xlabel('Output Power [%] \leftarrow \rightarrow')
ylabel('Switching Frequency [kHz]')
colormap(jet)
title('Efficiency difference [%] - 3L-NPC vs 2LC')
caxis([0 1.5]);
h.GridVisible = 'off' ;
%%
npc=eff_NPCC_P;
ttype=eff_TTYPE_P;
twolc=eff_2LC_P;
twolc_SiC=eff_2LC_SiC_P;

% save('eff.mat','npc','ttype','twolc','twolc_SiC','P_range','f_range')
%%
figure;
subplot(1,3,1)
h = heatmap(P_range/1000,f_range/1000,eff_TTYPE_P-eff_2LC_P)
xlabel('Output Power [%]')
ylabel('Switching Frequency [kHz]')
colormap(jet);
title('Efficiency difference [%] - 3L-TType vs 2LC')
caxis([0 1.5]);
h.GridVisible = 'off' ;

subplot(1,3,3)
h = heatmap(P_range/1000,f_range/1000,eff_NPCC_P-eff_2LC_P)
xlabel('Output Power [%] \leftarrow \rightarrow')
ylabel('Switching Frequency [kHz]')
colormap(jet)
title('Efficiency difference [%] - 3L-NPC vs 2LC')
caxis([0 1.5]);
h.GridVisible = 'off' ;

subplot(1,3,2)
h = heatmap(P_range/1000,f_range/1000,eff_NPCC_P-eff_TTYPE_P)
xlabel('Output Power [%] \leftarrow \rightarrow')
ylabel('Switching Frequency [kHz]')
colormap(jet)
title('Efficiency difference [%] - 3L-NPC vs 3L-TType')
caxis([-0.25 0.5]);
h.GridVisible = 'off' ;

%%
diff = eff_NPCC_P-eff_TTYPE_P ;

[a,b]=find(abs(diff)<0.01);
% 
% figure;
% subplot(1,2,1)
% h = heatmap(P_range/1000,f_range/1000,eff_TTYPE_P-eff_2LC_P)
% xlabel('Output Power [%]')
% ylabel('Switching Frequency [kHz]')
% colormap(jet);
% title('Efficiency difference [%] - 3L-TType vs 2LC')
% caxis([0 1.5]);
% h.GridVisible = 'off' ;
% 
% subplot(1,2,2)
% h = heatmap(P_range/1000,f_range/1000,eff_NPCC_P-eff_2LC_P)
% xlabel('Output Power [%]')
% ylabel('Switching Frequency [kHz]')
% colormap(jet)
% title('Efficiency difference [%] - 3L-NPC vs 2LC')
% caxis([0 1.5]);
% h.GridVisible = 'off' ;

%%
% 
% figure;
% h = heatmap(P_range/1000,f_range/1000,eff_TTYPE_P-eff_2LC_P)
% xlabel('Output Power [%]')
% ylabel('Switching Frequency [kHz]')
% colormap(jet);
% title('Efficiency difference [%] - 3L-TType vs 2LC')
% caxis([0 1.5]);
% h.GridVisible = 'off' ;
% 
% 
% figure;
% h = heatmap(P_range/1000,f_range/1000,eff_NPCC_P-eff_2LC_P)
% xlabel('Output Power [%]')
% ylabel('Switching Frequency [kHz]')
% colormap(jet)
% title('Efficiency difference [%] - 3L-NPC vs 2LC')
% caxis([0 1.5]);
% h.GridVisible = 'off' ;
