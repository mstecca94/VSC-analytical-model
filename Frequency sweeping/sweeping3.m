%% Design parameters
ccc

U__DC=900;
Vac=400;
ddeg=0;
figindx=1;
Pac=100000;

%% Sweeping frequency - P

f__s1=4000;
f__s2=16000;
f_range=f__s1:1000:f__s2;
Paux = 0 ;
i=0;

for f__s=f_range
i=i+1;

%%%%%%%%%%%% NPCC
[Pc_comp_NPCC,Ps_comp_NPCC,Ptot_comp_NPCC,Ptot_NPCC]=npcc_freq(U__DC,Vac,Pac,ddeg,f__s);
eff_NPCC_P(i)=(Pac-Paux-Ptot_NPCC)/Pac*100;

Semicond.Type = 'SEMiX305MLI07E4_NPCC_650V' ;
Semicond.topology = '3L_NPC' ;
[ Rs ]= heat_sink_design ( Ptot_comp_NPCC , Semicond );
Rs_NPCC_P(i)=Rs;

%%%%%%%%%%%% TTYPE
[Pc_comp_TTYPE,Ps_comp_TTYPE,Ptot_comp_TTYPE,Ptot_TTYPE]=ttype2_freq(U__DC,Vac,Pac,ddeg,f__s);
eff_TTYPE_P(i)=(Pac-Paux-Ptot_TTYPE)/Pac*100;

Semicond.Type = 'SEMiX305TMLI12E4B' ;
Semicond.topology = '3L_TTYPE' ;
[ Rs ]= heat_sink_design ( Ptot_comp_TTYPE , Semicond );
Rs_TTYPE_P(i)=Rs;

%%%%%%%%%% 2LC
[Pc_comp_2LC,Ps_comp_2LC,Ptot_comp_2LC,Ptot_2LC]=TwoLC_freq(U__DC,Vac,Pac,ddeg,f__s);
eff_2LC_P(i)=(Pac-Paux-Ptot_2LC)/Pac*100;

Semicond.Type = 'FF300R12KE4' ;
Semicond.topology = '2LC' ;
[ Rs ]= heat_sink_design ( Ptot_comp_2LC , Semicond );
Rs_2LC_P(i)=Rs;

%%%%%%%%%%% SiC
[Pc_comp_2LC_SiC,Ps_comp_2LC_SiC,Ptot_comp_2LC_SiC,Ptot_2LC_SiC]...
    =TwoLC_freqSiC(U__DC,Vac,Pac,ddeg,f__s);
eff_2LC_SiC_P(i)=(Pac-Paux-Ptot_2LC_SiC)/Pac*100;

Semicond.Type = 'cas300m12bm2' ;
Semicond.topology = '2LC' ;
[ Rs ]= heat_sink_design ( Ptot_comp_2LC_SiC , Semicond );
Rs_2LC_SiC_P(i)=Rs;

%%%%%%%%%%% SiC HP
[Pc_comp_2LC_SiC_HP,Ps_comp_2LC_SiC_HP,Ptot_comp_2LC_SiC_HP,Ptot_2LC_SiC_HP]...
    =TwoLC_freqSiC_HP(U__DC,Vac,Pac,ddeg,f__s);
eff_2LC_SiC_P(i)=(Pac-Paux-Ptot_2LC_SiC_HP)/Pac*100;

Semicond.Type = 'cas300m12bm2' ;
Semicond.topology = '2LC' ;
[ Rs ]= heat_sink_design ( Ptot_comp_2LC_SiC_HP , Semicond );
Rs_2LC_SiC_HP_P(i)=Rs;

end

%%
width=2;  

figure;
hold on
plot(f_range/1000,Rs_NPCC_P,'r','Linewidth',width)
plot(f_range/1000,Rs_TTYPE_P,'k','Linewidth',width)
plot(f_range/1000,Rs_2LC_P,'b','Linewidth',width)
plot(f_range/1000,Rs_2LC_SiC_P,'g','Linewidth',width)
legend('NPCC','TTYPE','2LC-Si','2LC-SiC')
xlabel('Switching frequency [kHz]')
ylabel('Rs [K/W]')
grid on
title('HeatSink Thermal Resistance')
%%
load('fs_Vs_eff.mat')
load('Losses_fs.mat')

figure;
hold on
plot(f_range/1000,eff_NPCC_P,'r','Linewidth',width)
plot(f_range/1000,eff_TTYPE_P,'k','Linewidth',width)
plot(f_range/1000,eff_2LC_P,'b','Linewidth',width)
plot(f_range/1000,eff_2LC_SiC_P,'g','Linewidth',width)
plot(fs/1000,n(1,:)*100,'b--','Linewidth',width)
plot(fs/1000,n(2,:)*100,'m--x','Linewidth',width)
plot(fs/1000,n(3,:)*100,'m--','Linewidth',width)
plot(fs/1000,n(4,:)*100,'g--','Linewidth',width)
legend('NPCC','TTYPE','2LC-Si','2LC-SiC','2LC-2xSi','2LC-Si+4xSiC','2LC-SiC+4xSi','2LC-2xSiC')
xlabel('Switching frequency [kHz]')
ylabel('Efficienccy [%]')
grid on
title('Efficiency')