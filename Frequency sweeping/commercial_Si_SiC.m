%% Design parameters
% ccc
clear all
U__DC=900;
Vac=400;
ddeg=0;
figindx=1;

%% Sweeping frequency - P

f__s1=4000;
f__s2=32000;
P=100000;
f_range=f__s1:1000:f__s2;
Paux=0;

i=1;
ii=0;
for f__s=f_range
ii=ii+1;

Pac=P;


[Pc_comp_2LC,Ps_comp_2LC,Ptot_comp_2LC,Ptot_2LC]=TwoLC_freq(U__DC,Vac,Pac,ddeg,f__s);

eff_2LC_P(ii,i)=(Pac-Paux-Ptot_2LC)/Pac*100;

[Pc_comp_2LC_SiC,Ps_comp_2LC_SiC,Ptot_comp_2LC_SiC,Ptot_2LC_SiC]...
    =TwoLC_freqSiC(U__DC,Vac,Pac,ddeg,f__s);

eff_2LC_SiC_P(ii,i)=(Pac-Paux-Ptot_2LC_SiC)/Pac*100;

end

save('SiC_commercial.mat','eff_2LC_SiC_P')
save('Si_commercial.mat','eff_2LC_P')
%%

figure;
hold on
plot(f_range/1000,eff_2LC_P(:,1))
plot(f_range/1000,eff_2LC_SiC_P(:,1))
% plot(f__s/1000,n*100)
legend('2LC-Si','2LC-SiC')

xlabel('Switching frequency [kHz]')
ylabel('Efficienccy [%]')
grid on