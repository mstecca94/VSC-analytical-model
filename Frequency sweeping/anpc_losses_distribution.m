ccc
%% Design parameters
U__DC=900;
Vac=400;
Pac=100000; %% 100 kW
ddeg=0;
figindx=1;

f__s=12000; % f at 20 kHz


t__a=0.5;
U__s=1;
L__s=1;


[Pc_comp,Ps_comp,~,~]=anpcc_switching_patterns(U__DC,Vac,Pac,ddeg,f__s,t__a,U__s,L__s)


figure(figindx)
subplot(3,4,1)
hold on
bar(transpose([Pc_comp; Ps_comp]),'stacked')
xticks(1:12)
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','T5','D5','T6','D6'})
title(sprintf('t_{a}=%.1d U_{s}=%d L_{s}=%d',t__a,U__s,L__s))
grid on
ylim([0 125])
t__a=0.5;
U__s=0;
L__s=1;


[Pc_comp,Ps_comp,~,~]=anpcc_switching_patterns(U__DC,Vac,Pac,ddeg,f__s,t__a,U__s,L__s)

subplot(3,4,2)
hold on
bar(transpose([Pc_comp; Ps_comp]),'stacked')
xticks(1:12)
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','T5','D5','T6','D6'})
title(sprintf('t_{a}=%.1d U_{s}=%d L_{s}=%d',t__a,U__s,L__s))
grid on
ylim([0 125])

t__a=0.5;
U__s=1;
L__s=0;


[Pc_comp,Ps_comp,~,~]=anpcc_switching_patterns(U__DC,Vac,Pac,ddeg,f__s,t__a,U__s,L__s)

subplot(3,4,3)
hold on
bar(transpose([Pc_comp; Ps_comp]),'stacked')
xticks(1:12)
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','T5','D5','T6','D6'})
title(sprintf('t_{a}=%.1d U_{s}=%d L_{s}=%d',t__a,U__s,L__s))
grid on
ylim([0 125])

t__a=0.5;
U__s=0;
L__s=0;


[Pc_comp,Ps_comp,~,~]=anpcc_switching_patterns(U__DC,Vac,Pac,ddeg,f__s,t__a,U__s,L__s)

subplot(3,4,4)
hold on
bar(transpose([Pc_comp; Ps_comp]),'stacked')
xticks(1:12)
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','T5','D5','T6','D6'})
title(sprintf('t_{a}=%.1d U_{s}=%d L_{s}=%d',t__a,U__s,L__s))
grid on
ylim([0 125])

t__a=0;
U__s=1;
L__s=1;


[Pc_comp,Ps_comp,~,~]=anpcc_switching_patterns(U__DC,Vac,Pac,ddeg,f__s,t__a,U__s,L__s)


subplot(3,4,5)
hold on
bar(transpose([Pc_comp; Ps_comp]),'stacked')
xticks(1:12)
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','T5','D5','T6','D6'})
title(sprintf('t_{a}=%.1d U_{s}=%d L_{s}=%d',t__a,U__s,L__s))
grid on
ylim([0 125])

t__a=0;
U__s=0;
L__s=1;


[Pc_comp,Ps_comp,~,~]=anpcc_switching_patterns(U__DC,Vac,Pac,ddeg,f__s,t__a,U__s,L__s)

subplot(3,4,6)
hold on
bar(transpose([Pc_comp; Ps_comp]),'stacked')
xticks(1:12)
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','T5','D5','T6','D6'})
title(sprintf('t_{a}=%.1d U_{s}=%d L_{s}=%d',t__a,U__s,L__s))
grid on
ylim([0 125])

t__a=0;
U__s=1;
L__s=0;


[Pc_comp,Ps_comp,~,~]=anpcc_switching_patterns(U__DC,Vac,Pac,ddeg,f__s,t__a,U__s,L__s)

subplot(3,4,7)
hold on
bar(transpose([Pc_comp; Ps_comp]),'stacked')
xticks(1:12)
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','T5','D5','T6','D6'})
title(sprintf('t_{a}=%.1d U_{s}=%d L_{s}=%d',t__a,U__s,L__s))
grid on
ylim([0 125])

t__a=0;
U__s=0;
L__s=0;


[Pc_comp,Ps_comp,~,~]=anpcc_switching_patterns(U__DC,Vac,Pac,ddeg,f__s,t__a,U__s,L__s)

subplot(3,4,8)
hold on
bar(transpose([Pc_comp; Ps_comp]),'stacked')
xticks(1:12)
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','T5','D5','T6','D6'})
title(sprintf('t_{a}=%.1d U_{s}=%d L_{s}=%d',t__a,U__s,L__s))
grid on
ylim([0 125])

t__a=1;
U__s=1;
L__s=1;


[Pc_comp,Ps_comp,~,~]=anpcc_switching_patterns(U__DC,Vac,Pac,ddeg,f__s,t__a,U__s,L__s)


subplot(3,4,9)
hold on
bar(transpose([Pc_comp; Ps_comp]),'stacked')
xticks(1:12)
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','T5','D5','T6','D6'})
title(sprintf('t_{a}=%.1d U_{s}=%d L_{s}=%d',t__a,U__s,L__s))
grid on
ylim([0 125])

t__a=1;
U__s=0;
L__s=1;


[Pc_comp,Ps_comp,~,~]=anpcc_switching_patterns(U__DC,Vac,Pac,ddeg,f__s,t__a,U__s,L__s)

subplot(3,4,10)
hold on
bar(transpose([Pc_comp; Ps_comp]),'stacked')
xticks(1:12)
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','T5','D5','T6','D6'})
title(sprintf('t_{a}=%.1d U_{s}=%d L_{s}=%d',t__a,U__s,L__s))
grid on
ylim([0 125])

t__a=1;
U__s=1;
L__s=0;


[Pc_comp,Ps_comp,~,~]=anpcc_switching_patterns(U__DC,Vac,Pac,ddeg,f__s,t__a,U__s,L__s)

subplot(3,4,11)
hold on
bar(transpose([Pc_comp; Ps_comp]),'stacked')
xticks(1:12)
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','T5','D5','T6','D6'})
title(sprintf('t_{a}=%.1d U_{s}=%d L_{s}=%d',t__a,U__s,L__s))
grid on
ylim([0 125])

t__a=1;
U__s=0;
L__s=0;


[Pc_comp,Ps_comp,~,~]=anpcc_switching_patterns(U__DC,Vac,Pac,ddeg,f__s,t__a,U__s,L__s)

subplot(3,4,12)
hold on
bar(transpose([Pc_comp; Ps_comp]),'stacked')
xticks(1:12)
xticklabels({'T1','D1','T2','D2','T3','D3','T4','D4','T5','D5','T6','D6'})
title(sprintf('t_{a}=%.1d U_{s}=%d L_{s}=%d',t__a,U__s,L__s))
grid on
ylim([0 125])
