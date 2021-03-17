ccc
%% Design parameters
U__DC=900;
Vac=400;
Pac=100000; %% 100 kW
ddeg=0;
figindx=1;
f__s=10000;

[Pc_comp_ANPCC,Ps_comp_ANPCC,Ptot_comp_ANPCC,Ptot_ANPCC]=anpcc_comp(U__DC,Vac,Pac,ddeg,f__s);



[Pc_comp_NPCC,Ps_comp_NPCC,Ptot_comp_NPCC,Ptot_NPCC]=npcc_comp(U__DC,Vac,Pac,ddeg,f__s);


[Pc_comp_TTYPE,Ps_comp_TTYPE,Ptot_comp_TTYPE,Ptot_TTYPE]=ttype2_comp(U__DC,Vac,Pac,ddeg,f__s);



[Pc_comp,Ps_comp,Ptot_comp,Ptot]=TwoLC_comp(U__DC,Vac,Pac,ddeg,f__s);

