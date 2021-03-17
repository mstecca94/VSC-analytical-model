ccc
%% Design parameters
U__DC=900;
Vac=400;
ddeg=0;
figindx=1;
f__s=10000; % f at 10 kHz

% Pv=[0  0  0  0  0  0  0 2.7 8.8 25.1 41.6 55.6 65.7 71.9 74 71.3 65.5 55.6 41.8 26.7 11.7 03.9 0.3 0  0];
%   % 1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21  22 23
% lpr=-[30 27 25 22 20 20 21 30 40 42 46 50 50 55 60 60 55 60 65 85 90 85 75 55 30];
      % 1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21  22 23    
q_bat=[30 63 69 64 50 25  0  0  0  0   0  0  0  0  0  0  0  0  0  0  0   0 10 20];

p_bat=[ 0  0  0  0  0  0  0  0  0  0  30  60 80  85  60 25 0 0 -20 -45 -15 0 0 0];



x=1/60/60:1/60/60:24/1;

xxx=1:24;

yyP_bat=pchip(xxx,p_bat,x);
yyQ_bat=pchip(xxx,q_bat,x);
yyS_bat=sqrt(yyP_bat^2+yyQ_bat^2);


figure(figindx);
plot(x,yyP_bat,x,yyQ_bat,x,yyS_bat)
figindx=figindx+1;
%% Battery parameters
load('VOC_BOSTON.mat')

V_max_cell=max(VOC_SoC(:,1));
U__DC_max=1004; % parameters of Samusung M3-R089
n_cells=round(U__DC_max/4.2); %% 4.2 max voltage of boston cells

VOC_SoC(:,1)=VOC_SoC(:,1)*n_cells;

xq=0:0.00001:1;
vq = interp1(VOC_SoC(:,2),VOC_SoC(:,1),xq);
OCV_SoC=[round(xq*100,3) ; vq];

% 
% figure(figindx);
% plot(OCV_SoC(1,:),OCV_SoC(2,:))
% figindx=figindx+1;

Ebattery=89000; %% 89 kWh



P_nom_batt=100000; % also battery 100 kW nominal power

% eta_battery=1-0.02*P_batt/P_nom_batt; % losses 2% at 1C-


SoC_0=round((Ebattery*1.5/2)/Ebattery*100,3);

U__DC_0=OCV_SoC(2,find(OCV_SoC(1,:)==SoC_0));
%% Sweeping P for frequency regulation


time_span=24*60*60; %%% this is 1 day --> 24h*60min*60s=86400 

Psp=Psp2(1:time_span)*1000000; %% in W


eff_ANPCC=zeros(time_span,1);
eff_NPCC=zeros(time_span,1);
eff_TTYPE=zeros(time_span,1);
eff_2LC=zeros(time_span,1);

U_dc_ANPCC=zeros(time_span,1);
U_dc_NPCC=zeros(time_span,1);
U_dc_TTYPE=zeros(time_span,1);
U_dc_2LC=zeros(time_span,1);

U_dc_ANPCC(1)=U__DC_0;
U_dc_NPCC(1)=U__DC_0;
U_dc_TTYPE(1)=U__DC_0;
U_dc_2LC(1)=U__DC_0;

SoC_ANPCC=zeros(time_span,1);
SoC_NPCC=zeros(time_span,1);
SoC_TTYPE=zeros(time_span,1);
SoC_2LC=zeros(time_span,1);

SoC_ANPCC(1)=SoC_0;
SoC_NPCC(1)=SoC_0;
SoC_TTYPE(1)=SoC_0;
SoC_2LC(1)=SoC_0;

Qsp=zeros(time_span,1);

%% ANPCC
for i=1:1:size(Psp,1) 
    
Pneg=0;

% initialize time values

Pac=Psp(i);

Qac=Qsp(i);

Sac=sqrt(Pac^2+Qac^2);

Ssp(i)=Sac;

ddeg=atand(Qac/Pac);

if isnan(ddeg) || (SoC_ANPCC(i)<10 && Pac>0) || (SoC_ANPCC(i)>90 && Pac<0)

    SoC_ANPCC(i+1)=SoC_ANPCC(i);

    U_dc_ANPCC(i+1)=U_dc_ANPCC(i);
    
    Pb_ANPCC(i)=0;

    Losses_ANPCC(i)=0;

    eff_ANPCC(i)=0;
    
    Pbatt_ANPCC(i)=0;
    
    eta_batt_ANPCC(i)=0;

continue

end

if Pac<0
    Pac=-Pac;
    ddeg=ddeg+180;
    
    Pneg=1; % rectifier
    
end

%% find inverter losses

% ANPCC

U_DC_ANPCC=U_dc_ANPCC(i);

[Pc_comp_ANPCC,Ps_comp_ANPCC,Ptot_comp_ANPCC,Ptot_ANPCC]=anpcc_freq(U_DC_ANPCC,Vac,Sac,ddeg,f__s);

Losses_ANPCC(i)=Ptot_ANPCC;

eff_ANPCC(i)=(Pac-Ptot_ANPCC)/Pac*100;

%% battery 

% power between battery and converter

if Pneg==1 % rectifier
    
    Pbatt_ANPCC(i)=Pac*eff_ANPCC(i)/100;
    
elseif Pneg==0 % inverter
    
    Pbatt_ANPCC(i)=Pac/eff_ANPCC(i)/100;

end

% battery efficicency

eta_batt_ANPCC(i)=1-0.02*abs(Pbatt_ANPCC(i))/P_nom_batt; % losses 2% at 1C-

% power from the battery and new State of charge

if Pneg==1 % rectifier
    
    Pb_ANPCC(i)=Pac*eta_batt_ANPCC(i);
    
    SoC_ANPCC(i+1)=SoC_ANPCC(i)+round((Pb_ANPCC(i)/3600)/Ebattery*100,3);
    
elseif Pneg==0 % inverter
    
    Pb_ANPCC(i)=Pac/eta_batt_ANPCC(i);

    SoC_ANPCC(i+1)=round(SoC_ANPCC(i)-(Pb_ANPCC(i)/3600)/Ebattery*100,3);
    
end

U_dc_ANPCC(i+1)=OCV_SoC(2,round(SoC_ANPCC(i+1)*1000-1)); 

end

%% NPCC
for i=1:1:size(Psp,1)
    
Pneg=0;

% initialize time values

Pac=Psp(i);

Qac=Qsp(i);

Sac=sqrt(Pac^2+Qac^2);

Ssp(i)=Sac;

ddeg=atand(Qac/Pac);

if isnan(ddeg) || (SoC_NPCC(i)<10 && Pac>0) || (SoC_NPCC(i)>90 && Pac<0)

    SoC_NPCC(i+1)=SoC_NPCC(i);
    
    U_dc_NPCC(i+1)=U_dc_NPCC(i);
    
    Pb_NPCC(i)=0;
    
    Losses_NPCC(i)=0;

    eff_NPCC(i)=0;

    Pbatt_NPCC(i)=0;

    eta_batt_NPCC(i)=0;

continue

end

if Pac<0
    Pac=-Pac;
    ddeg=ddeg+180;
    
    Pneg=1; % rectifier
    
end

%% find inverter losses

% NPCC

U_DC_NPCC=U_dc_NPCC(i);

[Pc_comp_NPCC,Ps_comp_NPCC,Ptot_comp_NPCC,Ptot_NPCC]=npcc_freq(U_DC_NPCC,Vac,Sac,ddeg,f__s);

Losses_NPCC(i)=Ptot_NPCC;

eff_NPCC(i)=(Pac-Ptot_NPCC)/Pac*100;

%% battery 

% power between battery and converter

if Pneg==1 % rectifier
    
    Pbatt_NPCC(i)=Pac*eff_NPCC(i)/100;
    
elseif Pneg==0 % inverter
    
    Pbatt_NPCC(i)=Pac/eff_NPCC(i)/100;
    
end

% battery efficicency

eta_batt_NPCC(i)=1-0.02*abs(Pbatt_NPCC(i))/P_nom_batt; % losses 2% at 1C-

% power from the battery and new State of charge

if Pneg==1 % rectifier
    
    Pb_NPCC(i)=Pac*eta_batt_NPCC(i);
    
    SoC_NPCC(i+1)=SoC_NPCC(i)+round((Pb_NPCC(i)/3600)/Ebattery*100,3);
elseif Pneg==0 % inverter
    
    Pb_NPCC(i)=Pac/eta_batt_NPCC(i);
    
    SoC_NPCC(i+1)=round(SoC_NPCC(i)-(Pb_NPCC(i)/3600)/Ebattery*100,3);
    
end

U_dc_NPCC(i+1)=OCV_SoC(2,round(SoC_NPCC(i+1)*1000-1)); 

end

%% T-TYPE

for i=1:1:size(Psp,1)
    
Pneg=0;

% initialize time values

Pac=Psp(i);

Qac=Qsp(i);

Sac=sqrt(Pac^2+Qac^2);

Ssp(i)=Sac;    

ddeg=atand(Qac/Pac);

if isnan(ddeg) || (SoC_TTYPE(i)<10 && Pac>0) || (SoC_TTYPE(i)>90 && Pac<0)

    SoC_TTYPE(i+1)=SoC_TTYPE(i);
    
    U_dc_TTYPE(i+1)=U_dc_TTYPE(i);
    
    Pb_TTYPE(i)=0;
    
    Losses_TTYPE(i)=0;
    
    eff_TTYPE(i)=0;
    
    Pbatt_TTYPE(i)=0;
    
    eta_batt_TTYPE(i)=0;

continue

end

if Pac<0
    Pac=-Pac;
    ddeg=ddeg+180;
    
    Pneg=1; % rectifier
    
end

%% find inverter losses / ttype

U_DC_TTYPE=U_dc_TTYPE(i);

[Pc_comp_TTYPE,Ps_comp_TTYPE,Ptot_comp_TTYPE,Ptot_TTYPE]=ttype2_freq(U_DC_TTYPE,Vac,Sac,ddeg,f__s);

Losses_TTYPE(i)=Ptot_TTYPE;

eff_TTYPE(i)=(Pac-Ptot_TTYPE)/Pac*100;


%% battery 

% power between battery and converter

if Pneg==1 % rectifier
    
    Pbatt_TTYPE(i)=Pac*eff_TTYPE(i)/100;

elseif Pneg==0 % inverter
    
    Pbatt_TTYPE(i)=Pac/eff_TTYPE(i)/100;
    
end

% battery efficicency

eta_batt_TTYPE(i)=1-0.02*abs(Pbatt_TTYPE(i))/P_nom_batt; % losses 2% at 1C-

% power from the battery and new State of charge

if Pneg==1 % rectifier
    
    Pb_TTYPE(i)=Pac*eta_batt_TTYPE(i);
    
    SoC_TTYPE(i+1)=SoC_TTYPE(i)+round((Pb_TTYPE(i)/3600)/Ebattery*100,3);
    
elseif Pneg==0 % inverter
    
    Pb_TTYPE(i)=Pac/eta_batt_TTYPE(i);
    
    SoC_TTYPE(i+1)=round(SoC_TTYPE(i)-(Pb_TTYPE(i)/3600)/Ebattery*100,3);
    
end

U_dc_TTYPE(i+1)=OCV_SoC(2,round(SoC_TTYPE(i+1)*1000-1)); 

end

%% 2LC
for i=1:1:size(Psp,1)
    
Pneg=0;

% initialize time values

Pac=Psp(i);

Qac=Qsp(i);

Sac=sqrt(Pac^2+Qac^2);

Ssp(i)=Sac;

ddeg=atand(Qac/Pac);

if isnan(ddeg) || (SoC_2LC(i)<10 && Pac>0) || (SoC_2LC(i)>90 && Pac<0)

    SoC_2LC(i+1)=SoC_2LC(i); 

    U_dc_2LC(i+1)=U_dc_2LC(i);     

    Pb_2LC(i)=0;  

    Losses_2LC(i)=0;

    eff_2LC(i)=0;

    Pbatt_2LC(i)=0;   

    eta_batt_2LC(i)=0;

continue

end

if Pac<0
    Pac=-Pac;
    ddeg=ddeg+180;
    
    Pneg=1; % rectifier
    
end

%% find inverter losses / 2LC

U_DC_2LC=U_dc_2LC(i);

[Pc_comp_2LC,Ps_comp_2LC,Ptot_comp_2LC,Ptot_2LC]=TwoLC_freq(U_DC_2LC,Vac,Sac,ddeg,f__s);

Losses_2LC(i)=Ptot_2LC;

eff_2LC(i)=(Pac-Ptot_2LC)/Pac*100;


%% battery 

% power between battery and converter

if Pneg==1 % rectifier
    
    Pbatt_2LC(i)=Pac*eff_2LC(i)/100;

elseif Pneg==0 % inverter

    Pbatt_2LC(i)=Pac/eff_2LC(i)/100;   
    
end

% battery efficicency

eta_batt_2LC(i)=1-0.02*(abs(Pbatt_2LC(i))/P_nom_batt); % losses 2% at 1C-

% power from the battery and new State of charge

if Pneg==1 % rectifier
    
    Pb_2LC(i)=Pac*eta_batt_2LC(i);

    SoC_2LC(i+1)=SoC_2LC(i)+round((Pb_2LC(i)/3600)/Ebattery*100,3);
    
elseif Pneg==0 % inverter
    
    Pb_2LC(i)=Pac/eta_batt_2LC(i);   

    SoC_2LC(i+1)=round(SoC_2LC(i)-(Pb_2LC(i)/3600)/Ebattery*100,3);

end


U_dc_2LC(i+1)=OCV_SoC(2,round(SoC_2LC(i+1)*1000-1)); 

end

%% matrix manipulation

% full system efficiency

eta_batt_2LC=eta_batt_2LC*100;
eta_batt_TTYPE=eta_batt_TTYPE*100;
eta_batt_NPCC=eta_batt_NPCC*100;
eta_batt_ANPCC=eta_batt_ANPCC*100;


eff_2LC_system=eff_2LC.*transpose(eta_batt_2LC)/100;
eff_TTYPE_system=eff_TTYPE.*transpose(eta_batt_TTYPE)/100;
eff_NPCC_system=eff_NPCC.*transpose(eta_batt_NPCC)/100;
eff_ANPCC_system=eff_ANPCC.*transpose(eta_batt_ANPCC)/100;

eff_2LC_system=[transpose(1:1:time_span) eff_2LC_system ];
eff_TTYPE_system=[transpose(1:1:time_span) eff_TTYPE_system ];
eff_NPCC_system=[transpose(1:1:time_span) eff_NPCC_system];
eff_ANPCC_system=[transpose(1:1:time_span) eff_ANPCC_system ];

eff_2LC_system( ~any(eff_2LC_system(:,2),2), : ) = [];  %rows
eff_TTYPE_system( ~any(eff_TTYPE_system(:,2),2), : ) = [];  %rows
eff_NPCC_system( ~any(eff_NPCC_system(:,2),2), : ) = [];  %rows
eff_ANPCC_system( ~any(eff_ANPCC_system(:,2),2), : ) = [];  %rows

% mean full y=system efficiency efficiency

eta_syst_2LC_mean=mean(eff_2LC_system(:,2));
eta_syst_TTYPE_mean=mean(eff_TTYPE_system(:,2));
eta_syst_NPCC_mean=mean(eff_NPCC_system(:,2));
eta_syst_ANPCC_mean=mean(eff_ANPCC_system(:,2));


Energy_cycled=sum(abs(Psp));

timerange=0:1:time_span;
timerange3=1:1:time_span;

timerange2=transpose(1:1:time_span);

% battery efficiency

eta_batt_2LC=[(1:1:time_span) ;eta_batt_2LC ];
eta_batt_TTYPE=[(1:1:time_span) ;eta_batt_TTYPE ];
eta_batt_NPCC=[(1:1:time_span) ;eta_batt_NPCC];
eta_batt_ANPCC=[(1:1:time_span); eta_batt_ANPCC ];

eta_batt_2LC(:, ~any(eta_batt_2LC(2,:),1) ) = [];  %columns
eta_batt_TTYPE(:, ~any(eta_batt_TTYPE(2,:),1) ) = [];  %columns
eta_batt_NPCC(:, ~any(eta_batt_NPCC(2,:),1) ) = [];  %columns
eta_batt_ANPCC(:, ~any(eta_batt_ANPCC(2,:),1) ) = [];  %columns

% mean battery efficiency

eta_batt_2LC_mean=mean(eta_batt_2LC(2,:));
eta_batt_TTYPE_mean=mean(eta_batt_TTYPE(2,:));
eta_batt_NPCC_mean=mean(eta_batt_NPCC(2,:));
eta_batt_ANPCC_mean=mean(eta_batt_ANPCC(2,:));

% converter efficiency

eff_2LC=[transpose(1:1:time_span) eff_2LC ];
eff_TTYPE=[transpose(1:1:time_span) eff_TTYPE ];
eff_NPCC=[transpose(1:1:time_span) eff_NPCC];
eff_ANPCC=[transpose(1:1:time_span) eff_ANPCC ];

eff_2LC( ~any(eff_2LC(:,2),2), : ) = [];  %rows
eff_TTYPE( ~any(eff_TTYPE(:,2),2), : ) = [];  %rows
eff_NPCC( ~any(eff_NPCC(:,2),2), : ) = [];  %rows
eff_ANPCC( ~any(eff_ANPCC(:,2),2), : ) = [];  %rows

% mean converter efficiency

eff_2LC_mean=mean(eff_2LC(:,2));
eff_TTYPE_mean=mean(eff_TTYPE(:,2));
eff_NPCC_mean=mean(eff_NPCC(:,2));
eff_ANPCC_mean=mean(eff_ANPCC(:,2));

%%% efficiency

eta_2lC=[eff_2LC_mean  eta_batt_2LC_mean eta_syst_2LC_mean];
eta_NPCC=[eff_NPCC_mean  eta_batt_NPCC_mean eta_syst_NPCC_mean];
eta_ANPCC=[eff_ANPCC_mean  eta_batt_ANPCC_mean eta_syst_ANPCC_mean];
eta_TTYPE=[eff_TTYPE_mean  eta_batt_TTYPE_mean eta_syst_TTYPE_mean];

eta=[eta_2lC ;eta_TTYPE ;  eta_NPCC ; eta_ANPCC];

%% Plotting
figure(figindx)
plot(timerange,SoC_2LC,timerange,SoC_TTYPE,timerange,SoC_NPCC,timerange,SoC_ANPCC)
title('State of Charge')
legend('ANPCC','NPCC','TTYPE','2LC')
figindx=figindx+1;

figure(figindx)
plot(timerange,U_dc_2LC,timerange,U_dc_TTYPE,timerange,U_dc_NPCC,timerange,U_dc_ANPCC)
title('U_{DC}')
legend('ANPCC','NPCC','TTYPE','2LC')
figindx=figindx+1;


figure(figindx)
hold on
scatter(eff_ANPCC(:,1),eff_ANPCC(:,2),1,'filled')
scatter(eff_NPCC(:,1),eff_NPCC(:,2),1,'filled')
scatter(eff_TTYPE(:,1),eff_TTYPE(:,2),1,'filled')
scatter(eff_2LC(:,1),eff_2LC(:,2),1,'filled')
title('Converter effieciency')
legend('ANPCC','NPCC','TTYPE','2LC')
figindx=figindx+1;

figure(figindx)
hold on
scatter(eta_batt_ANPCC(1,:),eta_batt_ANPCC(2,:),1,'filled')
scatter(eta_batt_NPCC(1,:),eta_batt_NPCC(2,:),1,'filled')
scatter(eta_batt_TTYPE(1,:),eta_batt_TTYPE(2,:),1,'filled')
scatter(eta_batt_2LC(1,:),eta_batt_2LC(2,:),1,'filled')
title('Battery effieciency')
legend('ANPCC','NPCC','TTYPE','2LC')
figindx=figindx+1;

figure(figindx)
hold on
scatter(eff_ANPCC_system(:,1),eff_ANPCC_system(:,2),1,'filled')
scatter(eff_NPCC_system(:,1),eff_NPCC_system(:,2),1,'filled')
scatter(eff_TTYPE_system(:,1),eff_TTYPE_system(:,2),1,'filled')
scatter(eff_2LC_system(:,1),eff_2LC_system(:,2),1,'filled')
title('Full system effieciency')
legend('ANPCC','NPCC','TTYPE','2LC')
figindx=figindx+1;


figure(figindx)
plot(timerange3,Pbatt_2LC,timerange3,Pbatt_TTYPE,timerange3,Pbatt_NPCC,timerange3,Pbatt_ANPCC)
title('Power from the battery')
legend('ANPCC','NPCC','TTYPE','2LC')
figindx=figindx+1;


figure(figindx)
bar(eta)
ylim([95 100])
title('Mean efficiency')
xticklabels({'2LC','T-Type','NPCC','ANPCC'})
legend('Converter','Battery','System')
figindx=figindx+1;

