ccc
%%%%  LCL FILTER DESIGN  %%%%
%%%  Two level converter  %%%
%        Marco Stecca       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global vg fs Ig fc  vdcmin vdcmax Pn Inom

%% Set parameters
vg=400/sqrt(3);             % phase voltage
fc=50;                      % grid frequency
Pn=100000;                  % Phase nominal current
Ig=Pn/sqrt(3)/400;          % Phase nominal current
vdcmin=774;                 % DC voltage - minimum 
vdcmax=1004;                % DC voltage - minimum maximum
Inom=300;                   % seminconductor nominal current

Lmin=1e-6;                  % min inductance

Lb=(vg/(2*pi*fc*Ig));       % base inductance
Cb=(Ig/(2*pi*fc*vg));       % base capacitance

%% Variables for the filter design
fs=15000;                                               % switching frequency
factor=0.1;                                            % position of fres respect to fs
qfact=0.05;                                             % max Q absorption
irpp=0.5;                                               % max Peak to peak ripple / p.u.

%% define L/C constraints 

[att,Ltot_min,L1_min,Cmax,Ltot_max,wres]=constraints_def(factor,qfact,irpp);

if Ltot_min>Ltot_max
    
    disp('Impossible Filter Design / Iductance too high')
    
end
%%
% %% LC filter 
%     
% [L1lc,Clc,irppactlc]=filter_lc(L1_min,wres,Cmax,Lmin,Ltot_min);
% 
% %% LCL filter  
% 
% [L1lcl,L2lcl,Clcl,irppactlcl]=filter_lcl(L1_min,wres,Cmax,Lmin,Ltot_min);

%% Find min factor for each switching frequency and LC/LCL filter

fsrange=5000:500:15000;
factorrange=0.05:0.01:0.5;
ii=0;
for fs=fsrange
    ii=ii+1;
    i=0;
    for factor=factorrange
        i=i+1;
        
        [att,Ltot_min,L1_min,Cmax,Ltot_max,wres]=constraints_def(factor,qfact,irpp);
        
        L1_min1(ii,i)=L1_min;
        Ltot_min1(ii,i)=Ltot_min;
        Cmax1(ii,i)=Cmax;
        % Check inductance

        if Ltot_min>Ltot_max

            disp('Impossible Filter Design / Iductance too high')
            
            LC_L(ii+1,i+1)=0.3/1000;
            LC_C(ii+1,i+1)=0.3/1000000;
            LC_Ir(ii+1,i+1)=0.3;
            
            LCL_L1(ii+1,i+1)=0.3/1000;
            LCL_L2(ii+1,i+1)=0.3/1000;
            LCL_C(ii+1,i+1)=0.3/1000000;
            LCL_Ir(ii+1,i+1)=0.3;
            
            LC_feasibility(ii+1,i+1)=2;

            LCL_feasibility(ii+1,i+1)=2; 
            
        elseif L1_min>=Ltot_min %% I am fine with a LC Filter
        
            [L1lc,Clc,irppactlc,lc_feas]=sel_filter_lc(L1_min,wres,Cmax,Lmin,Ltot_max);
            
            LC_feasibility(ii+1,i+1)=max(lc_feas);
            
            LC_L(ii+1,i+1)=min(L1lc);
            LC_C(ii+1,i+1)=max(Clc);
            LC_Ir(ii+1,i+1)=max(irppactlc);

            LCL_feasibility(ii+1,i+1)=-1;
            
            LCL_L1(ii+1,i+1)=0.1/1000;
            LCL_L2(ii+1,i+1)=0.1/1000;
            LCL_C(ii+1,i+1)=0.1/1000000;
            LCL_Ir(ii+1,i+1)=0.1;
            
        else %% I can design a LCL Filter
            
            LC_L(ii+1,i+1)=0.1/1000;
            LC_C(ii+1,i+1)=0.1/1000000;
            LC_Ir(ii+1,i+1)=0.1;
            LC_feasibility(ii+1,i+1)=-1;

           
            [L1lcl,L2lcl,Clcl,irppactlcl,lcl_feas]=sel_filter_lcl(L1_min,wres,Cmax,Lmin,Ltot_min,Ltot_max);
            
            LCL_feasibility(ii+1,i+1)=max(lcl_feas);
            
            LCL_L1(ii+1,i+1)=min(L1lcl);
            LCL_L2(ii+1,i+1)=min(L2lcl);
            LCL_C(ii+1,i+1)=max(Clcl);
            LCL_Ir(ii+1,i+1)=max(irppactlcl); 

        end         
    
    end
    
end
%%
figure;
sgtitle('LC filter - Feasibility')
surf(factorrange,fsrange,LC_feasibility(2:end,2:end));%,'FaceAlpha',0.5,'EdgeColor','none');
xlabel('w_{res}/w_{fs}')
ylabel('fs')


figure;
sgtitle('LCL filter - Feasibility')
surf(factorrange,fsrange,LCL_feasibility(2:end,2:end));%,'FaceAlpha',0.5,'EdgeColor','none');
xlabel('w_{res}/w_{fs}')
ylabel('fs')
zlabel('L_{1,min}')

figure;
surf(factorrange,fsrange,L1_min1);%,'FaceAlpha',0.5,'EdgeColor','none');
xlabel('w_{res}/w_{fs}')
ylabel('fs')
zlabel('L_{1,min}')

figure;
surf(factorrange,fsrange,Ltot_min1);%,'FaceAlpha',0.5,'EdgeColor','none');
xlabel('w_{res}/w_{fs}')
ylabel('fs')
zlabel('L_{tot,min}')

figure;
title('LC filter - L')
surf(factorrange,fsrange,LC_L(2:end,2:end));%,'FaceAlpha',0.5,'EdgeColor','none');
xlabel('w_{res}/w_{fs}')
ylabel('fs [Hz]')
zlabel('L [H]')

figure;
title('LC filter - C')
surf(factorrange,fsrange,LC_C(2:end,2:end));%,'FaceAlpha',0.5,'EdgeColor','none');
xlabel('w_{res}/w_{fs}')
ylabel('fs')
zlabel('C [F]')

figure;
title('LC filter - Ripple')
surf(factorrange,fsrange,LC_Ir(2:end,2:end));%,'FaceAlpha',0.5,'EdgeColor','none');
xlabel('w_{res}/w_{fs}')
ylabel('fs')
zlabel('I_{r,p} [p.u]')

figure;
title('LCL filter - L_{converter}')
surf(factorrange,fsrange,LCL_L1(2:end,2:end));%,'FaceAlpha',0.5,'EdgeColor','none');
xlabel('w_{res}/w_{fs}')
ylabel('fs [Hz]')
zlabel('L_{1} [H]')

figure;
title('LCL filter - L_{grid}')
surf(factorrange,fsrange,LCL_L2(2:end,2:end));%,'FaceAlpha',0.5,'EdgeColor','none');
xlabel('w_{res}/w_{fs}')
ylabel('fs [Hz]')
zlabel('L_{2} [H]')

figure;
title('LCL filter - C')
surf(factorrange,fsrange,LCL_C(2:end,2:end));%,'FaceAlpha',0.5,'EdgeColor','none');
xlabel('w_{res}/w_{fs}')
ylabel('fs')
zlabel('C [F]')

figure;
title('LCL filter - Ripple')
surf(factorrange,fsrange,LCL_Ir(2:end,2:end));%,'FaceAlpha',0.5,'EdgeColor','none');
xlabel('w_{res}/w_{fs}')
ylabel('fs')
zlabel('I_{r,p} [p.u]')

LC_L(1,:)=[0 factorrange];
LC_C(1,:)=[0 factorrange];
LC_Ir(1,:)=[0 factorrange];

LC_L(:,1)=[0 ; transpose(fsrange)];
LC_C(:,1)=[0 ; transpose(fsrange)];
LC_Ir(:,1)=[0 ; transpose(fsrange)];

LCL_L1(1,:)=[0 factorrange];
LCL_L2(1,:)=[0 factorrange];
LCL_C(1,:)=[0 factorrange];
LCL_Ir(1,:)=[0 factorrange];

LCL_L1(:,1)=[0 ; transpose(fsrange)];
LCL_L2(:,1)=[0 ; transpose(fsrange)];
LCL_C(:,1)=[0 ; transpose(fsrange)];
LCL_Ir(:,1)=[0 ; transpose(fsrange)];

%% Per unit transformations
    
Ltot_max_pu=Ltot_max/Lb;

L1pulc=L1lc/Lb;
Cpulc=Clc/Cb;

L1pulcl=L1lcl/Lb;
L2pulcl=L2lcl/Lb;
Cpulcl=Clcl/Cb;

%% Plotting

figure;
sgtitle('2 Level converter possible LCL filter design')
subplot(2,1,1)
hold on
plot(L1pulcl)
plot(L2pulcl)
plot(Cpulcl)
legend('L_{c}','L_{g}','C')
subplot(2,1,2)
plot(irppactlcl)

figure;
sgtitle('2 Level converter possible LC filter design')
subplot(2,1,1)
hold on
plot(L1pulc)
plot(Cpulc)
legend('L','C')
subplot(2,1,2)
plot(irppactlc)