%% Parameters
ccc

m_f=99;

t=0:0.0000001:.02;

I_p=75;
Ip=I_p;

xa = sawtooth(2*pi*50*t*m_f+pi/2,1/2);

I_R=200*0.05*2; % max +/- 5% ripple of the full power RMS current

vp=400/sqrt(3)*sqrt(2);

vdc=800;
m=2*vp/vdc;

%% Reference voltage

va=sin(2*pi*50*t);
vb=sin(2*pi*50*t-pi*2/3);
vc=sin(2*pi*50*t+pi*2/3);


%% Sinusoidal AC currents 
phi = 1 ;
i_a=I_p*sin(2*pi*50*t+phi);
i_b=I_p*sin(2*pi*50*t-pi*2/3+phi);
i_c=I_p*sin(2*pi*50*t+pi*2/3+phi);

%% Ripple AC currents

ia_rip_real=Ip.*sin(2*pi*50*t)+(1/2+m/2.*sin(2*50*pi.*t)).*(1-...
    m.*sin(2.*50.*pi.*t)).*(I_R*0.5).*sawtooth(2.*pi.*50.*t.*m_f,1/2);


ib_rip_real=Ip.*sin(2*pi*50*t-2*pi/3)+(1/2+m/2.*sin(2*50*pi.*t-2*pi/3)).*(1-...
    m.*sin(2.*50.*pi.*t-2*pi/3)).*(I_R*0.5).*sawtooth(2.*pi.*50.*t.*m_f-2*pi/3,1/2);

ic_rip_real=Ip.*sin(2*pi*50*t+2*pi/3)+(1/2+m/2.*sin(2*50*pi.*t+2*pi/3)).*(1-...
    m.*sin(2.*50.*pi.*t+2*pi/3)).*(I_R*0.5).*sawtooth(2.*pi.*50.*t.*m_f+2*pi/3,1/2);

%% Plotting voltage and currents

figure;
plot(t,va,t,vb,t,vc,t,xa)
grid on

figure;
plot(t,i_a,t,i_b,t,i_c)
grid on

figure;
plot(t,ia_rip_real,t,ib_rip_real,t,ic_rip_real)
grid on

%% Gate signals

for i=1:1:size(va,2)
    
    if va(i)>=xa(i)
        gate_a(i)=1;
    else
        gate_a(i)=0;
    end

    if vb(i)>=xa(i)
        gate_b(i)=1;
    else
        gate_b(i)=0;
    end

    if vc(i)>=xa(i)
        gate_c(i)=1;
    else
        gate_c(i)=0;
    end
    
end


%% Leg currents

i_leg_a=i_a.*gate_a;
i_leg_b=i_b.*gate_b;
i_leg_c=i_c.*gate_c;

i_leg_a_RIPPLE=ia_rip_real.*gate_a;
i_leg_b_RIPPLE=ib_rip_real.*gate_b;
i_leg_c_RIPPLE=ic_rip_real.*gate_c;

%% DC Current

i_DC=i_leg_a+i_leg_b+i_leg_c;

i_DC_RIPPLE=i_leg_a_RIPPLE+i_leg_b_RIPPLE+i_leg_c_RIPPLE;

I_DC=mean(i_DC)

I_DC_RMS=rms(i_DC)

I_DC_RIPPLE=mean(i_DC_RIPPLE)

I_DC_RMS_RIPPLE=rms(i_DC_RIPPLE)

%% Plotting

figure;
subplot(3,1,1)
plot(t,i_leg_a)
grid on
subplot(3,1,2)
plot(t,i_leg_b)
grid on
subplot(3,1,3)
plot(t,i_leg_c)
grid on

figure;
subplot(3,1,1)
plot(t,i_leg_a_RIPPLE)
grid on
subplot(3,1,2)
plot(t,i_leg_b_RIPPLE)
grid on
subplot(3,1,3)
plot(t,i_leg_c_RIPPLE)
grid on

figure;
hold on
plot(t,i_DC)
plot(t,i_DC_RIPPLE)


%%
[m,delta]=meshgrid(0.05:0.05:0.95,0:1/180*pi:pi);
I__ac=1;
i_dc_ms=2*sqrt(3)/pi.*m.*(1/4+(cos(delta)).^2);

i_dc_hf_rms=sqrt(2*sqrt(3)/pi.*m.*(1/4+(cos(delta)).^2)-(3.*I__ac.*m.*cos(delta).*(1/4)).^2);

figure;
surf(m,delta,i_dc_hf_rms,'FaceAlpha',0.5,'EdgeColor','none');
xlabel('m')
ylabel('\delta')
zlabel('Ripple')

%%

Fs = 10000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 2000;             % Length of signal
t = (0:L-1)*T;        % Time vector

Y = fft(i_DC,L)/L*2;
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
figure;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%%

t = i_DC(:,1)-min(i_DC(:,1));
I_grid = i_DC(:,2);
Ts = t(2)-t(1);               % Sampling period
Fs = 1/Ts;
N=length(I_grid);             % Number of samples
f=round((0:N-1)*(Fs/N));      % Frequency Range
X =fft(I_grid,N)/N*2;         % FFT of a given current signal  %%% Divide by N/2 to get it normalised // RMS of FFT signal => sqrt(X(1)^2+sum(abs(X(2:N/2)).^2)/2);
X(1)=(X(1)/2);
% Core.I1max = abs(X(2));
Ph_X = unwrap(angle(X(1:N/2))); 
Xa = sort(X(1:N/2),'descend');

