ccc

m_f=200;
fs=50*m_f;
t1=0:0.0000001:.02;

Ip=100;
I_R=200*0.2*2; % max +/- 5% ripple of the full power RMS current


vp=400/sqrt(3)*sqrt(2);

vdc=800;
m=0.65;%2*vp/vdc;

ia=Ip.*sin(2*pi*50*t1);
ib=Ip.*sin(2*pi*50*t1-2*pi/3);
ic=Ip.*sin(2*pi*50*t1+2*pi/3);



% [m2,t2]=meshgrid(0.05:0.05:0.95,0:0.001:0.02);

ripple=(vp.*abs(sin(2*50*pi.*t1))-m.*vp.*(sin(2*50*pi.*t1)).^2);
figure;
hold on
plot(t1,ripple)

grid on

% 
% (1/2+m2/2.*sin(2*50*pi.*t2)).*(m2.*sin(2.*50.*pi.*t2));

% figure;
% surf(m2,t2,ripple,'FaceAlpha',0.5,'EdgeColor','none');
% % colormap(1-gray)
% xlabel('m')
% ylabel('t')
% zlabel('Ripple')

%% Phase A
ia_rip_real=Ip.*sin(2*pi*50*t1)+((abs(sin(2*50*pi.*t1))-m.*(sin(2*50*pi.*t1)).^2))*(I_R.*m*0.5)...
    .*sawtooth(2.*pi.*50.*t1.*m_f,1/2);


ia_rip_up=Ip.*sin(2*pi*50*t1)+(abs(sin(2*50*pi.*t1))-m.*(sin(2*50*pi.*t1)).^2)...
    .*(I_R*0.5);

ia_rip_down=Ip.*sin(2*pi*50*t1)-(abs(sin(2*50*pi.*t1))-m.*(sin(2*50*pi.*t1)).^2)...
    .*(I_R*0.5);

%% Phase B

ib_rip_real=Ip.*sin(2*pi*50*t1-2*pi/3)+(1/2-m/2.*abs(sin(2*50*pi.*t1-2*pi/3))).*(...
    m.*sin(2.*50.*pi.*t1-2*pi/3)).*(I_R*0.5).*sawtooth(2.*pi.*50.*t1.*m_f-2*pi/3,1/2);

ib_rip_up=Ip.*sin(2*pi*50*t1-2*pi/3)+(1/2-m/2.*abs(sin(2*50*pi.*t1-2*pi/3))).*(...
    m.*sin(2.*50.*pi.*t1-2*pi/3)).*(I_R*0.5);

ib_rip_down=Ip.*sin(2*pi*50*t1-2*pi/3)-(1/2-m/2.*abs(sin(2*50*pi.*t1-2*pi/3))).*(...
    m.*sin(2.*50.*pi.*t1-2*pi/3)).*(I_R*0.5);

%% Phase C

ic_rip_real=Ip.*sin(2*pi*50*t1+2*pi/3)+(1/2-m/2.*abs(sin(2*50*pi.*t1+2*pi/3))).*(...
    m.*sin(2.*50.*pi.*t1+2*pi/3)).*(I_R*0.5).*sawtooth(2.*pi.*50.*t1.*m_f+2*pi/3,1/2);

ic_rip_up=Ip.*sin(2*pi*50*t1+2*pi/3)+(1/2-m/2.*abs(sin(2*50*pi.*t1+2*pi/3))).*(...
    m.*sin(2.*50.*pi.*t1+2*pi/3)).*(I_R*0.5);

ic_rip_down=Ip.*sin(2*pi*50*t1+2*pi/3)-(1/2-m/2.*abs(sin(2*50*pi.*t1+2*pi/3))).*(...
    m.*sin(2.*50.*pi.*t1+2*pi/3)).*(I_R*0.5);


%% Plotting

figure;
hold on
plot(t1,ia,t1,ib,t1,ic)

plot(t1,ia_rip_real,t1,ib_rip_real,t1,ic_rip_real)

plot(t1,ia_rip_real,t1,ia_rip_up,t1,ia_rip_down)

plot(t1,ib_rip_real,t1,ib_rip_up,t1,ib_rip_down)

plot(t1,ic_rip_real,t1,ic_rip_up,t1,ic_rip_down)
grid on

%%
figure;
hold on
plot(t1,ia)

plot(t1,ia_rip_real,t1,ia_rip_up,t1,ia_rip_down)

grid on


figure;
hold on
plot(t1,ia_rip_real)

grid on

figure;
hold on
plot(t1,((abs(sin(2*50*pi.*t1))-m.*(sin(2*50*pi.*t1)).^2))*(I_R.*m*0.5)...
    .*sawtooth(2.*pi.*50.*t1.*m_f,1/2))
grid on