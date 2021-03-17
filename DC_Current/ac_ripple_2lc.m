ccc

m_f=39;
fs=50*m_f;
t1=0:0.0000001:.02;

Ip=10;
I_R=200*0.1*2; % max +/- 5% ripple of the full power RMS current


vp=400/sqrt(3)*sqrt(2);

vdc=800;
m=1;2*vp/vdc;

ia=Ip.*sin(2*pi*50*t1);
ib=Ip.*sin(2*pi*50*t1-2*pi/3);
ic=Ip.*sin(2*pi*50*t1+2*pi/3);

vp2=400/sqrt(3)*sqrt(2);

[vdc2,t2]=meshgrid(655:5:1100,0:0.001:0.02);

ripple=(vdc2/2+vp2.*sin(2*50*pi.*t2)).*(1/2-vp2./vdc2.*sin(2.*50.*pi.*t2));
ripple3l=(vdc2/2.*(abs(2*vp2./vdc2.*sin(2.*50.*pi.*t2))+2*vp2./vdc2.*vp2.*(sin(2*50*pi.*t2)).^2));


figure;
surf(vdc2,t2,ripple,'FaceAlpha',0.5,'EdgeColor','none');
% colormap(1-gray)
xlabel('V_{DC}')
ylabel('t')
zlabel('Ripple')

figure;
surf(vdc2,t2,ripple3l,'FaceAlpha',0.5,'EdgeColor','none');
% colormap(1-gray)
xlabel('V_{DC}')
ylabel('t')
zlabel('Ripple 3L')

ta2=0:0.0000001:.02/6;
ta1=0.02/6:0.0000001:.02/3;
figure;
hold on
plot(ta1,(5/6-m/2.*sin(2*50*pi.*ta1)).*(1/2+...
    m/2.*sin(2.*50.*pi.*ta1)))
plot(ta2,(1/3-m/2.*sin(2*50*pi.*ta2)).*(1/2+...
    m/2.*sin(2.*50.*pi.*ta2)))%.*(I_R*0.5).*sawtooth(2.*pi.*50.*t1.*m_f,1/2))
% plot(t1,(1/2-m/2.*sin(2*50*pi.*t1)).*(m.*sin(2.*50.*pi.*t1)))
% plot(t1,(-1/2-m/2.*sin(2*50*pi.*t1)).*(m.*sin(2.*50.*pi.*t1+pi)))
hold off
%% Phase A
ia_rip_real=Ip.*sin(2*pi*50*t1)+(1/2+m/2.*sin(2*50*pi.*t1)).*(1/2-...
    m/2.*sin(2.*50.*pi.*t1)).*(I_R*0.5).*sawtooth(2.*pi.*50.*t1.*m_f,1/2);

ia_rip_up=Ip.*sin(2*pi*50*t1)+(1/2+m/2.*sin(2*50*pi.*t1)).*(1/2-...
    m/2.*sin(2.*50.*pi.*t1)).*(I_R*0.5);

ia_rip_down=Ip.*sin(2*pi*50*t1)-(1/2+m/2.*sin(2*50*pi.*t1)).*(1/2-...
    m/2.*sin(2.*50.*pi.*t1)).*(I_R*0.5);

%% Phase B

ib_rip_real=Ip.*sin(2*pi*50*t1-2*pi/3)+(1/2+m/2.*sin(2*50*pi.*t1-2*pi/3)).*(1/2-...
    m/2.*sin(2.*50.*pi.*t1-2*pi/3)).*(I_R*0.5).*sawtooth(2.*pi.*50.*t1.*m_f-2*pi/3,1/2);

ib_rip_up=Ip.*sin(2*pi*50*t1-2*pi/3)+(1/2+m/2.*sin(2*50*pi.*t1-2*pi/3)).*(1/2-...
    m/2.*sin(2.*50.*pi.*t1-2*pi/3)).*(I_R*0.5);

ib_rip_down=Ip.*sin(2*pi*50*t1-2*pi/3)-(1/2+m/2.*sin(2*50*pi.*t1-2*pi/3)).*(1/2-...
    m/2.*sin(2.*50.*pi.*t1-2*pi/3)).*(I_R*0.5);

%% Phase C

ic_rip_real=Ip.*sin(2*pi*50*t1+2*pi/3)+(1/2+m/2.*sin(2*50*pi.*t1+2*pi/3)).*(1/2-...
    m/2.*sin(2.*50.*pi.*t1+2*pi/3)).*(I_R*0.5).*sawtooth(2.*pi.*50.*t1.*m_f+2*pi/3,1/2);

ic_rip_up=Ip.*sin(2*pi*50*t1+2*pi/3)+(1/2+m/2.*sin(2*50*pi.*t1+2*pi/3)).*(1/2-...
    m/2.*sin(2.*50.*pi.*t1+2*pi/3)).*(I_R*0.5);

ic_rip_down=Ip.*sin(2*pi*50*t1+2*pi/3)-(1/2+m/2.*sin(2*50*pi.*t1+2*pi/3)).*(1/2-...
    m/2.*sin(2.*50.*pi.*t1+2*pi/3)).*(I_R*0.5);

%% Plotting

figure;
hold on
plot(t1,ia,t1,ib,t1,ic)

plot(t1,ia_rip_real,t1,ib_rip_real,t1,ic_rip_real)

plot(t1,ia_rip_real,t1,ia_rip_up,t1,ia_rip_down)

plot(t1,ib_rip_real,t1,ib_rip_up,t1,ib_rip_down)

plot(t1,ic_rip_real,t1,ic_rip_up,t1,ic_rip_down)
grid on

figure;
hold on
plot(t1,ia)

plot(t1,ia_rip_real)

plot(t1,ia_rip_real,t1,ia_rip_up,t1,ia_rip_down)


grid on

figure;
hold on
plot(t1,(1/2+m/2.*sin(2*50*pi.*t1)).*(1/2-...
    m/2.*sin(2.*50.*pi.*t1)).*(I_R*0.5).*sawtooth(2.*pi.*50.*t1.*m_f,1/2));
