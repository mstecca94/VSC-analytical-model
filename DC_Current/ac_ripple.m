ccc
syms t
m_f=99;
fs=50*m_f;
t1=0:0.0000001:.02;

Ip=200;
I_R=200*0.05*2; % max +/- 5% ripple of the full power RMS current

L=0.4*1e-3;

vp=400/sqrt(3)*sqrt(2);

va_ref=vp*sin(2*pi*50*t);

vdc=800;
m=2*vp/vdc;

i_ripple_2l=(1/2+m/2*sin(2*50*pi*t))*(1-m*sin(2*50*pi*t))*I_R;
i_ripple_2l_b=(1/2+m/2*sin(2*50*pi*t-2*pi/3))*(1-m*sin(2*50*pi*t-2*pi/3))*I_R;
i_ripple_2l_c=(1/2+m/2*sin(2*50*pi*t+2*pi/3))*(1-m*sin(2*50*pi*t+2*pi/3))*I_R;
i_ripple_3l=(m*sin(2*50*pi*t))*(1-m*sin(2*50*pi*t))*I_R;

i2=vdc/2-vp*sin(2*pi*50*t);
ia=Ip*sin(2*pi*50*t);
ib=Ip*sin(2*pi*50*t-2*pi/3);
ic=Ip*sin(2*pi*50*t+2*pi/3);

ia_rip_up=ia+0.5*i_ripple_2l;
ia_rip_down=ia-0.5*i_ripple_2l;
ia_rip_real=Ip.*sin(2*pi*50*t1)+(1/2+m/2.*sin(2*50*pi.*t1)).*(1-...
    m.*sin(2.*50.*pi.*t1)).*(I_R*0.5).*sawtooth(2.*pi.*50.*t1.*m_f,1/2);


ib_rip_up=ib+0.5*i_ripple_2l_b;
ib_rip_down=ib-0.5*i_ripple_2l_b;

ic_rip_up=ic+0.5*i_ripple_2l_c;
ic_rip_down=ic-0.5*i_ripple_2l_c;

i_ab_up=ia_rip_up-ib_rip_up;
i_ab_down=ia_rip_down-ib_rip_down;

I_p=1;

figure;
% plot(t,xa_2l,t,va_ref)
fplot(va_ref,[0 0.02])
grid on

%%
figure;
subplot(3,1,1)
hold on
fplot(i2,[0 0.02])
% fplot(ia,[0 0.02])
subplot(3,1,2)
hold on l
fplot(0.5*i_ripple_2l,[0 0.02])
fplot(-0.5*i_ripple_2l,[0 0.02])
% fplot(i_ripple_3l,[0 0.02])
hold off
subplot(3,1,3)
hold on
fplot(ia_rip_up,[0 0.02])
fplot(ia_rip_down,[0 0.02])
fplot(ia,[0 0.02])
plot(t1,ia_rip_real)
% fplot(i_ab_down,[0 0.02])
% fplot(i_ab_up,[0 0.02])
grid on
hold off
% figure;
% subplot(3,1,1)
% plot(t,ia)
% subplot(3,1,2)
% plot(t,i_ripple_up,t,i_ripple_down)
% subplot(3,1,3)
% plot(t,ia_rip_up,t,ia_rip_down)
% grid on

