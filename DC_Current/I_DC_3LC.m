
ccc

m_f=99;

t=0:0.0000001:.02;


xa_3l_up = 0.5+0.5*sawtooth(2*pi*50*t*m_f,1/2);
xa_3l_down = -0.5+0.5*sawtooth(2*pi*50*t*m_f,1/2);


% figure;
% plot(t,xa_3l_up,t,xa_3l_down)
% grid on
%%
va_ref=sin(2*pi*50*t);
vb_ref=sin(2*pi*50*t+pi*2/3);
vc_ref=sin(2*pi*50*t-pi*2/3);

I_p=1;

% figure;
% plot(t,xa_3l_up,t,xa_3l_down,t,va_ref,t,vb_ref,t,vc_ref)
% grid on

d=0;
tt=0;


for d=-180:1:180
d=d/180*pi;
tt=tt+1;

i_a=I_p*sin(2*pi*50*t+d);
i_b=I_p*sin(2*pi*50*t+pi*2/3+d);
i_c=I_p*sin(2*pi*50*t-pi*2/3+d);
i=1;
for i=1:1:size(va_ref,2)
    
    if va_ref(i)>0
        if va_ref(i)>=xa_3l_up(i)
            t1_a(i)=1;
            t4_a(i)=0;
        else
            t1_a(i)=0;
            t4_a(i)=0;
        end
    else
        if va_ref(i)>=xa_3l_down(i)
            t1_a(i)=0;
            t4_a(i)=0;
        else
            t1_a(i)=0;
            t4_a(i)=1;
        end
    end
    
    if vb_ref(i)>0

        if vb_ref(i)>=xa_3l_up(i)
            t1_b(i)=1;
            t4_b(i)=0;
        else
            t1_b(i)=0;
            t4_b(i)=0;
        end
    else
        if vb_ref(i)>=xa_3l_down(i)
            t1_b(i)=0;
            t4_b(i)=0;
        else
            t1_b(i)=0;
            t4_b(i)=1;
       end
    end
            
    if vc_ref(i)>0

        if vc_ref(i)>=xa_3l_up(i)
            t1_c(i)=1;
%             t2_c(i)=1;
%             t3_c(i)=1;
            t4_c(i)=0;
        else
            t1_c(i)=0;
            t4_c(i)=0;
        end
    
    else
        if vc_ref(i)>=xa_3l_down(i)
            t1_c(i)=0;
            t4_c(i)=0;
        else
            t1_c(i)=0;
            t4_c(i)=1;
        end
    
    end            
            
end

i_leg_a_p=i_a.*t1_a;
i_leg_b_p=i_b.*t1_b;
i_leg_c_p=i_c.*t1_c;

i_leg_a_n=i_a.*t4_a;
i_leg_b_n=i_b.*t4_b;
i_leg_c_n=i_c.*t4_c;

for i=1:1:size(va_ref,2)
    
    if va_ref(i)>0
        i_leg_a_0(i)=i_a(i)-i_leg_a_p(i);
    else
        i_leg_a_0(i)=i_a(i)-i_leg_a_n(i);
    end
    
    if vb_ref(i)>0
        i_leg_b_0(i)=i_b(i)-i_leg_b_p(i);
    else
        i_leg_b_0(i)=i_b(i)-i_leg_b_n(i);
    end
    
    if vc_ref(i)>0
        i_leg_c_0(i)=i_c(i)-i_leg_c_p(i);
    else
        i_leg_c_0(i)=i_c(i)-i_leg_c_n(i);
    end
 

end

i_p=i_leg_a_p+i_leg_b_p+i_leg_c_p;
i_0=i_leg_a_0+i_leg_b_0+i_leg_c_0;
i_n=i_leg_a_n+i_leg_b_n+i_leg_c_n;

i_dc=i_p+i_0/2;

% figure;
% subplot(3,1,1)
% plot(t,t1_a,t,t4_a)
% grid on
% subplot(3,1,2)
% plot(t,t1_b,t,t4_b)
% grid on
% subplot(3,1,3)
% plot(t,t1_c,t,t4_c)
% grid on
% 
% 
% figure;
% sgtitle('Leg A')
% subplot(3,1,1)
% plot(t,i_leg_a_p)
% grid on
% subplot(3,1,2)
% plot(t,i_leg_a_n)
% grid on
% subplot(3,1,3)
% plot(t,i_leg_a_0)
% grid on
% 
% 
% figure;
% sgtitle('Leg B')
% subplot(3,1,1)
% plot(t,i_leg_b_p)
% grid on
% subplot(3,1,2)
% plot(t,i_leg_b_n)
% grid on
% subplot(3,1,3)
% plot(t,i_leg_b_0)
% grid on
% 
% figure;
% sgtitle('Leg C')
% subplot(3,1,1)
% plot(t,i_leg_c_p)
% grid on
% subplot(3,1,2)
% plot(t,i_leg_c_n)
% grid on
% subplot(3,1,3)
% plot(t,i_leg_c_0)
% grid on
% 
% figure;
% sgtitle('I_{p} I_{n} I_{0}')
% subplot(4,1,1)
% plot(t,i_p)
% grid on
% subplot(4,1,2)
% plot(t,i_n)
% grid on
% subplot(4,1,3)
% plot(t,i_0)
% grid on
% subplot(4,1,4)
% plot(t,i_dc)
% grid on


i_dc_avg(tt)=mean(i_dc);
i_dc_ms(tt)=rms(i_dc)^2;
end
figure;
subplot(3,1,1)
plot(t,t1_a,t,t4_a)
grid on
subplot(3,1,2)
plot(t,t1_b,t,t4_b)
grid on
subplot(3,1,3)
plot(t,t1_c,t,t4_c)
grid on


figure;
sgtitle('Leg A')
subplot(4,1,1)
plot(t,va_ref,t,i_a)
grid on
subplot(4,1,2)
plot(t,i_leg_a_p)
grid on
subplot(4,1,3)
plot(t,i_leg_a_n)
grid on
subplot(4,1,4)
plot(t,i_leg_a_0,t,i_a)
grid on


figure;
sgtitle('Leg B')
subplot(4,1,1)
plot(t,vb_ref,t,i_b)
grid on
subplot(4,1,2)
plot(t,i_leg_b_p)
grid on
subplot(4,1,3)
plot(t,i_leg_b_n)
grid on
subplot(4,1,4)
plot(t,i_leg_b_0,t,i_b)
grid on

figure;
sgtitle('Leg C')
subplot(4,1,1)
plot(t,vc_ref,t,i_c)
grid on
subplot(4,1,2)
plot(t,i_leg_c_p)
grid on
subplot(4,1,3)
plot(t,i_leg_c_n)
grid on
subplot(4,1,4)
plot(t,i_leg_c_0,t,i_c)
grid on

figure;
sgtitle('I_{p} I_{n} I_{0}')
subplot(4,1,1)
plot(t,i_p)
grid on
subplot(4,1,2)
plot(t,i_n)
grid on
subplot(4,1,3)
plot(t,i_0)
grid on
subplot(4,1,4)
plot(t,i_dc)
%% functions of maple

I__ac=1;
m=1;

delta=-pi:1/180*pi:pi;

i_avg_f=3*I__ac*m*cos(delta)*(1/4);


i_ms_tot=I__ac^2.*(8*m*sqrt(3).*cos(d)^2+6*sqrt(3).*cos(d)^2-18*m.*cos(d).^2+...
    2.*m.*sqrt(3)+2*pi-3*sqrt(3))/(4*pi);



%%

grid on
figure;
subplot(2,1,1)
plot(delta,i_dc_avg,delta,i_avg_f)
subplot(2,1,2)
plot(delta,i_dc_ms,delta,i_ms_tot)



%%
m=1;
delta=-pi:1/180*pi:pi;
i_hf_ms=I__ac^2.*(8*m*sqrt(3).*cos(d)^2+6*sqrt(3).*cos(d)^2-18*m.*cos(d).^2+...
    2.*m.*sqrt(3)+2*pi-3*sqrt(3))/(4*pi)-(3*I__ac*m.*cos(delta)*(1/4)).^2;
figure;
plot(delta,i_hf_ms)




%%
[m,delta]=meshgrid(0.05:0.05:0.95,-pi:1/180*pi:pi);

i_hf_rms=sqrt(I__ac^2.*(8.*m.*sqrt(3).*cos(delta).^2+6*sqrt(3).*cos(delta).^2-18*m.*cos(delta).^2+...
    2.*m.*sqrt(3)+2*pi-3*sqrt(3))/(4*pi)-(3*I__ac*m.*cos(delta).*(1/4)).^2);

i_hf_rms2=sqrt(((32*m+24).*cos(delta).^2+8.*m-12)*sqrt(3)+(-9*pi.*m.^2-72.*m).*cos(delta).^2+...
    8*pi)*abs(I__ac)/(4*sqrt(pi));

figure;
surf(m,delta,i_hf_rms2,'FaceAlpha',0.5,'EdgeColor','none');
% colormap(1-gray)
xlabel('m')
ylabel('\delta')
zlabel('Ripple')

figure;
surf(m,delta,i_hf_rms2,'FaceAlpha',0.5,'EdgeColor','none');
% colormap(1-gray)
xlabel('m')
ylabel('\delta')
zlabel('Ripple')