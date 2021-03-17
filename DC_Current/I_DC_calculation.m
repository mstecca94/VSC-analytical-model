
ccc

m_f=39;

t=0:0.0000001:.02;



xa = sawtooth(2*pi*50*t*m_f+pi/2,1/2);
xb = sawtooth(2*pi*50*t*m_f+pi*2/3+pi/2,1/2);
xc = sawtooth(2*pi*50*t*m_f-pi*2/3+pi/2,1/2);

% figure;
% plot(t,xa,t,xb,t,xc)
% grid on

xa_3l_up = 0.5+0.5*sawtooth(2*pi*50*t*m_f,1/2);
xa_3l_down = -0.5+0.5*sawtooth(2*pi*50*t*m_f,1/2);


% figure;
% plot(t,xa_3l_up,t,xa_3l_down)
% grid on
%%
ma=sin(2*pi*50*t);
mb=sin(2*pi*50*t+pi*2/3);
mc=sin(2*pi*50*t-pi*2/3);


figure;
plot(t,xa_3l_up,t,xa_3l_down,t,ma,t,mb,t,mc)
grid on

%%
% for i=1:1:size(ma,2)
%     
%     if ma(i)>=xa(i)
%         gate_a(i)=1;
%     else
%         gate_a(i)=0;
%     end
% 
%     if mb(i)>=xb(i)
%         gate_b(i)=1;
%     else
%         gate_b(i)=0;
%     end
% 
%     if mc(i)>=xc(i)
%         gate_c(i)=1;
%     else
%         gate_c(i)=0;
%     end
%     
% end

for i=1:1:size(ma,2)
    
    if ma(i)>0
        if ma(i)>=xa_3l_up(i)
            gate_a2(i)=1;
            gate_a2_0(i)=0;
        else
            gate_a2(i)=0;
            gate_a2_0(i)=1;
        end
    else
        if ma(i)>=xa_3l_down(i)
            gate_a2(i)=1;
            gate_a2_0(i)=0;
        else
            gate_a2(i)=0;
            gate_a2_0(i)=1;
        end
    end
    
    if mb(i)>0

        if mb(i)>=xa_3l_up(i)
            gate_b2(i)=1;
            gate_b2_0(i)=0;
        else
            gate_b2(i)=0;
            gate_b2_0(i)=1;
        end
    else
        if mb(i)>=xa_3l_down(i)
            gate_b2(i)=1;
            gate_b2_0(i)=0;
        else
            gate_b2(i)=0;
            gate_b2_0(i)=1;
       end
    end
            
    if mc(i)>0

        if mc(i)>=xa_3l_up(i)
            gate_c2(i)=1;
            gate_c2_0(i)=0;
        else
            gate_c2(i)=0;
            gate_c2_0(i)=1;
        end
    
    else
        if mc(i)>=xa_3l_down(i)
            gate_c2(i)=1;
            gate_c2_0(i)=0;
        else
            gate_c2(i)=0;
            gate_c2_0(i)=1;
        end
    
    end            
            
end


% figure;
% subplot(2,1,1)
% plot(t,xa,t,ma)
% grid on
% subplot(2,1,2)
% plot(t,gate_a)
% grid on
% 
% figure;
% subplot(2,1,1)
% plot(t,xb,t,mb)
% grid on
% subplot(2,1,2)
% plot(t,gate_b)
% grid on
% 
% figure;
% subplot(2,1,1)
% plot(t,xc,t,mc)
% grid on
% subplot(2,1,2)
% plot(t,gate_c)
% grid on

d=0;
I_p=1;

i_a=I_p*sin(2*pi*50*t+d);
i_b=I_p*sin(2*pi*50*t+pi*2/3+d);
i_c=I_p*sin(2*pi*50*t-pi*2/3+d);


i_leg_a3_0=i_a.*gate_a2_0;
i_leg_b3_0=i_b.*gate_b2_0;
i_leg_c3_0=i_c.*gate_c2_0;

i_leg_a3=i_a.*gate_a2;
i_leg_b3=i_b.*gate_b2;
i_leg_c3=i_c.*gate_c2;

i_DC=i_leg_a3+i_leg_b3+i_leg_c3+i_leg_a3_0/2+i_leg_b3_0/2+i_leg_c3_0/2;

figure;
plot(i_DC)

figure;
subplot(3,1,1)
plot(t,gate_a2,t,gate_a2_0)
grid on
subplot(3,1,2)
plot(t,gate_b2,t,gate_b2_0)
grid on
subplot(3,1,3)
plot(t,gate_c2,t,gate_c2_0)
grid on

figure;
subplot(3,1,1)
plot(t,i_leg_a3)
grid on
subplot(3,1,2)
plot(t,i_leg_b3)
grid on
subplot(3,1,3)
plot(t,i_leg_c3)
grid on


i_dc_avg=mean(i_DC)
i_dc_rms=rms(i_DC)

% I_p=1;
% 
% for d=-pi:pi/12:pi
% i_a=I_p*sin(2*pi*50*t+d);
% i_b=I_p*sin(2*pi*50*t+pi*2/3+d);
% i_c=I_p*sin(2*pi*50*t-pi*2/3+d);
% 
% i_leg_a=i_a.*gate_a;
% i_leg_b=i_b.*gate_b;
% i_leg_c=i_c.*gate_c;
% 
% i_leg_a2=i_a.*(1/2+1/2*sin(2*pi*50*t));
% i_leg_b2=i_b.*(1/2+1/2*sin(2*pi*50*t+pi*2/3));
% i_leg_c2=i_c.*(1/2+1/2*sin(2*pi*50*t-pi*2/3));
% 
% i_leg_a3=i_a.*gate_a2;
% i_leg_b3=i_b.*gate_b2;
% i_leg_c3=i_c.*gate_c2;
% 
% i_DC=i_leg_a+i_leg_b+i_leg_c;
% i_DC2=i_leg_a2+i_leg_b2+i_leg_c2;
% i_DC3=i_leg_a3+i_leg_b3+i_leg_c3;
% 
% figure;
% subplot(3,1,1)
% plot(t,i_leg_a)
% grid on
% subplot(3,1,2)
% plot(t,i_leg_b)
% grid on
% subplot(3,1,3)
% plot(t,i_leg_c)
% grid on
% 
% figure;
% subplot(3,1,1)
% hold on
% plot(t,i_leg_a)
% plot(t,i_leg_a2)
% plot(t,i_leg_a3)
% grid on
% subplot(3,1,2)
% hold on
% plot(t,i_leg_b)
% plot(t,i_leg_b2)
% plot(t,i_leg_b3)
% grid on
% subplot(3,1,3)
% hold on
% plot(t,i_leg_c)
% plot(t,i_leg_c2)
% plot(t,i_leg_c3)
% grid on
% 
% 
% 
% I_DC_avg=mean(i_DC)
% I_DC_avg2=mean(i_DC3)
% 
% mean_a_p=mean(i_leg_a(1,10000:20000));
% mean_a_n=mean(i_leg_a(1,1:10000));
% mean_a=mean(i_leg_a);
% mean_b=mean(i_leg_b);
% mean_c=mean(i_leg_c);
% 
% I_DC_RMS=rms(i_DC)
% 
% I_DC_RMS3=rms(i_DC3)
% 
% figure(100)
% hold on
% scatter(d,I_DC_avg)
% scatter(d,I_DC_avg2)
% % subplot(2,1,1)
% % hold on
% % plot(t,i_DC)
% % plot(t,i_DC2)
% % subplot(2,1,2)
% % plot(t,i_DC3)
% 
% end
% %%
% a=@(d1) 3/4*cos(d1);
% figure(100)
% hold on
% fplot(a,[-pi 2*pi])