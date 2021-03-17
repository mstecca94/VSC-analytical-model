
close all
T = 10*(1/50);

fs = 100000;
t_s = 0:1/fs:T-1/fs;

xa = sawtooth(2*pi*50*t_s+pi/2,1/2);
xb = sawtooth(2*pi*50*t_s+pi*2/3+pi/2,1/2);
xc = sawtooth(2*pi*50*t_s-pi*2/3+pi/2,1/2);

figure;
plot(t_s,xa,t_s,xb,t_s,xc)
grid on

fs = 1000000;
T_f = (1/50);
t_f = 0:1/fs:T_f-1/fs;

ma=sin(2*pi*50*t_f);
mb=sin(2*pi*50*t_f+pi*2/3);
mc=sin(2*pi*50*t_f-pi*2/3);

figure;
plot(t_f,ma,t_f,mb,t_f,mc)
grid on


for i=1:1:size(ma,2)
    
    if ma(i)>=xa(i)
        gate_a(i)=1;
    else
        gate_a(i)=0;
    end

    if mb(i)>=xb(i)
        gate_b(i)=1;
    else
        gate_b(i)=0;
    end

    if mc(i)>=xc(i)
        gate_c(i)=1;
    else
        gate_c(i)=0;
    end
    
end

figure;
subplot(2,1,1)
plot(t_s,xa,t_f,ma)
grid on
subplot(2,1,2)
plot(t_s,gate_a)
grid on

figure;
subplot(2,1,1)
plot(t_s,xb,t_f,mb)
grid on
subplot(2,1,2)
plot(t_s,gate_b)
grid on

figure;
subplot(2,1,1)
plot(t_s,xc,t_f,mc)
grid on
subplot(2,1,2)
plot(t_s,gate_c)
grid on

I_p=1;

i_a=I_p*sin(2*pi*50*t_f);
i_b=I_p*sin(2*pi*50*t_f+pi*2/3);
i_c=I_p*sin(2*pi*50*t_f-pi*2/3);

i_t1_a=i_a.*gate_a;
i_t1_b=i_b.*gate_b;
i_t1_c=i_c.*gate_c;

i_DC=i_t1_a+i_t1_b+i_t1_c;

figure;
subplot(3,1,1)
plot(t_s,i_t1_a)
grid on
subplot(3,1,2)
plot(t_f,i_t1_a)
grid on
subplot(3,1,3)
plot(t_s,i_t1_a)
grid on

figure;
plot(t_s,i_DC)