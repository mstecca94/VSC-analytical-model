ccc
% d=0;%-pi/3;


for d=0:pi/12:pi
i_a=@(x) sin(x)*(1/2+1/2*sin(x+d));
i_b=@(x) sin(x-2*pi/3)*(1/2+1/2*sin(x-2*pi/3+d));
i_c=@(x) sin(x+2*pi/3)*(1/2+1/2*sin(x+2*pi/3+d));



i_dc=@(x) i_a(x) + i_b(x) + i_c(x);
figure(1)
hold on
fplot(i_dc,[-pi 2*pi])
grid on

figure(2)
hold on
grid on
fplot(i_a,[-pi 2*pi])

figure(3)
hold on
grid on
fplot(i_b,[-pi 2*pi])
figure(4)
hold on
grid on
fplot(i_c,[-pi 2*pi])


end

% figure;
% hold on
% grid on
% fplot(i_a,[-pi 2*pi])
% 
% figure;
% hold on
% fplot(i_b,[-pi 2*pi])
% 
% grid on
% 
% figure;
% hold on
% fplot(i_c,[-pi 2*pi])
% grid on
% 
% figure;
% hold on
% 
% fplot(i_a,[-pi 2*pi])
% fplot(i_b,[-pi 2*pi])
% fplot(i_c,[-pi 2*pi])
% grid on

