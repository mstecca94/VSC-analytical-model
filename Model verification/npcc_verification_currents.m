function [Irms,Iavg,Irms_ripple,I_rms_diff]=npcc_verification_currents(m,I__ac,ddeg,I__rp)

d=ddeg/180.*pi;

%% RMS currents

irms_t23=sqrt(3).*sqrt(-2.*m.*cos(d).^2+4.*m.*cos(d)+3.*pi-2.*m).*abs(I__ac)/(6.*sqrt(pi));

irms_t14=sqrt(2).*sqrt(3).*(cos(d)+1).*sqrt(m).*abs(I__ac)/(6.*sqrt(pi));

irms_d1234=-sqrt(2).*sqrt(3).*(cos(d)-1).*sqrt(m).*abs(I__ac)/(6.*sqrt(pi));

irms_d56=sqrt(3).*sqrt(-(2.*(m.*cos(2.*d)+3.*m-3.*pi.*(1/2))).*I__ac^2)/(6.*sqrt(pi));

%% AVG currents

iavg_t23=-I__ac.*(-m.*cos(abs(d)).*abs(d)+m.*sin(abs(d))-4)/(4.*pi);

iavg_t14=-I__ac.*((-pi+abs(d)).*cos(abs(d))-sin(abs(d))).*m/(4.*pi);

iavg_d1234=m.*I__ac.*(-cos(abs(d)).*abs(d)+sin(abs(d)))/(4.*pi);

iavg_d56=I__ac.*(m.*(abs(d)-(1/2).*pi).*cos(abs(d))-m.*sin(abs(d))+2)/(2.*pi);

%% All currents

Irms=[irms_t14 irms_t23 irms_d1234 irms_d56];

Iavg=[iavg_t14 iavg_t23 iavg_d1234 iavg_d56];


%% with ripple and m / delta / Iac 

irms_t14_ripple2=sqrt((4*I__rp^2*cos(abs(d))^5*m^4*(1/15)-(1/9)*(8*(m^2-3*m*sin(abs(d))*(1/4)+...
    1/2))*I__rp^2*m^2*cos(abs(d))^3+64*I__ac^2*cos(abs(d))^2*(1/3)+(4*I__rp^2*m^4*(1/3)-...
    5*m^3*I__rp^2*sin(abs(d))*(1/3)+4*I__rp^2*m^2*(1/3)+128*I__ac^2*(1/3))*cos(abs(d))+...
    32*I__rp^2*m^4*(1/45)+I__rp^2*(-pi+abs(d))*m^3+8*I__rp^2*m^2*(1/9)+64*I__ac^2*(1/3))*m/(128*pi));

irms_t23_ripple2=sqrt((48*I__rp^2*cos(abs(d))^5*m^5-160*I__rp^2*(m-3*sin(abs(d))*(1/16))...
    *m^4*cos(abs(d))^3-1920*cos(abs(d))^2*m+((-75*I__rp^2*m^4-60*I__rp^2*m^2)*sin(abs(d))+...
    240*I__rp^2*m^5+3840*m)*cos(abs(d))+45*I__rp^2*(abs(d)-2*pi)*m^4+160*I__rp^2*m^3+...
    60*I__rp^2*m^2*abs(d)-1920*m+2880*pi)/(11520*pi));

irms_d1234_ripple2=sqrt(abs((-4*I__rp^2*cos(abs(d))^5*m^4*(1/15)+8*m^2*(m^2+3*m*sin(abs(d))*(1/4)+...
    1/2)*I__rp^2*cos(abs(d))^3*(1/9)+64*I__ac^2*cos(abs(d))^2*(1/3)+(-4*I__rp^2*m^4*(1/3)-...
    5*m^3*I__rp^2*sin(abs(d))*(1/3)-4*I__rp^2*m^2*(1/3)-128*I__ac^2*(1/3))*cos(abs(d))+...
    I__rp^2*abs(d)*m^3+32*I__rp^2*m^4*(1/45)+8*I__rp^2*m^2*(1/9)+64*I__ac^2*(1/3))*m^2/(128*pi)));

irms_d56_ripple2=sqrt((-120*I__rp^2*m^3*(m*sin(abs(d))-4/3)*cos(abs(d))^3-3840*I__ac^2*cos(abs(d))^2*m+...
    300*I__rp^2*m^3*(m*sin(abs(d))-8/5)*cos(abs(d))-128*I__rp^2*m^5-180*I__rp^2*(abs(d)-...
    3*pi*(1/4))*m^4-160*I__rp^2*m^3+60*I__rp^2*m^2*pi-3840*I__ac^2*m+2880*pi*I__ac^2)/(11520*pi));

Irms_ripple=[irms_t14_ripple2 irms_t23_ripple2 irms_d1234_ripple2 irms_d56_ripple2];

I_rms_diff=[ (irms_t14_ripple2-irms_t14)/irms_t14*100 (irms_t23_ripple2-irms_t23)/irms_t23*100 ...
(irms_d1234_ripple2-irms_d1234)/irms_d1234*100 (irms_d56_ripple2-irms_d56)/irms_d56*100];

end