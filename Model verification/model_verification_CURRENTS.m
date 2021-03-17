ccc
%% Design parameters
figindx=1;
I__ac=1;
m=1;

%% Sweeping delta 

delta_1=-0;
delta_2=180;

i=0;

for d=delta_1:0.1:delta_2
i=i+1;
ddeg=d;
I__rp=1;

[Irms,Iavg,Irms_ripple,I_rms_diff]=npcc_verification_currents(m,I__ac,ddeg,I__rp);

%%% t14 t23 d1234 d56%%%
iavg(i,:)=Iavg;

irms(i,:)=Irms;

irms_ripple(i,:)=Irms_ripple;

i_rms_diff(i,:)=I_rms_diff;

end

figure(figindx)
hold on
subplot(5,1,1)
plot(irms(:,1))
plot(irms_ripple(:,1))

subplot(5,1,2)
hold on
plot(irms(:,2))
plot(irms_ripple(:,2))

subplot(5,1,3)
hold on
plot(irms(:,3))
plot(irms_ripple(:,3))

subplot(5,1,4)
hold on
plot(irms(:,4))
plot(irms_ripple(:,4))

subplot(5,1,5)
hold on
plot(i_rms_diff)

figindx=figindx+1;
%% data from plex

ddd_avg=[0 30 60 90 120 180];

it1rms=[0.460659 0.429799 0.345488 0.230317 0.115146 0];

it2rms=[0.5 0.499046 0.486554 0.443785 0.36143 0.194405];

id2rms=[0 0.0308347 0.115146 0.230317 0.345488 0.460659 ];

id5rms=[0.194405 0.253614 0.342597 0.379341 0.342597 0.194405];

it1avg=[0.25 0.220207 0.152244 0.07957 0.027244  0.0];

it2avg=[0.31831 0.314603 0.291057 0.238729 0.166057 0.0683099];

id2avg=abs([0 -0.00370034 -0.027244 -0.0795716 -0.152244 -0.25]);

id5avg=[0.0683099 0.0943963 0.138813 0.159158 0.138813 0.0683099];

ddd2=round(ddd_avg,1);
drange2=ddd2(1):0.01:ddd2(end);
iavg=[  transpose(delta_1:0.1:delta_2)  iavg];
irms=[  transpose(delta_1:0.1:delta_2)  irms];

%% inserting ripple

ddd_avg2=[0 0 0 0];

it1rms2=[0.460659 0.455396 0.455968 0.457144];

it2rms2=[0.5 0.494771 0.496397 0.500366];

id2rms2=[0 0 0.000388 0.000444092];

id5rms2=[0.194405 0.193426 0.196222 0.20336];

it1avg2=[0.25 0.247108 0.2407061 0.246899];

it2avg2=[0.31831 0.315113  0.316728 0.320605];

id2avg2=abs([0  0 0 0]);

id5avg2=[0.068399 0.0680043  0.069667 0.0737061];

ripple_max=[0 0.222785 0.451834 0.73715];

ripple_rms=[0 0.0417707 0.0832839 0.139709];

m_r=[1 0.6677 0.688971 0.732511];

Lg=[0 0.0005 0.00025 0.00015]; 

delt=[0 -4 -4 -4 ];
    
it1rms2=(it1rms2-it1rms2(1))/it1rms2(1)*100;
it2rms2=(it2rms2-it2rms2(1))/it2rms2(1)*100;
id2rms2=(id2rms2-id2rms2(1))/id2rms2(1)*100;
id5rms2=(id5rms2-id5rms2(1))/id5rms2(1)*100;

it1avg2=(it1avg2-it1avg2(1))/it1avg2(1)*100;
it2avg2=(it2avg2-it2avg2(1))/it2avg2(1)*100;
id2avg2=(id2avg2-id2avg2(1))/id2avg2(1)*100;
id5avg2=(id5avg2-id5avg2(1))/id5avg2(1)*100;
%% Plotting current with ripple



ddeg=0;

i=0;

ripple_range=[0 60 120 200];

for a=1:4
    
i=i+1;

ddeg=delt(i);

m=m_r(i);

I__rp=a;
[Irms,Iavg,Irms_ripple,Irms_ripple2]=npcc_verification_currents(m,I__ac,ddeg,I__rp);

%%% t14 t23 d1234 d56%%%
iavg_r(i,:)=round(Iavg,3);
irms_r(i,:)=round(Irms,3);

irms_t14_ripple1(i,1)=round(Irms_ripple2(1),3);
irms_t23_ripple1(i,1)=round(Irms_ripple2(2),3);
irms_d1234_ripple1(i,1)=round(Irms_ripple2(3),3);
irms_d56_ripple1(i,1)=round(Irms_ripple2(4),3);

irms_t14_ripple1_var(i,1)=round((irms_t14_ripple1(i,1)-irms_r(i,1))/irms_r(i,1)*100,3);
irms_t23_ripple1_var(i,1)=round((irms_t23_ripple1(i,1)-irms_r(i,2))/irms_r(i,2)*100,3);
irms_d1234_ripple1_var(i,1)=round((irms_d1234_ripple1(i,1)-irms_r(i,3))/irms_r(i,3)*100,3);
irms_d56_ripple1_var(i,1)=round((irms_d56_ripple1(i,1)-irms_r(i,4))/irms_r(i,4)*100,3);


end

%%
figure(figindx);
hold on
title('I_{rms}')
plot(ripple_range,irms_t14_ripple1_var)
plot(ripple_range,irms_t23_ripple1_var)
plot(ripple_range,irms_d1234_ripple1_var)
plot(ripple_range,irms_d56_ripple1_var)
legend('T1r','T2r','D2r','D5r')
figindx=figindx+1;



figure(figindx);
hold on
title('I_{rms} - Switches')
plot(ripple_max,it1rms2,'--s')
plot(ripple_max,it2rms2,'--s')
plot(ripple_max,id2rms2,'--s')
plot(ripple_max,id5rms2,'--s')
plot(ripple_max,irms_t14_ripple1_var)
plot(ripple_max,irms_t23_ripple1_var)
plot(ripple_max,irms_d1234_ripple1_var)
plot(ripple_max,irms_d56_ripple1_var)
legend('T1-s','T2-s','D2-s','D5-s','T1-c','T2-c','D2-c','D5-c')
figindx=figindx+1;


figure(figindx);
hold on
title('I_{avg} - Switches')
plot(ripple_max,it1avg2,'--s')
plot(ripple_max,it2avg2,'--s')
plot(ripple_max,id2avg2,'--s')
plot(ripple_max,id5avg2,'--s')
legend('T1','T2','D2','D5')
figindx=figindx+1;

%% Plotting current without ripple
figure(figindx);
hold on
title('I_{rms}')
plot(ddd_avg,it1rms,'--s')
plot(ddd_avg,it2rms,'--s')
plot(irms(:,1),irms(:,2:3))
plot(ddd_avg,id2rms,'--s')
plot(ddd_avg,id5rms,'--s')
plot(irms(:,1),irms(:,4:5))
legend('T1-s','T2-s','T1-c','T2-c','D2-s','D5-s','D2-c','D5-c')

figindx=figindx+1;

figure(figindx);
hold on
title('I_{avg}')
plot(ddd_avg,it1avg,'--s')
plot(ddd_avg,it2avg,'--s')
plot(iavg(:,1),iavg(:,2:3))
plot(ddd_avg,id2avg,'--s')
plot(ddd_avg,id5avg,'--s')
plot(iavg(:,1),iavg(:,4:5))
legend('T1-s','T2-s','T1-c','T2-c','D2-s','D5-s','D2-c','D5-c')

figindx=figindx+1;

%% plecs data 2

% ddd_avg=[-0.456647 26.1995 44.7711 56.1694 75.5061 85.8894];
% 
% ilinerms=[14.1427 15.7614 19.9214 25.4052 29.158 33.3048];
% 
% it1rms=[7.52217 8.0589 9.52925 11.9646 3.84894 3.98608];
% 
% it2rms=[9.98631 11.0562 13.7054 16.9 20.4076 23.2598];
% 
% id2rms=[0.277797 1.37662 3.30546 6.13824 2.86717 3.61181];
% 
% id5rms=[6.56836 7.56928 9.85046 11.9357 20.0413 22.9157];
% 
% it1avg=[3.35502 3.56059 4.11922 5.18805 0.505214 0.49381];
% 
% it2avg=[6.34398 6.86275 8.8893 9.60045 12.7906 14.5489];
% 
% id2avg=abs([-0.0204574 -0.226659 -0.78852 -1.85699 -0.318512 -0.422192]);
% 
% id5avg=[2.98896 3.30215 4.06972 4.4124 12.2854 14.055];
% 
% ripplemax=[0.4429242 0.451433 0.480604 0.492748 1.75007 1.88518];
% 
% ripplerms=[0.140694 0.132701 0.132926 0.152765 0.495842 0.537357];
% 
% for i=1:1:6
%     
%     it1rms(i)=it1rms(i)/ilinerms(i)/sqrt(2);
%     it2rms(i)=it2rms(i)/ilinerms(i)/sqrt(2);
%     id2rms(i)=id2rms(i)/ilinerms(i)/sqrt(2);
%     id5rms(i)=id5rms(i)/ilinerms(i)/sqrt(2);
%     
%     it1avg(i)=it1avg(i)/ilinerms(i)/sqrt(2);
%     it2avg(i)=it2avg(i)/ilinerms(i)/sqrt(2);
%     id2avg(i)=id2avg(i)/ilinerms(i)/sqrt(2);
%     id5avg(i)=id5avg(i)/ilinerms(i)/sqrt(2);
% 
% end
% 
% %%
% ddd2=round(ddd_avg,1);
% drange2=ddd2(1):0.01:ddd2(end);
% iavg=[  transpose(delta_1:0.1:delta_2)  iavg];
% irms=[  transpose(delta_1:0.1:delta_2)  irms];
% 
% %%
% figure(figindx);
% hold on
% title('I_{rms} - Switches')
% plot(ddd_avg,it1rms,'--s')
% plot(ddd_avg,it2rms,'--s')
% plot(irms(:,1),irms(:,2:3))
% legend('T1-s','T2-s','T1-c','T2-c')
% 
% figindx=figindx+1;
% 
% figure(figindx);
% hold on
% title('I_{rms} - Diodes')
% plot(ddd_avg,id2rms,'--s')
% plot(ddd_avg,id5rms,'--s')
% plot(irms(:,1),irms(:,4:5))
% legend('D2-s','D5-s','D2-c','D5-c')
% 
% figindx=figindx+1;
% 
% figure(figindx);
% hold on
% title('I_{avg} - Switches')
% plot(ddd_avg,it1avg,'--s')
% plot(ddd_avg,it2avg,'--s')
% plot(iavg(:,1),iavg(:,2:3))
% legend('T1-s','T2-s','T1-c','T2-c')
% 
% figindx=figindx+1;
% 
% figure(figindx);
% hold on
% title('I_{avg}- Diodes')
% plot(ddd_avg,id2avg,'--s')
% plot(ddd_avg,id5avg,'--s')
% plot(iavg(:,1),iavg(:,4:5))
% legend('D2-s','D5-s','D2-c','D5-c')
% 
% figindx=figindx+1;
