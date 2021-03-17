close all
I__ac=10;
d=-pi:1/180*pi:pi;
dn=-pi:1/180*pi:0;
dp=0:1/180*pi:pi;

figure(1)
xlim([-pi pi])
hold on
for m=0.05:0.05:0.95
    it1n=sqrt(3)*sqrt((8*m+3*sin(dn)).*cos(dn)-3*dn)*abs(I__ac)/(6*sqrt(pi));

    plot(dn,it1n)
    
end
xlim([-pi pi])

figure(2)
xlim([-pi pi])
hold on
for m=0.05:0.05:0.95
    id1n=    sqrt(3)*sqrt((-8*m-3*sin(dn)).*cos(dn)+3*dn+3*pi)*abs(I__ac)/(6*sqrt(pi));

    plot(dn,id1n)
    
end


figure(1)
xlim([-pi pi])
hold on
for m=0.05:0.05:0.95
    it1p=sqrt(3)*sqrt((8*m-3*sin(dp)).*cos(dp)+3*dp)*abs(I__ac)/(6*sqrt(pi));

    plot(dp,it1p)
    
end


figure(2)
xlim([-pi pi])
hold on
for m=0.05:0.05:0.95
    id1p=sqrt(3)*sqrt((-8*m+3*sin(dp)).*cos(dp)-3*dp+3*pi)*abs(I__ac)/(6*sqrt(pi));

    plot(dp,id1p)
    
end

figure(4)
xlim([-pi pi])
hold on
for m=0.05:0.05:0.95
    it1p=sqrt(3)*sqrt((8*m-3*sin(abs(d))).*cos(abs(d))+3*abs(d))*abs(I__ac)/(6*sqrt(pi));

    plot(d,it1p)
    
end


figure(3)
xlim([-pi pi])
hold on
for m=0.05:0.05:0.95
    id1p=sqrt(3)*sqrt((-8*m+3*sin(abs(d))).*cos(abs(d))-3*abs(d)+3*pi)*abs(I__ac)/(6*sqrt(pi));

    plot(d,id1p)
    
end