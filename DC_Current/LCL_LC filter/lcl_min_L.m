function [L1lcl,L2lcl,Clcl,irppactlcl]=lcl_min_L(L1_min,wres,Ltot_min,Ltot_max,Cmax)

global vg fs Ig Lmin

i=0;

for t=L1_min:0.0000001:Ltot_min
    
    i=i+1;

    L1lcl(i)=t;

    L2lcl(i)=Ltot_min-L1lcl(i);

    Clcl(i)=(L1lcl(i)+L2lcl(i))/(L1lcl(i)*L2lcl(i)*wres^2);

    irppactlcl(i)=vg/(2*sqrt(6)*fs*L1lcl(i))/Ig;
    
    if L1lcl(i) <= Lmin || L2lcl(i) <= Lmin || Clcl(i)>=Cmax
        
        break
        
    end
     
end

end

