function [L1lc,Clc,irppactlc,lc_feas]=filter_lc(L1_min,wres,Cmax,Lmin,Ltot_max)

global vg fs Ig

i2=0;

for t=L1_min:0.0000001:Ltot_max
    i2=i2+1;

    L1lc(i2)=t;

    Clc(i2)=1/(wres^2*L1lc(i2));

    irppactlc(i2)=vg/(2*sqrt(6)*fs*L1lc(i2))/Ig;
    
    if  Clc(i2)<=Cmax 
        
        lc_feas(i2)=1;
        
        continue
   
    end

    if Clc(i2)>=Cmax 
        
%         disp('Maximum Capacitance Reached')

        Clc(i2)=Cmax*0.9;
        
        L1lc(i2)=1/(wres^2*Clc(i2));

        irppactlc(i2)=vg/(2*sqrt(6)*fs*L1lc(i2))/Ig;
        
        lc_feas(i2)=0;
        
    end
    
    
end

end