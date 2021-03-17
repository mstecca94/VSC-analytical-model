function [L1lcl,L2lcl,Clcl,irppactlcl,lcl_feas]=sel_filter_lcl(L1_min,wres,Cmax,Lmin,Ltot_min,Ltot_max)

global vg fs Ig

i=0;

for t=L1_min:0.0000001:Ltot_max
    
    i=i+1;

    L1lcl(i)=t;

    L2lcl(i)=Ltot_min-L1lcl(i);

    Clcl(i)=(L1lcl(i)+L2lcl(i))/(L1lcl(i)*L2lcl(i)*wres^2);

    irppactlcl(i)=vg/(2*sqrt(6)*fs*L1lcl(i))/Ig;
    
    if  Clcl(i)<=Cmax 
        
        lcl_feas(i)=1;
        
        continue
        
    end

    if  Clcl(i)>=Cmax 
        
        i3=1;
        
        Ltot=L2lcl(i)+L1lcl(i);
        
        while Ltot<=Ltot_max 
            
            
            L2lcl(i)=Ltot_min+(Ltot_max-Ltot_min)*i3/100-L1lcl(i);
            
            Clcl(i)=(L1lcl(i)+L2lcl(i))/(L1lcl(i)*L2lcl(i)*wres^2);

            irppactlcl(i)=vg/(2*sqrt(6)*fs*L1lcl(i))/Ig;
            
            Ltot=L2lcl(i)+L1lcl(i);
            i3=i3+1;
            if Clcl(i)<=Cmax 
                       
                disp('Design not possible')
        
                L1lcl(i)=Ltot_max*1.1;

                L2lcl(i)=Ltot_max*1.1;

                Clcl(i)=Cmax*0.9;

                irppactlcl(i)=1;
                
                lcl_feas(i)=0;
        
                continue
                
            end
            
    end
    
%     if L1lcl(i)<=Lmin || L2lcl(i)<=Lmin
%         
%         disp('Minimum Inductance Step Reached')
%         
%         break
%     end

end

end

