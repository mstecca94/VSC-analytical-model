function [ Rs ]= heat_sink_design ( Ptot_comp , Semicond )

%% Load data

[ Semicond ] = semiconductor_data ( Semicond  ) ;

%% here I can add different switches and switch to compare them

Tjmax = 150;
Tamb = 40 ;


topology = Semicond.topology ;
%% switch topology

switch topology
    case '2LC'   
        % Find lenght of heat sink 
        P_Tot =  sum(Ptot_comp)*3 ;
        P_T =  Ptot_comp(1) ;
        P_D =  Ptot_comp(2) ;
        RjcT = Semicond.Param.RjcT ;
        RjcD = Semicond.Param.RjcD ;
        switch Semicond.Type 
            case 'FF300R12KE4'
                RcsT = Semicond.Param.RcsT ;
                RcsD = Semicond.Param.RcsD ;                
                RsT = ( Tjmax - Tamb - P_T*(RjcT + RcsT) ) / P_Tot ;
                RsD = ( Tjmax - Tamb - P_D*(RjcD + RcsD) ) / P_Tot ;
            case 'cas300m12bm2' 
                P_Module =  sum(Ptot_comp) ;
                Rcs = Semicond.Param.Rcs ;                
                RsT = ( Tjmax - Tamb - P_T*RjcT - P_Module*Rcs ) / P_Tot ;
                RsD = ( Tjmax - Tamb - P_D*RjcD - P_Module*Rcs ) / P_Tot ;
        end

        Rs = min (RsT,RsD) ;

    case '3L_NPC'
        RjcT1 = Semicond.Param.RjcT1 ;
        RjcT2 = Semicond.Param.RjcT2 ;
        RjcD1 = Semicond.Param.RjcD1 ;
        RjcD2 = Semicond.Param.RjcD2 ;
        RjcD5 = Semicond.Param.RjcD5 ;
        RcsT1 = Semicond.Param.RcsT1 ;
        RcsT2 = Semicond.Param.RcsT2 ;
        RcsD1 = Semicond.Param.RcsD1 ;
        RcsD2 = Semicond.Param.RcsD2 ;
        RcsD5 = Semicond.Param.RcsD5 ;
        % Find lenght of heat sink 
        P_Tot =  sum(Ptot_comp)*3 ;
        P_T1 =  Ptot_comp(1) ;
        P_T2 =  Ptot_comp(3) ;
        P_D1 =  Ptot_comp(2) ;
        P_D2 =  Ptot_comp(4) ;
        P_D5 =  Ptot_comp(9) ;
        

        RsT1 = ( Tjmax - Tamb - P_T1*(RjcT1 + RcsT1) ) / P_Tot ;
        RsT2 = ( Tjmax - Tamb - P_T2*(RjcT2 + RcsT2) ) / P_Tot ;
        RsD1 = ( Tjmax - Tamb - P_D1*(RjcD1 + RcsD1) ) / P_Tot ;
        RsD2 = ( Tjmax - Tamb - P_D2*(RjcD2 + RcsD2) ) / P_Tot ;
        RsD5 = ( Tjmax - Tamb - P_D5*(RjcD5 + RcsD5) ) / P_Tot ;
        Rs = min ( [ RsT1 RsT2 RsD1 RsD2 RsD5 ]) ;


    case '3L_TTYPE'
        RjcT1 = Semicond.Param.RjcT1 ;
        RjcT2 = Semicond.Param.RjcT2 ;
        RjcD1 = Semicond.Param.RjcD1 ;
        RjcD2 = Semicond.Param.RjcD2 ;
        RcsT1 = Semicond.Param.RcsT1 ;
        RcsT2 = Semicond.Param.RcsT2 ;
        RcsD1 = Semicond.Param.RcsD1 ;
        RcsD2 = Semicond.Param.RcsD2 ;
        % Find lenght of heat sink 
        P_Tot =  sum(Ptot_comp)*3 ;
        P_T1 =  Ptot_comp(1) ;
        P_T2 =  Ptot_comp(3) ;
        P_D1 =  Ptot_comp(2) ;
        P_D2 =  Ptot_comp(4) ;
        
        RsT1 = ( Tjmax - Tamb - P_T1*(RjcT1 + RcsT1) ) / P_Tot ;
        RsT2 = ( Tjmax - Tamb - P_T2*(RjcT2 + RcsT2) ) / P_Tot ;
        RsD1 = ( Tjmax - Tamb - P_D1*(RjcD1 + RcsD1) ) / P_Tot ;
        RsD2 = ( Tjmax - Tamb - P_D2*(RjcD2 + RcsD2) ) / P_Tot ;
        Rs = min ([ RsT1 RsT2 RsD1 RsD2]) ;

end

end