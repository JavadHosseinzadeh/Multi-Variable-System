function Reduced_Sys=Order_Reduction(A,B,C,D,k,option) 
warning off 
Method_of_Order_Reduction=option; 
Sys_Canonical=canon(ss(A,B,C,D)); 
Ac=Sys_Canonical.A;Bc=Sys_Canonical.B;Cc=Sys_Canonical.C;Dc=Sys_Canonical.D; 
system=pck(Ac,Bc,Cc,Dc); 
sys1 = strans(system); 
switch Method_of_Order_Reduction 
    case 'Truncation' 
        Sys_Trun = strunc(sys1,k);                                         % Truncation 
        [At,Bt,Ct,Dt] = unpck(Sys_Trun); 
        Reduced_Sys = ss(At,Bt,Ct,Dt); 
    case 'Residualization' 
        Sys_Res = sresid(sys1,k);                                          % Residualization 
        [Ar,Br,Cr,Dr] = unpck(Sys_Res); 
        Reduced_Sys = ss(Ar,Br,Cr,Dr); 
    case 'Balanced Truncation' 
        sys2 = sysbal(system);                                             % balanced model 
        btru = strunc(sys2,k);                                             % balanced Truncation  
        [Abt,Bbt,Cbt,Dbt] = unpck(btru); 
        Reduced_Sys = ss(Abt,Bbt,Cbt,Dbt); 
    case 'Balanced Residualization' 
        sys2 = sysbal(system);                                             % balanced model 
        bres= sresid(sys2,k);                                              % balanced Residualization  
        [Abr,Bbr,Cbr,Dbr] = unpck(bres); 
        Reduced_Sys = ss(Abr,Bbr,Cbr,Dbr); 
end 
end 
%%                                                               The End of Program. 
% Author: Peyman Bagheri 
% email: peyman.bk@gmail.com 
% Created: September 2011;  
% Last revision: 09-November-2011; 
