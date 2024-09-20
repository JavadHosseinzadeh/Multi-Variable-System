function [Alfa,idz,odz,iodz,Poles_G]=Decoupling_Zeros(P,Q,R) 
  
%%                                                        Finding Poles of System: |P(s)|=0 
s=sym('s');                                                                % Symbolic parameter 
Alfa=roots(sym2poly(simplify(det(P))));          
 
%%                                                  Finding Input Decoupling Zeros of System: idz 
S_idz=[P Q]; 
[SMMilan_idz,P_idz,Z_idz]=SmithMacMilan(S_idz); 
idz=roots(sym2poly(Z_idz)); 
 
%%                                                 Finding Output Decoupling Zeros of System: odz 
S_odz=[P ; R]; 
[SMMilan_odz,P_odz,Z_odz]=SmithMacMilan(S_odz); 
odz=roots(sym2poly(Z_odz)); 
 
%%                                                 Finding Poles of G(s) using Smith-MacMilan Form 
G=simplify(R*inv(P)*Q); 
[SMMilan_G,P_G,Z_G]=SmithMacMilan(G); 
sym2poly(P_G); Poles_G=roots(sym2poly(P_G)); 
 
%%                                               Finding Input-Output Decoupling Zeros of System: iodz 
iodz=Poles_G; 
L_iodz=length(iodz); L_idz=length(idz); L_odz=length(odz); 
iodz(L_iodz+1:L_iodz+L_idz,1)=idz; 
L_iodz=length(iodz); iodz(L_iodz+1:L_iodz+L_odz,1)=odz; 
iodz=sort(iodz); Alfa=sort(Alfa); 
iodz2=iodz; Alfa2=Alfa; Alfa_backup=Alfa; 
clear iodz; clear Alfa 
 if length(iodz2)~=length(Alfa2) 
    while length(Alfa2)>0 
        j=0; 
        while j<length(Alfa2) 
            j=j+1; 
            if iodz2(j,1)==Alfa2(1,1) 
                clear iodz; clear Alfa 
                iodz(1:j-1,1)=iodz2(1:j-1); iodz(j:length(iodz2)-1,1)=iodz2(j+1:length(iodz2)); 
                Alfa(1:j-1,1)=Alfa2(1:j-1); Alfa(j:length(Alfa2)-1,1)=Alfa2(j+1:length(Alfa2));               
                clear iodz2; clear ALfa2 
                iodz2=iodz; Alfa2=Alfa; j=0; 
            end 
        end 
    end 
else 
    iodz=[]; 
end 
 Alfa=Alfa_backup; 
end 
%%                                                               The End of Program. 
