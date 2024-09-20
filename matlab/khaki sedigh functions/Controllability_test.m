function [Mu,Controllability_Index,Partial_Controllability_Matix]=Controllability_test(A,B) 
  
%%                                        Cheking Controllabillity of System using Controllability Matrix 
Fi_c=ctrb(A,B);                                                            % Controllability matrix 
if rank(Fi_c)==size(A,1) 
    disp('The System is Controllable!'); flag=1; 
else 
    Rank_Deffiency_C=size(A,1)-rank(Fi_c); 
    disp(['The System is not Controllable! ,  Rank defficiency = ',num2str(Rank_Deffiency_C)]); 
    flag=0; 
end 
 
%%                                Controllability Indices, Controllability Index and Partial Controllability Matrix 
if flag==1 
    [n,m]=size(B); count=0; 
    Mu(1:m)=zeros(m,1);                                                    % Mu:Controllability Indices 
    for i=1:n 
        for j=1:m 
            count=count+1; 
            if count==1 
                F=Fi_c(:,1); Mu(j)= Mu(j)+1; 
            else 
                if rank([F Fi_c(:,count)]) > rank(F) 
                    F=[F Fi_c(:,count)]; Mu(j)= Mu(j)+1; 
                end 
            end 
        end 
    end 
    Controllability_Index=max(Mu);                                         % Controllability Index 
    Partial_Controllability_Matix=B;                                       % Partial Controllability Matix 
    for i=1:Controllability_Index-1 
        Partial_Controllability_Matix=[Partial_Controllability_Matix (A^i)*B]; 
    end 
else                                                                       % The case of uncontrollability 
   Mu=[]; Controllability_Index=[]; Partial_Controllability_Matix=[];  
end 
end 
%%                                                               The End of Program. 
