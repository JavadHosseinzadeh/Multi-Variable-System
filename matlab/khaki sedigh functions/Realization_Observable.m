function [A,B,C,D]=Realization_Observable(G) 
 
 %%                                                            Creating Symbolic Form of G(s) 
s=sym('s'); 
switch class(G) 
    case 'tf' 
        G_TF=G; [num,den]=tfdata(G_TF);                                    % Transfer Function Form 
        for i=1:size(num,1)                                                % Number of Outputs 
            for j=1:size(num,2)                                            % Number of Inputs 
                num_Sym(i,j)=poly2sym(num{i,j},'s'); 
                den_Sym(i,j)=poly2sym(den{i,j},'s'); 
            end 
        end 
        G_Sym=simplify(num_Sym./den_Sym); 
    case 'sym'  
        G_Sym=simplify(G);                                                 % Symbolic Form 
end 
[num_Sym,den_Sym]=numden(G_Sym); 
 
%%                                                             Creating D Matrix of G(s) 
[l,m]=size(G_Sym); D=zeros(l,m); 
for i=1:l 
    for j=1:m 
        [q,r]=quorem(num_Sym(i,j),den_Sym(i,j)); 
        if q~=0 
            D(i,j)=q; G_Sym(i,j)=simplify(r/den_Sym(i,j)); 
        end 
    end 
end 
 
%%                                                          Creating A, B and C of G(s) 
for j=1:m 
    den_LCM(1,j)=den_Sym(1,j); 
    for i=1:size(den_Sym,1) 
        den_LCM(1,j)=lcm(den_LCM(1,j),den_Sym(i,j)); 
    end 
    Poly_den{j}=sym2poly(den_LCM(1,j));
    n_A{j}=size(Poly_den{j},2); 
    x=zeros(size(den_Sym,1),n_A{j}); 
    for i=1:size(den_Sym,1) 
        num_LCM(i,j)=simplifyFraction((den_LCM(1,j)/den_Sym(i,j))*num_Sym(i,j)); 
        L=size(sym2poly(num_LCM(i,j)),2);
        x(i,n_A{j}-L+1:n_A{j})=sym2poly(num_LCM(i,j)); 
    end 
    Poly_num{j}=x; a=Poly_den{j}; 
    A{j}=[zeros(n_A{j}-2,1) eye(n_A{j}-2) ; -a(end:-1:2)]; b{j}=[zeros(n_A{j}-2,1);1]; 
    c=Poly_num{j}; Cc{j}=c(:,end:-1:2); 
end 
AA=blkdiag(A{1:end}); BB=blkdiag(b{1:end}); CC=[Cc{1:end}];
A=AA';
B=CC';
C=BB';
%%                                                                 Minimality check 
if rank(ctrb(A,B))<size(A,2) 
    disp('The Realization is not Minimal!') 
else 
    disp('The Realization is Minimal!') 
end 
end 
%%                                                               The End of Program. 
