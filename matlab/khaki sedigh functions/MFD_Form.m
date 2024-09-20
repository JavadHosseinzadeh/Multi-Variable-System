function [G_Sym,DL,NL,Deg_L,Minimality_L,DR,NR,Deg_R,Minimality_R]=MFD_Form(G)
%symbolic form of G(s)
s=sym('s');
switch class(G)
    case 'tf'
        G_TF=G
        [num,den]=tfdata(G_TF);
        for i=1:size(num,1)
            for j=1:size(num,2)
                num_Sym(i,j)=  poly2sym(num(i,j),s);
                den_Sym(i,j)=  poly2sym(den(i,j),s);
            end
        end
        G_Sym=simplify(num_Sym./den_Sym);
    case 'sym'
        G_Sym=simplify(G);
end
G_Sym_backup=G_Sym;

%% Finding Least common multiple in denominator and representing G(s) in the form of: G=(1/d)*N 
[num_Sym,den_Sym] =numden(G_Sym);
den_LCM=den_Sym(1,1); 
for i=1:size(den_Sym,1)
    for j=1:size(den_Sym,2)
        den_LCM=lcm(den_LCM,den_Sym(i,j));
    end
end
for i=1:size(den_Sym,1)
    for j=1:size(den_Sym,2)
        num_LCM(i,j)=simplify((den_LCM/den_Sym(i,j))*num_Sym(i,j));
    end
end
%%	Creating Left and Right MFDs
DL=den_LCM*eye(size(num_Sym,1));	% Left MFD: DL
NL=num_LCM;	% Left MFD: NL
Deg_L=size(sym2poly(det(DL)),2)-1;	% Order of Right MFD
DR=den_LCM*eye(size(num_Sym,2));	% Right MFD: DR
NR=num_LCM;	% Right MFD: NR
Deg_R=size(sym2poly(det(DR)),2)-1;	% Order of Right MFD
%%	Minimality of MFDs
PL=[NL DL];
[Smith_MacMilan,P,Z]=SmithMacMilan(PL);	% Smith Formf of [NL DL]
Min_S=min(size(Smith_MacMilan));
Min_DL=size(sym2poly(det(Smith_MacMilan(1:Min_S,1:Min_S))),2)-1;
if Min_DL==0
Minimality_L='left MFD is MINIMAL';
else
Minimality_L='Left MFD is not MINIMAL';
end
PR=[NR ; DR];
[Smith_MacMilan,P,Z]=SmithMacMilan(PR);	% Smith Formf of [NR;DR]
Min_S=min(size(Smith_MacMilan));
Min_DR=size(sym2poly(det(Smith_MacMilan(1:Min_S,1:Min_S))),2)-1;
if Min_DR==0
Minimality_R='Right MFD is MINIMAL';
else
Minimality_R='Right MFD is not MINIMAL';
end
Minimality_R=Minimality_R;
Minimality_L=Minimality_L;
%%	The End of Program.
end

