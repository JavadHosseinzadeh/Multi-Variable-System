function [ NR, DR ] = minRMFD( gs )
%Copyright (C) year 2012   Ahmadreza Saadatkhah  <ar.saadatkhah@gmail.com>
% it's a function that calculate a minimal right matrix fraction
% description for any MIMO transfer function.
%the input must be symbolic transfer function
%Revised on October 2019
syms s;
[~,d]=numden(gs);
[l,m]=size(gs);
L=1;
for i=1:l
    for j=1:m
        L=lcm(d(i,j),L);
    end
end
DL=L*eye(l);
NL=simplify(DL*gs);
M=length(sym2poly(NL(1,1)));
for i=1:l
    for j=1:m
        M=max(length(sym2poly(NL(i,j))),M);
    end
end
sMatsz=max(length(sym2poly(L)),M);
fn=zeros(l,m,sMatsz);
fd=zeros(l,l,sMatsz);
for r=1:l
    for c=1:m
        aux=sym2poly(NL(r,c));
        fn(r,c,1:(sMatsz-length(aux))) = 0;
        fn(r,c,(sMatsz-length(aux)+1):sMatsz)=aux;
    end
end
for r=1:l
    for c=1:l
        aux=sym2poly(DL(r,c));
        fd(r,c,1:(sMatsz-length(aux))) = 0;
        fd(r,c,(sMatsz-length(aux)+1):sMatsz)=aux;
    end
end
sMat=zeros(2*l*sMatsz-l,(l+m)*sMatsz);
for i=1:l:l*sMatsz
    for j=1:(l+m):(l+m)*(((i-1)/l)+1)
        sMat(i:i+l-1,j:j+l+m-1)=[fd(:,:,(sMatsz-((i-1)/l))+(j-1)/(l+m)) fn(:,:,(sMatsz-((i-1)/l))+(j-1)/(l+m))];
    end
end
for i=(l*sMatsz+1):l:2*l*(sMatsz)-l
    for j=(l+m)*((i-1-l*sMatsz)/l +1)+1:(l+m):(l+m)*sMatsz
        sMat(i:i+l-1,j:j+l+m-1)=[fd(:,:,(sMatsz-((i-1)/l))+(j-1)/(l+m)) fn(:,:,(sMatsz-((i-1)/l))+(j-1)/(l+m))];
    end
end
z=null(sMat,'r');
[r,c]=size(z);
flag=1;
while (flag~=0)
    NR=zeros(l,m);
    DR=zeros(m,m);
    for i=1:l+m:r+1-l-m
        NR=-z(i:i+l-1,1:m)*(s^((i-1)/(l+m)))+NR;
        DR=z(i+l:i+l+m-1,1:m)*(s^((i-1)/(l+m)))+DR;
    end
    if det(DR)==0
        if flag>(c-m)
            disp('Oops! it seems that I could not find a non-singular polynomial matrix for a minimal DR. Please report to me (ar.saadatkhah@gmail.com) if this message displays, so I will find a soultion for this problem, thanks.')
            flag=0;
        else
            z(:,m)=[];
            flag=flag+1;
        end
    else
        flag=0;        
    end
end
NR=simplify(NR);
DR=simplify(DR);
end