clc
clear all
close all
syms s
g11=-1/(s+1)
g12=1/(s+1)
g21=-1/(s+1)
g22=-1/(s+1)
g=[g11 g12;
   g21 g22]
g=simplifyFraction(g)
eigs=eig(g)
pretty(eigs)
nyquist(sym2tf(simplify(eigs(1))))
hold on 
nyquist(sym2tf(simplify(eigs(2))))
%%  eig orthogonal
clc
clear all 
close all
A=(1/6)*[3 4;4.5 6]
[V,E]=eig(A,'vector')
[V,~]=qr(V)
[VV,EE]=eig(A','vector')
[VV,~]=qr(VV)
%% 
clc
clear all 
close all
A=(1/6)*[3 4;4.5 6]
[V,E]=eig(A,'vector')
[V,~]=qr(V)
[VV,EE]=eig(A','vector')
[VV,~]=qr(VV)
%%
A=(1/6)*[3 4;4.5 6]
[U,S,V] = svd(A)
%%
syms s
Q=s^2-2*s-16==0
zeroes=solve(Q,s)
%%
clc
clear all
close all
% G=[tf([1],[1 2]) tf([1 1],[1 0]);
%    tf([1 0],[1 4]) tf([2],[1 4])]
s=sym('s'); 
% G=input('G = '); 
g11=2/(s+1)
g12=s/(s+2)
g21=-s/(s+1)
g22=3/(s+1)
G=[g11 g12;
   g21 g22]
[M,poles,zeros] = smform(G)