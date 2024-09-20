clear all
clc
close all
%% unstructured uncertainty
%% taarife system ha va adame ghateiat
A=[-0.0019 0.0024 -0.0009;-0.0034 -0.0043 -0.0008;-0.0009 -0.0052 -0.0055];
B=[0.0003 0.0013 -0.0566; -0.0005 0.0019 -0.0822; 0.0025 0.0067 -0.2818];
C=[-12.8258 4.0787 0.1685;-3.7520 -0.8266 -1.9955;-6.2674 -17.6857 -0.7129];
D=zeros(3,3);
K=[-0.0440 -0.0058 -0.0114;0.0172 0.0046 -0.0430;0.0147 -0.1111 -0.0099];
syms s
Gns = (C*inv(s*eye(3)-A)*B+D);
Gn = sym2tf(Gns);
Gds = (C*inv(s*eye(3)-A)*K)+eye(3);
Gd = sym2tf(Gds);
Delta = ultidyn('Delta',[3 3],'Bound',1);
%% adame ghateiate jami bedune feedback
clc
Ga = Gn+Delta;
systemnames = 'Ga Gd'; 
inputvar = '[r{3}; d{3}]';
outputvar = '[Ga+Gd]';
input_to_Ga = '[r]';
input_to_Gd = '[d]';
cleanupsysic = 'yes';
Ta = sysic;
[stabmarg,destabunc,report] = robuststab(Ta)
Ta.NominalValue
[Ma,Delta] = lftdata(Ta);
Ma;
[perfmarg,perfmargunc,report,info] = robustperf(Ta);
step(Ta)
%% adame ghateiate jami ba feedback
clc
Gaf = Gn+Delta;
systemnames = 'Gaf Gd'; 
inputvar = '[r{3}; d{3}]';
outputvar = '[r-Gaf-Gd]';
input_to_Gaf = '[r-Gaf-Gd]';
input_to_Gd = '[d]';
cleanupsysic = 'yes';
Taf = sysic;
[stabmarg,destabunc,report] = robuststab(Taf)
Taf.NominalValue
[Maf,Delta] = lftdata(Taf);
Maf;
[perfmarg,perfmargunc,report,info] = robustperf(Taf);
step(Taf)
%% adame ghateiate zarbi vorudi bedune feedback
clc
Gi = Gn*(eye(3)+Delta);
systemnames = 'Gi Gd'; 
inputvar = '[r{3}; d{3}]';
outputvar = '[Gi+Gd]';
input_to_Gi = '[r]';
input_to_Gd = '[d]';
cleanupsysic = 'yes';
Ti = sysic;
[stabmarg,destabunc,report] = robuststab(Ti)
Ti.NominalValue
[Mi,Delta] = lftdata(Ti);
Mi;
[perfmarg,perfmargunc,report,info] = robustperf(Ti);
step(Ti)
%% adame ghateiate zarbi vorudi ba feedback
clc
Gif = Gn*(eye(3)+Delta);
systemnames = 'Gif Gd'; 
inputvar = '[r{3}; d{3}]';
outputvar = '[r-Gif-Gd]';
input_to_Gif = '[r-Gif-Gd]';
input_to_Gd = '[d]';
cleanupsysic = 'yes';
Tif = sysic;
[stabmarg,destabunc,report] = robuststab(Tif)
Tif.NominalValue
[Mif,Delta] = lftdata(Tif);
Mif;
%% check kardane karaei moghavem baraye in halat (8.11)
[perfmarg,perfmargunc,report,info] = robustperf(Tif);
step(Tif)
%% adame ghateiate zarbi khoruji bedune feedback
Go = (eye(3)+Delta)*Gn;
systemnames = 'Go Gd'; 
inputvar = '[r{3}; d{3}]';
outputvar = '[Go+Gd]';
input_to_Go = '[r]';
input_to_Gd = '[d]';
cleanupsysic = 'yes';
To = sysic;
[stabmarg,destabunc,report] = robuststab(To)
To.NominalValue
[Mo,Delta] = lftdata(To);
Mo;
[perfmarg,perfmargunc,report,info] = robustperf(To);
step(To)
%% adame ghateiate zarbi khoruji ba feedback
Gof = (eye(3)+Delta)*Gn;
systemnames = 'Gof Gd'; 
inputvar = '[r{3}; d{3}]';
outputvar = '[r-Gof-Gd]';
input_to_Gof = '[r-Gof-Gd]';
input_to_Gd = '[d]';
cleanupsysic = 'yes';
Tof = sysic;
[stabmarg,destabunc,report] = robuststab(Tof)
Tof.NominalValue
[Mof,Delta] = lftdata(Tof);
Mof;
[perfmarg,perfmargunc,report,info] = robustperf(Tof);
step(Tof)
%% adame ghateiate jami maakoos bedune feedback
clc
Gia = Gn*inv(eye(3)-Delta*Gn);
systemnames = 'Gia Gd'; 
inputvar = '[r{3}; d{3}]';
outputvar = '[Gia+Gd]';
input_to_Gia = '[r]';
input_to_Gd = '[d]';
cleanupsysic = 'yes';
Tia = sysic;
[stabmarg,destabunc,report] = robuststab(Tia)
Tia.NominalValue
[Mia,Delta] = lftdata(Tia);
Mia;
[perfmarg,perfmargunc,report,info] = robustperf(Tia);
step(Tia)
%% adame ghateiate jami maakoos ba feedback
clc
Giaf = Gn*inv(eye(3)-Delta*Gn);
systemnames = 'Giaf Gd'; 
inputvar = '[r{3}; d{3}]';
outputvar = '[r-Giaf-Gd]';
input_to_Giaf = '[r-Giaf-Gd]';
input_to_Gd = '[d]';
cleanupsysic = 'yes';
Tiaf = sysic;
[stabmarg,destabunc,report] = robuststab(Tiaf)
Tiaf.NominalValue
[Miaf,Delta] = lftdata(Tiaf);
Miaf;
[perfmarg,perfmargunc,report,info] = robustperf(Tiaf);
step(Tiaf)
%% adame ghateiate zarbi vorudi maakoos bedune feedback
clc
Gii = Gn*inv(eye(3)-Delta);
systemnames = 'Gii Gd'; 
inputvar = '[r{3}; d{3}]';
outputvar = '[Gii+Gd]';
input_to_Gii = '[r]';
input_to_Gd = '[d]';
cleanupsysic = 'yes';
Tii = sysic;
[stabmarg,destabunc,report] = robuststab(Tii)
Tii.NominalValue
[Mii,Delta] = lftdata(Tii);
Mii;
[perfmarg,perfmargunc,report,info] = robustperf(Tii);
step(Tii)
%% adame ghateiate zarbi vorudi maakoos ba feedback
clc
Giif = Gn*inv(eye(3)-Delta);
systemnames = 'Giif Gd'; 
inputvar = '[r{3}; d{3}]';
outputvar = '[r-Giif-Gd]';
input_to_Giif = '[r-Giif-Gd]';
input_to_Gd = '[d]';
cleanupsysic = 'yes';
Tiif = sysic;
[stabmarg,destabunc,report] = robuststab(Tiif)
Tiif.NominalValue
[Miif,Delta] = lftdata(Tiif);
Miif;
[perfmarg,perfmargunc,report,info] = robustperf(Tiif);
step(Tiif)
%% adame ghateiate zarbi khoruji maakoos bedune feedback
Gio = inv(eye(3)-Delta)*Gn;
systemnames = 'Gio Gd'; 
inputvar = '[r{3}; d{3}]';
outputvar = '[Gio+Gd]';
input_to_Gio = '[r]';
input_to_Gd = '[d]';
cleanupsysic = 'yes';
Tio = sysic;
[stabmarg,destabunc,report] = robuststab(Tio)
Tio.NominalValue
[Mio,Delta] = lftdata(Tio);
Mio;
[perfmarg,perfmargunc,report,info] = robustperf(Tio);
step(Tio)
%% adame ghateiate zarbi khoruji maakoos ba feedback
Giof = inv(eye(3)-Delta)*Gn;
systemnames = 'Giof Gd'; 
inputvar = '[r{3}; d{3}]';
outputvar = '[r-Giof-Gd]';
input_to_Giof = '[r-Giof-Gd]';
input_to_Gd = '[d]';
cleanupsysic = 'yes';
Tiof = sysic;
[stabmarg,destabunc,report] = robuststab(Tiof)
Tiof.NominalValue
[Miof,Delta] = lftdata(Tiof);
Miof;
[perfmarg,perfmargunc,report,info] = robustperf(Tiof);
step(Tiof)
%% bakhshe dovom
clear all
clc
close all
%% structured uncertainty (adame ghateiate parametri)
a11 = ureal('a11',-0.0019,'Percentage',10);
a22 = ureal('a22',-0.0043,'Percentage',10);
a33 = ureal('a33',-0.0055,'Percentage',10);
a23 = ureal('a23',-0.0008,'Percentage',10);
A=[a11 0.0024 -0.0009;-0.0034 a22 a23;-0.0009 -0.0052 a33];
B=[0.0003 0.0013 -0.0566; -0.0005 0.0019 -0.0822; 0.0025 0.0067 -0.2818];
C=[-12.8258 4.0787 0.1685;-3.7520 -0.8266 -1.9955;-6.2674 -17.6857 -0.7129];
D=zeros(3,3);
A0=[-0.0019 0.0024 -0.0009;-0.0034 -0.0043 -0.0008;-0.0009 -0.0052 -0.0055];
Gss = ss(A,B,C,D);
Gss0 = ss(A0,B,C,D);
K=[-0.0440 -0.0058 -0.0114;0.0172 0.0046 -0.0430;0.0147 -0.1111 -0.0099];
Dd = eye(3);
Gdss = ss(A0,K,C,Dd);
%% paydarie moghaveme systeme halghe baz
systemnames = 'Gss Gdss'; 
inputvar = '[r{3}; d{3}]';
outputvar = '[Gss+Gdss]';
input_to_Gss = '[r]';
input_to_Gdss = '[d]';
cleanupsysic = 'yes';
T1 = sysic;
[M,Delta] = lftdata(T1)
[stabmarg,destabunc,report] = robuststab(T1)
%% paydarie moghaveme systeme halghe baste
systemnames = 'Gss Gdss'; 
inputvar = '[r{3}; d{3}]';
outputvar = '[r-Gss-Gdss]';
input_to_Gss = '[r-Gss-Gdss]';
input_to_Gdss = '[d]';
cleanupsysic = 'yes';
T2 = sysic;
[stabmarg,destabunc,report] = robuststab(T2)
%% karaeie moghaveme system ba adame ghateiate sakhtar yafte
[perfmarg,perfmargunc,report,info] = robustperf(T1);
step(T1)
%% taarife system ba vazneye karaei baraye tarahie controller
s = tf('s');
w = 500*(s+100)/(s+0.1);
W = blkdiag(w,w,w);
systemnames = 'Gss Gdss W';
inputvar = '[r{3}; d{3}; c{3}]';
outputvar = '[W; Gdss; r-Gss-Gdss]';
input_to_Gss = '[c]';
input_to_W = '[r-Gss-Gdss]';
input_to_Gdss = '[d]';
T3 = sysic;
%% tarahie controllere mu
n = 5;
opt=dksynOptions('NumberofAutoIterations',n,'DisplayWhileAutoIter','on','AutoIterSmartTerminate','off');
[cont,clp,bnd] = dksyn(T3,3,3,opt)
get(cont)
%% paydari o karaeie systeme halghe baste ba controller
clp_ic = lft(T3,cont);
[stabmarg,destabu,report,info] = robuststab(clp_ic)
[perfmarg,perfmargunc,report,info] = robustperf(clp_ic)
%% chek kardane natije
clc
cont=balred(cont,5);
systemnames = 'Gss cont Gdss'; 
inputvar = '[r{3}; d{3}]';
outputvar = '[Gdss+Gss]';
input_to_Gss = '[cont]'; 
input_to_cont = '[r-Gdss-Gss]';
input_to_Gdss = '[d]'; 
cleanupsysic = 'yes';
T4 = sysic; 
[stabmarg,destabu,report,info] = robuststab(T4)
[perfmarg,perfmargunc,report,info] = robustperf(T4)
step(T4)
%% singular value
s=tf('s')
[M,Delta] = lftdata(T1)
M_transfer_function=M.c*(s*eye(6)-M.a)*M.b+M.d
% Define the Frequency
w=logspace(-3,3,500);
% Singular Values
sv=sigma(M_transfer_function,w);
% Plotting the Singular Values
loglog(w,sv(1,:),'linewidth',1.2)
hold on
loglog(w,sv(2,:),'linewidth',1.2)
hold on
loglog(w,sv(3,:),'linewidth',1.2)
hold on
loglog(w,sv(4,:),'linewidth',1.2)
hold on
loglog(w,sv(5,:),'linewidth',1.2)
hold on
loglog(w,sv(6,:),'linewidth',1.2)
xlabel('\omega [rad/s]')
ylabel('Singular Values of G(s)')
%% ssv
%  BlockStructure=[-1 0;-1 0;-1 0;-1 0]
%  bounds = mussv(M,BlockStructure )
%%                                      End of program
% Author: Javad Hosseinzadeh
% email: javadhosseinzadeh701@gmail.com 
% Created: February 2021