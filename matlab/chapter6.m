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

systemnames = 'Gn Gd'; 
inputvar = '[r{3}; d{3}]';
outputvar = '[Gn+Gd]';
input_to_Gn = '[r]';
input_to_Gd = '[d]';
cleanupsysic = 'yes';
T = sysic;
%% adame ghateiate zarbi vorudi bedune feedback delta .1
clc
Delta = ultidyn('Delta',[3 3],'Bound',.1);
Gi = Gn*(eye(3)+Delta);
systemnames = 'Gi Gd'; 
inputvar = '[r{3}; d{3}]';
outputvar = '[Gi+Gd]';
input_to_Gi = '[r]';
input_to_Gd = '[d]';
cleanupsysic = 'yes';
Ti_1 = sysic;
figure,
step(Ti_1)
hold on
step(T)
figure,
bodemag(Ti_1)
hold on
bodemag(T)
%% delta 1
clc
Delta = ultidyn('Delta',[3 3],'Bound',1);
Gi = Gn*(eye(3)+Delta);
systemnames = 'Gi Gd'; 
inputvar = '[r{3}; d{3}]';
outputvar = '[Gi+Gd]';
input_to_Gi = '[r]';
input_to_Gd = '[d]';
cleanupsysic = 'yes';
Ti = sysic;
figure,
step(Ti)
hold on
step(T)
figure,
bodemag(Ti)
hold on
bodemag(T)
%% adame ghateiate zarbi khoruji bedune feedback delta .1
Delta = ultidyn('Delta',[3 3],'Bound',.1);
Go = (eye(3)+Delta)*Gn;
systemnames = 'Go Gd'; 
inputvar = '[r{3}; d{3}]';
outputvar = '[Go+Gd]';
input_to_Go = '[r]';
input_to_Gd = '[d]';
cleanupsysic = 'yes';
To_1 = sysic;
figure,
step(To_1)
hold on
step(T)
figure,
bodemag(To_1)
hold on
bodemag(T)
%% delta 1
Delta = ultidyn('Delta',[3 3],'Bound',.1);
Go = (eye(3)+Delta)*Gn;
systemnames = 'Go Gd'; 
inputvar = '[r{3}; d{3}]';
outputvar = '[Go+Gd]';
input_to_Go = '[r]';
input_to_Gd = '[d]';
cleanupsysic = 'yes';
To = sysic;
figure,
step(To)
hold on
step(T)
figure,
bodemag(To)
hold on
bodemag(T)
%%                                      End of program
% Author: Javad Hosseinzadeh
% email: javadhosseinzadeh701@gmail.com 
% Created: February 2021