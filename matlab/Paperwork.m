%% State Space Model Predictive Control of a Reactive
clc
clear all
close all
format shortg
syms s
A=[-0.0019 0.0024 -0.0009;-0.0034 -0.0043 -0.0008;-0.0009 -0.0052 -0.0055];
B=[0.0003 0.0013 -0.0566; -0.0005 0.0019 -0.0822; 0.0025 0.0067 -0.2818];
C=[-12.8258 4.0787 0.1685;-3.7520 -0.8266 -1.9955;-6.2674 -17.6857 -0.7129];
D=zeros(3,3);
K=[-0.0440 -0.0058 -0.0114;0.0172 0.0046 -0.0430;0.0147 -0.1111 -0.0099]
SYSTEM=ss(A,B,C,D)
%% 1

%% 2

%% 3
A 
B
C
D
%% 4
%% 5
%% 6
state_system=ss(A,B,C,D)
Top_Input=4*[zeros(100,1);ones(4900,1)];                                   
Reactive_Input=1*[zeros(100,1);ones(4900,1)];
Bottom_Input=0.56*[zeros(100,1);ones(4900,1)];                      
Input=[Top_Input Reactive_Input Bottom_Input];                                      
Time=(0:1:4999)';                                             
Output_state_space=lsim(state_system,Input,Time);

figure('Name','state-space Output','NumberTitle','off');
subplot(3,1,1);
plot(Time,Output_state_space(:,1));
p_x=xlabel('time(sec)');p_y=ylabel('Top Temp(C)');
set(p_y,'Color','k','fontsize',10,'Fontweight','b'); 
hold on
plot(Time,Top_Input);
hold off
title('Top Segment')

subplot(3,1,2);
plot(Time,Output_state_space(:,2));
p_x=xlabel('time(sec)');p_y=ylabel('Reative Temp(C)');
set(p_y,'Color','k','fontsize',10,'Fontweight','b'); 
hold on
plot(Time,Reactive_Input);
hold off
title('Reactive')

subplot(3,1,3);
plot(Time,Output_state_space(:,3));
p_x=xlabel('time(sec)');p_y=ylabel('Bottom Temp(C)');
set(p_y,'Color','k','fontsize',10,'Fontweight','b'); 
hold on
plot(Time,Bottom_Input);
hold off
title('Bottom segment')

%% 7
Transfer_function=C*(inv(s*eye(3 )-A))*B+D;
%Transfer_function=vpa(simplifyFraction(vpa(Transfer_function)))
Rosenbrock=[s*eye(3)-A B;-C D]
[G_Sym,DL,NL,Deg_L,Minimality_L,DR,NR,Deg_R,Minimality_R]=MFD_Form(Transfer_function)

%% 8
figure('Name','Transfer-function Output','NumberTitle','off');
Transfer_function_laplace=sym2tf(Transfer_function);
Output_transfer_function=lsim(ss(Transfer_function_laplace),Input,Time);
subplot(3,1,1);
plot(1:size(Output_transfer_function(:,1),1),Output_transfer_function(:,1));
p_x=xlabel('time(sec)');p_y=ylabel('Top Temp(C)');
set(p_y,'Color','k','fontsize',10,'Fontweight','b'); 
title('Top Segment')

subplot(3,1,2);
plot(1:size(Output_transfer_function(:,2),1),Output_transfer_function(:,2));
p_x=xlabel('time(sec)');p_y=ylabel('Reative Temp(C)');
set(p_y,'Color','k','fontsize',10,'Fontweight','b'); 
title('Reactive')

subplot(3,1,3);
plot(1:size(Output_transfer_function(:,3),1),Output_transfer_function(:,3));
p_x=xlabel('time(sec)');p_y=ylabel('Bottom Temp(C)');
set(p_y,'Color','k','fontsize',10,'Fontweight','b'); 
title('Bottom segment')
% figure,
% step(sym2tf(Transfer_function))

%% 9
figure('Name','diagonal controller Transfer-function Output','NumberTitle','off');
Transfer_function_tff=sym2tf(Transfer_function);
con11=(1.2/1)*(1/1.23)*tf([1 1],[2.5 10.8])
con22=(1.05/1.28)*tf([-1 -1],[2 2.8])
con33=(1/1.75)*tf([75 1],[3 720])
controller=blkdiag(con11,con22,con33)
%step(Transfer_function_tff*controller)
Top_Input=[zeros(100,1);ones(4900,1)];                                   
Reactive_Input=[zeros(300,1);ones(4700,1)];
Bottom_Input=[zeros(500,1);ones(4500,1)];                      
Input=[Top_Input Reactive_Input Bottom_Input];                                      
Time=(0:1:4999)';  
lsim(Transfer_function_tff*controller,Input,Time);
%% 10
% [Transfer_function, Smith_form, multiplicity_poles, multiplicity_zeros, Poles, Zeros] = SmithMac_Form(Transfer_function);
% System_degree=size(sym2poly(P),2)-1
Deg_L
Minimality_L
Deg_R
Minimality_R
%% 11
[ coprime_numinator_right, coprime_dominator_right ] = minRMFD( Transfer_function )
%% 12
disp('for state space we have:')
s*eye(3)-A
pole_polynomial=det(s*eye(3)-A)
State_poles=solve(pole_polynomial==0,s);
State_poles=double(State_poles)
disp('for transfer function we have:')
Transfer_function
[M,Poles,Zeros,multiplicity_zeros,multiplicity_poles] = smform(Transfer_function)
Poles
%% 13
multiplicity_poles
multiplicity_zeros
%% 14
disp('calculate type 0 of G(s)')
delta(s)=multiplicity_poles
G0(s)=Transfer_function
delta(0)*det(G0(0))
disp('calculate type 1 of G(s) we analyze sG(s)')
G1=s.*Transfer_function
[Smith_form11,Poles11,Zeros11,multiplicity_zeros11,multiplicity_poles11] = smform(G1)
delta(s)=multiplicity_poles11
G0(s)=G1
delta(0)*det(G0(0))
disp('cause the value get zero we have singularity so the type of system will be zero')
disp('system is vector type [0 0 0]')
%% 15
%element zeros
disp('element zeros:')
[num11,den11]=numden(Transfer_function(1,1));
element11_zeros=double(solve(num11==0,s))
[num12,den12]=numden(Transfer_function(1,2));
element12_zeros=double(solve(num12==0,s))
[num13,den13]=numden(Transfer_function(1,3));
element13_zeros=double(solve(num13==0,s))
[num21,den21]=numden(Transfer_function(2,1));
element21_zeros=double(solve(num21==0,s))
[num22,den22]=numden(Transfer_function(2,2));
element22_zeros=double(solve(num22==0,s))
[num23,den23]=numden(Transfer_function(2,3));
element23_zeros=double(solve(num23==0,s))
[num31,den31]=numden(Transfer_function(3,1));
element31_zeros=double(solve(num31==0,s))
[num32,den32]=numden(Transfer_function(3,2));
element32_zeros=double(solve(num32==0,s))
[num33,den33]=numden(Transfer_function(3,3));
element33_zeros=double(solve(num33==0,s))
disp('transmission zeros:(cause we have D=0) we dont have transmission zeros')
Zeros
disp('for inut decoupling zero we should analyze [s*I-A B]')
Input_transfer_function=Rosenbrock(1:3,:)
normal_rank_input=rank(Input_transfer_function)
Output_transfer_function=Rosenbrock(:,1:3)
normal_rank_output=rank(Output_transfer_function)
Invariant_transfer_function=Rosenbrock
normal_rank_invariant=rank(Invariant_transfer_function)
disp('Matrix B & C are numbers and dont relative to frequency')
disp('so our normal rank wont change in different frequencies so we dont have input decoupling zeros')
disp('and output decouling zeros and invariant zeros')
%% 16
disp('we dont have transmission zero so we dont have zero direction')
%% 17
disp('we want to place zeros at -2 for example')
Acl=-2*eye(3)
%Acl=A-B*Q*C
K_controller=inv(B)*(A-Acl)*inv(C)
zero_placement_function=inv(K_controller)+Transfer_function;
[zero_placement_Smith_form11,zero_placement_Poles11,zero_placement_Zeros11,zero_placement_multiplicity_zeros11,zero_placement_multiplicity_poles11] = smform(zero_placement_function)
%% 18
disp('<strong>Controllabiliy</strong>')
disp('1.The 𝒏×𝒏𝒎 controllability matrix 𝚽=[𝑩 𝑨𝑩 𝑨^𝟐𝑩] has rank')
phi_C=[B A*B A*A*B]
rank_phi_C=rank(phi_C)
disp('2.matrix W is nonsingular for ant t>0')
Wc = gram(SYSTEM,'c')
rank_control_grammian=rank(Wc)
disp('For every eigenvalue 𝝀 of 𝑨, the 𝒏×𝒏+𝒎 complex matrix 𝝀𝑰−𝑨𝑩 has rank 𝒏(𝒇𝒖𝒍𝒍 𝒓𝒐𝒘 𝒓𝒂𝒏𝒌)')
lambda=eig(A)
% for eig1
control_part_3_1=[lambda(1)*eye(3)-A B]
control_part_3_1_rank=rank(control_part_3_1)
control_part_3_2=[lambda(2)*eye(3)-A B]
control_part_3_2_rank=rank(control_part_3_2)
control_part_3_3=[lambda(3)*eye(3)-A B]
control_part_3_3_rank=rank(control_part_3_3)
disp('All eigenvalues of 𝑨 have negative real parts, and the unique solution 𝑾𝒄 is positive definite.')
A_eigenvalues=eig(A)
Wc_eigenvalues=eig(Wc)

disp('<strong>Observability</strong>')
disp('1.The 𝒏×𝒏𝒎 controllability matrix 𝚽=[C^t A^tC^t A^t^2C^2] has rank')
phi_O=[C' A'*C' A'*A'*C']
rank_phi_C=rank(phi_O)
disp('2.matrix W is nonsingular for ant t>0')
Wo = gram(SYSTEM,'o')
rank_control_grammian=rank(Wc)
disp('For every eigenvalue 𝝀 of 𝑨, the 𝒏×𝒏+𝒎 complex matrix 𝝀𝑰−𝑨𝑩 has rank 𝒏(𝒇𝒖𝒍𝒍 𝒓𝒐𝒘 𝒓𝒂𝒏𝒌)')
lambda=eig(A)
% for eig1
observe_part_3_1=[lambda(1)*eye(3)-A;C]
observe_part_3_1_rank=rank(observe_part_3_1)
observe_part_3_2=[lambda(2)*eye(3)-A;C]
observe_part_3_2_rank=rank(observe_part_3_2)
observe_part_3_3=[lambda(3)*eye(3)-A;C]
observe_part_3_3_rank=rank(observe_part_3_3)
disp('All eigenvalues of 𝑨 have negative real parts, and the unique solution 𝑾𝒄 is positive definite.')
A_eigenvalues=eig(A)
Wc_eigenvalues=eig(Wo)
%% 19
[Controllability_Mu,Controllability_Index,Partial_Controllability_Matix]=Controllability_test(A,B) 
[Observibility_Mu,Observibility_Index,Partial_Observibility_Matix]=Controllability_test(A',C') 
%% 20
disp('system is observable and controllable')
%% 21
phi_Out=[C*B C*A*B C*A*A*B  D]
rank_phi_Out=rank(phi_Out)
disp('for functional controllability if det(G(s))!=0 then matrix is functionally controllable')
det(Transfer_function)
%% 22
%% 23
[A_direct_realize_control,B_direct_realize_control,C_direct_realize_control,D_direct_realize_control]=Realization_Controllable(Transfer_function)
[A_direct_realize_observe,B_direct_realize_observe,C_direct_realize_observe,D_direct_realize_observe]=Realization_Observable(Transfer_function)
%% 24
Hankel_0=D;
Hankel_1=C*B;
Hankel_2=C*A^1*B;
Hankel_3=C*A^2*B;
Hankel_4=C*A^3*B;
Hankel_5=C*A^4*B;
Hankel_matrix=[Hankel_1 Hankel_2 Hankel_3;Hankel_2 Hankel_3 Hankel_4;Hankel_3 Hankel_4 Hankel_5]
[Hankel_input,Hankel_singularvalue,Hankel_output]=svd(Hankel_matrix)
Hankel_input_hat=Hankel_input*(Hankel_singularvalue.^(1/2))
Hankel_output_hat=(Hankel_singularvalue.^(1/2))*Hankel_output
disp('observable canonical form by hankel matrix:')
det(coprime_dominator_right)
alpha1=0.011726;
alpha2=5.8006e-05;
alpha3=9.4646e-08;
A_observe_hankel=[zeros(3,3) eye(3) zeros(3,3);zeros(3,3) zeros(3,3) eye(3);-alpha3.*eye(3) -alpha2.*eye(3) -alpha1.*eye(3)]
B_observe_hankel=[Hankel_1;Hankel_2;Hankel_3]
C_observe_hankel=[eye(3) zeros(3,3) zeros(3,3)]
D_observe_hankel=Hankel_0;
disp('controllable canonical form by hankel matrix:')
A_control_hankel=[zeros(3,3) zeros(3,3) -alpha3.*eye(3);eye(3) zeros(3,3) -alpha2.*eye(3);zeros(3,3) eye(3) -alpha1.*eye(3)]
B_control_hankel=[eye(3) zeros(3,3) zeros(3,3)]'
C_control_hankel=[Hankel_1;Hankel_2;Hankel_3]'
D_control_hankel=Hankel_0
%% 25
figure('Name','Transfer-function Singular values','NumberTitle','off');
%s=tf('s');
Transfer_function_tf=sym2tf(Transfer_function)
% Define the Frequency
w=logspace(-3,3,500);
% Singular Values
sv=sigma(Transfer_function_tf,w);
% Plotting the Singular Values
loglog(w,sv(1,:),'linewidth',1.2)
hold on
loglog(w,sv(2,:),'linewidth',1.2)
hold on
loglog(w,sv(3,:),'linewidth',1.2)
xlabel('\omega [rad/s]')
ylabel('Singular Values of G(s)')
grid
%% 26
Gilbert_realization=canon(ss(A,B,C,D),'modal')
%% 27
Truncation=Order_Reduction(A,B,C,D,2,'Truncation')
Residualization=Order_Reduction(A,B,C,D,2,'Residualization')
%% 28
Balanced_Truncation=Order_Reduction(A,B,C,D,2,'Balanced Truncation')
Balanced_Residualization=Order_Reduction(A,B,C,D,2,'Balanced Residualization')
%% 29
% by looking at transfer function we see that our relative degree of every
% row is 1 so we have D(s) as follows
Relative_degree=s*eye(3)
B_star=limit(Relative_degree*Transfer_function,s,inf)
Rank_B_Star=rank(B_star)
disp('Bstar matrix is singular so we cant have static state feedback')
% so we have state feedback decouplation
 B_bar=[C(1,:)*A;
        C(2,:)*A;
        C(3,:)*A];
A_new_decouple=A-B*inv(B_star)*B_bar
B_new_decouple=B*inv(B_star)
C_new_decouple=C-D*inv(B_star)*B_bar
D_new_decouple=D*inv(B_star)

% syms s
% Dd=[s^-d1 0;
%     0 s^-d2]
%% 30
K=eye(3);
A_CL=A-B*inv((inv(K)+D))*C
poles_for_BIBO_stability=eig(A_CL)
K=-eye(3);
A_CL=A-B*inv((inv(K)+D))*C
poles_for_BIBO_stability=eig(A_CL)

K=-1/250*eye(3);
A_CL=A-B*inv((inv(K)+D))*C
poles_for_BIBO_stability=eig(A_CL)
%% 31
k=eye(3);
H12_transfer=inv(eye(3)-k*Transfer_function)*Transfer_function
[Smith_form_H12,Poles_H12,Zeros_H12,multiplicity_zeros_H12,multiplicity_poles_H12] = smform(H12_transfer)
k1=-eye(3);
H12_transfer_1=inv(eye(3)-k1*Transfer_function)*Transfer_function
[Smith_form_H12_1,Poles_H12_1,Zeros_H12_1,multiplicity_zeros_H12_1,multiplicity_poles_H12_1] = smform(H12_transfer_1)
k2=-0.001.*eye(3);
H12_transfer_2=inv(eye(3)-k2*Transfer_function)*Transfer_function
[Smith_form_H12_2,Poles_H12_2,Zeros_H12_2,multiplicity_zeros_H12_2,multiplicity_poles_H12_2] = smform(H12_transfer_2)
%% 32
nyquist_plotter(Transfer_function)
%% 33
nyquist_plotter(pinv(Transfer_function))
%% 6.1
Disturbance_Transfer_function=C*(inv(s*eye(3 )-A))*K+eye(3)
Disturbance_Transfer_function=simplifyFraction(Disturbance_Transfer_function)
%% 6.2
Transfer_function_scaled=[Transfer_function(:,1)/9 Transfer_function(:,2)/5 Transfer_function(:,3)/0.6]
Sensitivity_function=inv(eye(3)+Transfer_function_scaled)
complementary_Sensitivity_function=Transfer_function*inv(eye(3)+Transfer_function_scaled)
figure('Name','Sensitivity Transfer-function Singular values','NumberTitle','off');
%s=tf('s');
Sensitivity_function_tf=sym2tf(Sensitivity_function)
% Define the Frequency
w=logspace(-3,3,500);
% Singular Values
sv=sigma(Sensitivity_function_tf,w);
% Plotting the Singular Values
loglog(w,sv(1,:),'linewidth',1.2)
hold on
loglog(w,sv(2,:),'linewidth',1.2)
hold on
loglog(w,sv(3,:),'linewidth',1.2)
xlabel('\omega [rad/s]')
ylabel('Singular Values of G(s)')
grid
figure('Name','Complementary Sensitivity Transfer-function Singular values','NumberTitle','off');
%s=tf('s');
complementary_Sensitivity_function_tf=sym2tf(complementary_Sensitivity_function)
% Define the Frequency
w=logspace(-3,3,500);
% Singular Values
sv=sigma(complementary_Sensitivity_function_tf,w);
% Plotting the Singular Values
loglog(w,sv(1,:),'linewidth',1.2)
hold on
loglog(w,sv(2,:),'linewidth',1.2)
hold on
loglog(w,sv(3,:),'linewidth',1.2)
xlabel('\omega [rad/s]')
ylabel('Singular Values of G(s)')
grid
norm(complementary_Sensitivity_function_tf,inf)
norm(Sensitivity_function_tf,inf)
%% 6.3
Functionally_control_ch6=rank(Transfer_function)
%% 6.4
disp('we dont have time delay to discuss this part')
%% 6.5
disp('we dont have unstable zero to discuss this part')
%% 6.6
disp('we dont have unstable pole to discuss this part')
%% 6.7
disp('we dont have unstable pole neither unstable zero this part')
%% 6.8
figure('Name','elements of Gd','NumberTitle','off')
bodemag(sym2tf(Disturbance_Transfer_function),logspace(-3,3,120))
xlabel('\omega [rad/s]')
ylabel('elements of Gd')
direction_of_first_disturbance=(1/norm(Disturbance_Transfer_function(:,1))*Disturbance_Transfer_function(:,1))
direction_of_second_disturbance=(1/norm(Disturbance_Transfer_function(:,2)))*Disturbance_Transfer_function(:,2)
direction_of_third_disturbance=(1/norm(Disturbance_Transfer_function(:,3)))*Disturbance_Transfer_function(:,3)
distt1=Sensitivity_function*Disturbance_Transfer_function(:,1)
value_distt1=norm(sym2tf(distt1),inf)
distt2=Sensitivity_function*Disturbance_Transfer_function(:,2)
value_distt2=norm(sym2tf(distt2),inf)
distt3=Sensitivity_function*Disturbance_Transfer_function(:,3)
value_distt3=norm(sym2tf(distt3),inf)


disp('and we dont have zero direction')
%% 6.9
figure('Name','elements of inv(G)Gd','NumberTitle','off')
bodemag(sym2tf(Disturbance_Transfer_function*Transfer_function_scaled),logspace(-3,3,120))
xlabel('\omega [rad/s]')
ylabel('element of inv(G)Gd')
%% 6.10 in chapter6 script
%% 6.11
singular_value_dis1_direction(s)=svd(Sensitivity_function*direction_of_first_disturbance)
singular_value_dis2_direction(s)=svd(Sensitivity_function*direction_of_second_disturbance)
singular_value_dis3_direction(s)=svd(Sensitivity_function*direction_of_third_disturbance)
max_singular_value_transfer=svd(Transfer_function)
max_singular_value_distrubance_transfer_1=svd(pinv(Transfer_function)*direction_of_first_disturbance)
max_singular_value_distrubance_transfer_2=svd(pinv(Transfer_function)*direction_of_second_disturbance)
max_singular_value_distrubance_transfer_3=svd(pinv(Transfer_function)*direction_of_third_disturbance)
Disturbance_condition_number_1(s)=max_singular_value_transfer(1)*max_singular_value_distrubance_transfer_1(1)
Disturbance_condition_number_2(s)=max_singular_value_transfer(1)*max_singular_value_distrubance_transfer_2(1)
Disturbance_condition_number_3(s)=max_singular_value_transfer(1)*max_singular_value_distrubance_transfer_3(1)
Disturbance_condition_number_1(0)
Disturbance_condition_number_2(0)
Disturbance_condition_number_3(0)

Transfer_function_s(s)=Transfer_function
RGA_steady=Transfer_function_s(0).*pinv(Transfer_function_s(0)).'
svddd=svd(Transfer_function_s(0))
Condition_number=svddd(1)/svddd(3)
% figure('Name','elements RGA','NumberTitle','off')
% bodemag(sym2tf(Transfer_function_s.*pinv(Transfer_function_s.')),logspace(-3,3,120))
% xlabel('\omega [rad/s]')
%% 8.1 in chapter8 script
%% 8.2 in chapter8 script
%% 8.3 in chapter8 script
%% 8.4 in chapter8 script
%% 8.5 in chapter8 script
%% 8.6 in chapter8 script
%% 8.7 in chapter8 script
%% 8.8 in chapter8 script
%% 8.9 in chapter8 script
%% 8.10 in chapter8 script
%% 8.11 in chapter8 script
%% 8.12 in chapter8 script

%% author comment
% vpa(G) for simplify
%%                                      End of program
% Author: Javad Hosseinzadeh
% email: javadhosseinzadeh701@gmail.com 
% Created: February 2021