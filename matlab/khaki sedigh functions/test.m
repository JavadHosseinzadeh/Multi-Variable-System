Sys=ss(A,B,C,D);                                                           % Define State-Space system 
Red_Sys=Order_Reduction(A,B,C,D,5,'Truncation');                           % Reducing system by Truncation 
[y1,t1]=step(Sys(1,1),30); [y2,t2]=step(Sys(1,2),30);                      % Step responces 
[y3,t3]=step(Sys(2,1),30); [y4,t4]=step(Sys(2,2),30); 
[y5,t5]=step(Red_Sys(1,1),30); [y6,t6]=step(Red_Sys(1,2),30); 
[y7,t7]=step(Red_Sys(2,1),30); [y8,t8]=step(Red_Sys(2,2),30); 
figure(1);                                                                 % Plot step responces 
subplot(2,2,2);plot(t1,y1,'r','linewidth',2.5);hold on;plot(t5,y5,'b','linewidth',2.5) 
p_y=ylabel('y_1'); set(p_y,'Color','k','fontsize',10,'Fontweight','b') 
p_x=xlabel('time (sec)'); set(p_x,'Color','k','fontsize',10,'Fontweight','b') 
subplot(2,2,1);plot(t2,y2,'r','linewidth',2.5);hold on;plot(t6,y6,'b','linewidth',2.5) 
p_y=ylabel('y_1'); set(p_y,'Color','k','fontsize',10,'Fontweight','b') 
p_x=xlabel('time (sec)'); set(p_x,'Color','k','fontsize',10,'Fontweight','b') 
subplot(2,2,4);plot(t3,y3,'r','linewidth',2.5);hold on;plot(t7,y7,'b','linewidth',2.5) 
p_y=ylabel('y_2'); set(p_y,'Color','k','fontsize',10,'Fontweight','b') 
p_x=xlabel('time (sec)'); set(p_x,'Color','k','fontsize',10,'Fontweight','b') 
subplot(2,2,3);plot(t4,y4,'r','linewidth',2.5);hold on;plot(t8,y8,'b','linewidth',2.5) 
p_y=ylabel('y_2'); set(p_y,'Color','k','fontsize',10,'Fontweight','b') 
p_x=xlabel('time (sec)'); set(p_x,'Color','k','fontsize',10,'Fontweight','b') 
hold off 
[a1,b1]=sigma(Sys); [a2,b2]=sigma(Red_Sys);                                % Singular values 
figure(2)                                                                  % Plot singular values 
semilogx(b1,20*log10(a1(1,:)),'r','linewidth',2.5); hold on 
semilogx(b1,20*log10(a1(2,:)),'b','linewidth',2.5) 
semilogx(b2,20*log10(a2(1,:)),'r--','linewidth',2.5) 
semilogx(b2,20*log10(a2(2,:)),'b--','linewidth',2.5) 
p_x=xlabel('Frequency (rad\sec)'); 
set(p_x,'Color','k','fontsize',10,'Fontweight','b') 
p_y=ylabel('Singular Values (dB)'); 
set(p_y,'Color','k','fontsize',10,'Fontweight','b'); hold off 
norm(Sys-Red_Sys,inf)                                                      % Error    
%%                                                               The End of Program. 
% Author: Peyman Bagheri 
% email: peyman.bk@gmail.com 
% Created: September 2011;  
% Last revision: 09-November-2011; 
