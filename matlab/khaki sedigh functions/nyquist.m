clc
clear all
close all
% G=[tf([1],[1 2]) tf([1 1],[1 0]);
%    tf([1 0],[1 4]) tf([2],[1 4])]
s=sym('s'); 
%G=input('G = '); 
s=sym('s'); 
% G=input('G = '); 
g11=0.5*(s+2)/(s^2+0.6*s+1)
g12=4*(s+0.1)/(s^2+0.6*s+1)
g21=.3*(s+8.3)/(s^2+0.6*s+1)
g22=3.2*(s+1.25)/(s^2+0.6*s+1)
G=[g11 g12;
   g21 g22]
% G=[1/(s+2) (s+1)/s;
%     s/(s+4) 2/(s+4)]
switch class(G) 
    case 'tf'                                                              % Transfer Function Form 
        G_TF=G; [num,den]=tfdata(G);                                                
         for i=1:size(num,1)                                                                     
            for j=1:size(num,2)                                                                  
                num_symbolic(i,j)=poly2sym(num{i,j},'s');
                den_symbolic(i,j)=poly2sym(den{i,j},'s');  
            end 
        end 
        G_S=simple(num_symbolic./den_symbolic);                            % Converting to Symbolic Form 
    case 'sym'                                                             % Symbolic Form 
        G_S=simplifyFraction(G); [a,b]=numden(G_S); 
        for i=1:size(a,1)                                                                       
            for j=1:size(a,2)                                                                   
                num{i,j}=sym2poly(a(i,j)); den{i,j}=sym2poly(b(i,j)); 
                G_TF(i,j)=tf( num{i,j}, den{i,j})                          % Converting to Transfer Function Form 
            end 
        end 
end 
g=G_TF; 
 
%%                                                                  Greshgorin Bands 
[n,m]=size(g);                                                             % Number of Inputs and Outputs 
w=logspace(-5,6,150);                                                      % Frequency Range 
q=0:(pi/100):(2*pi); 
figure 
for i=1:n 
    for j=1:m 
        if i==j 
            subplot(n,m,(i-1)*size(num,2)+j); nyquist(g(i,i));             % Nyquist Diagram of G(i,i) 
            title(['Nyquist Diagram of G(',num2str(i),',',num2str(j),')']) 
%%                                                                Row Diagonal Dominant 
            for iest=i:i                                                                                    
                R=zeros(1,length(w)); 
                for jest=1:m 
                    if iest~=jest 
                        hold on 
                        C=subs(G_S(i,j),complex(0,w));                     % Center of Circles 
                        R=R+abs(subs(G_S(iest,jest),complex(0,w)));        % Radius Circles 
                        if jest==m 
                            for k=1:length(C) 
                                plot((R(k)*cos(q))+real(C(k)),(R(k)*sin(q))+imag(C(k)),'r-')  
                            end 
                        end 
                        if (iest==m) && (jest==m-1) 
                            for k=1:length(C) 
                                plot((R(k)*cos(q))+real(C(k)),(R(k)*sin(q))+imag(C(k)),'r-')   
                            end 
                        end 
                        hold off 
                    end 
                end 
            end 
        else 
            subplot(n,m,(i-1)*size(num,2)+j); nyquist(g(i,j));             % Nyquist Diagram of G(i,j) 
            title(['Nyquist Diagram of G(',num2str(i),',',num2str(j),')']) 
        end 
    end 
end 
%%                                                               The End of Program. 
% Author: Peyman Bagheri 
% email: peyman.bk@gmail.com 
% Created: September 2011;  
% Last revision: 09-November-2011; 

