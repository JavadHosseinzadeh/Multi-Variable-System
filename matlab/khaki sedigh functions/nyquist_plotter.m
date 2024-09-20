function nyquist_plotter(g)
%NYQUIST_PLOTTER Summary of this function goes here
%   Detailed explanation goes here
step=10000;
wi=-12;
wf=12; 
e=eig(g);
m=length(e); 
freq1 =complex(0,-logspace(wf,wi,step));  
step=10000;
wi=-12;
wf=12; 
e=eig(g);
m=length(e); 
freq2 =complex(0,logspace(wi,wf,step)); 
freq=[freq1 freq2]; 
n=length(freq1); 
 
%%                                                                      MAPPING 
for k=1:m 
    D1=double(subs(e(k),freq1)); 
    D2=double(subs(e(k),freq2)); 
end 
switch class(D2) 
    case 'cell' 
    for k=1:m 
        D1=double(subs(e(k),freq1)); 
        D2=double(subs(e(k),freq2)); 
    end 
    for k=1:m 
        for w=1:n 
            nyq1(w,k)=D1{1,w}(1,k); 
            nyq2(w,k)=D2{1,w}(1,k); 
        end 
    end 
case 'double' 
    nyq1=zeros(n,m); 
    nyq2=zeros(n,m); 
    for k=1:m 
        nyq1(:,k)=subs(e(k),freq1); 
        nyq2(:,k)=subs(e(k),freq2); 
    end 
end 
 
%%                                                                     PLOTTING 
figure; hold on 
for k=1:m 
    switch k 
        case 1 
            co='b'; 
        case 2 
            co='r'; 
        case 3 
            co='g'; 
        case 4 
            co='m'; 
        case 5 
            co='y';   
        otherwise 
            co='k';      
    end 
    plot(real(nyq1(:,k)),imag(nyq1(:,k)),co);  
    plot(real(nyq2(:,k)),imag(nyq2(:,k)),co); 
    text(real(nyq2(fix(n/2.5),k)),imag(nyq2(fix(n/2.5),k)),[' \rightarrow \omega=',int2str(imag(freq2(fix(n/2.5))))],'FontSize',10,'fontweight','b') 
    text(real(nyq2(fix(n/4),k)),imag(nyq2(fix(n/4),k)),[' \rightarrow \omega=',int2str(imag(freq2(fix(n/4))))],'FontSize',10,'fontweight','b') 
end 
plot(-1,0,'m*') 
grid on

end

