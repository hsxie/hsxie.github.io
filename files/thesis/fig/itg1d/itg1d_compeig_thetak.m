% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2014-12-18 21:07
% companion matrix eig to solve 1D ITG eigen equation in ballooning space
% Cite: [1] H. S. Xie & Y. Xiao, Unconventional ballooning structures for
%           toroidal drift waves, Physics of Plasmas, 22, 090703 (2015).
%       [2] H. S. Xie, "Numerical Simulations of Micro-turbulence in
%           Tokamak Edge", PhD thesis, Chap 6, Zhejiang University, 2015.

% Default to calculate Fig.5 in Ref.[1]
% CPU time ~ O(Nt^2.6), Nt=1024, ~ 1000s

close all; clear; clc;

pltc=[0.0  1.0  0.0
    1.0  0.0  0.0
    0.2  0.2  1.0
    0.8  0.8  0.0
    1.0  0.6  0.0
    0.9  0.0  0.9
    0.0  0.8  0.8
    0.0  0.0  0.0
    0.6  0.0  0.0
    0.4  0.7  0.4 
    0.0  0.0  0.5 
    0.6  0.0  0.6 
    0.0  0.5  1.0];

sp=0; % sparse matrix
par=1;
if(par==1)
    % Dickinson14's case
%     parstr='n=50; s=2.0; ktrhoi=0.33; R=10; q=1.8; etai=5.0; epsn=0.03;';
%     parstr='n=50; s=2.0; ktrhoi=0.33; R=10; q=1.8; etai=5.0; epsn=0.90;';
    parstr='n=50; s=0.8; ktrhoi=0.4; R=10; q=1.0; etai=3.0; epsn=0.5;';
    eval(parstr);
    
%     sm=-0.5+1.5i;
    sm=-2.0+2.5i;
else
    % LuZX12's case
%     parstr='s=1.0; ktrhoi=0.3; q=2.0; etai=5.4; epsn=0.6; wn=1.0;';
    parstr='s=2.0; ktrhoi=0.1; q=2.0; etai=2.0; epsn=0.1; wn=1.0;';
    eval(parstr);
    sgm=1.0; 
    
    sm=0.5+3.0i;
end

t0=0.0; % t -> theta, t0 -> theta0

Nt=1024; tmax=20; 
% Nt=1024*3; tmax=40; 

tmin=-tmax;
dt=(tmax-tmin)/(Nt);
dt2=dt*dt;
t=(tmin+dt):dt:(tmax);
% t00=-pi:pi/4:pi;
% t00=1.0*pi:-0.2*pi:0.0*pi;
% t00=0.1*pi:-0.1*pi:0.0*pi;
t00=0.0*pi;
nt0=length(t00);
ww=[]; VV=[]; Xjt=[];
ip=1; ipl=[]; covv=[]; 
epsn00=[0.5,0.2];

h=figure('unit','normalized','Position',[0.01 0.17 0.5 0.65],...
    'DefaultAxesFontSize',15);
for jp=1:2
    
    epsn=epsn00(jp);
    sgm=epsn/(q*ktrhoi); 
    
    fa=ktrhoi^2*(1.0+s^2*(t-t0).^2);
    fb=2.0*epsn*(cos(t)+s*(t-t0).*sin(t));
%     fb=0.3*2.0*epsn*(cos(t)+s*(t-t0).*sin(t));
    
    II=1:Nt; JJ=II; SS=1+0.*II;
    I0=zeros(3*Nt-2,1); J0=I0; S0=I0; I0(1)=1; J0(1)=1; 
    I1=zeros(3*Nt-2,1); J1=I1; S1=I1; I1(1)=1; J1(1)=1; 
    I2=zeros(Nt,1); J2=I2; S2=I2;
    I3=zeros(Nt,1); J3=I3; S3=I3;
    j=1;
    S0jstr='S0(j)=-2.0*etai*sgm^2/(dt2);';
    S1jstr='S1(j)=-2.0*sgm^2/(dt2)+etai*fb(j);';
    
    if(par==1)
        S2jstr='S2(j)=fb(j)+etai*fa(j)-1.0;';
        S3jstr='S3(j)=fa(j)+1.0;';
    else
        S2jstr='S2(j)=fb(j)+etai*fa(j)-wn;';
        S3jstr='S3(j)=fa(j)+1.0;';
    end
    
    eval(S0jstr);
    eval(S1jstr);
    
    for j=2:Nt
       I0(j)=j; J0(j)=j; 
       eval(S0jstr);
       I0(Nt+j-1)=j; J0(Nt+j-1)=j-1; S0(Nt+j-1)=etai*sgm^2/(dt2);
       I0(2*Nt+j-2)=j-1; J0(2*Nt+j-2)=j; S0(2*Nt+j-2)=etai*sgm^2/(dt2);
       
       I1(j)=j; J1(j)=j; 
       eval(S1jstr);
       I1(Nt+j-1)=j; J1(Nt+j-1)=j-1; S1(Nt+j-1)=sgm^2/(dt2);
       I1(2*Nt+j-2)=j-1; J1(2*Nt+j-2)=j; S1(2*Nt+j-2)=sgm^2/(dt2);       
    end
    
    for j=1:Nt
       I2(j)=j; J2(j)=j;
       eval(S2jstr);
       
       I3(j)=j; J3(j)=j;
       eval(S3jstr);
    end

    MS0=sparse(I0,J0,S0,Nt,Nt);
    MS1=sparse(I1,J1,S1,Nt,Nt);
    MS2=sparse(I2,J2,S2,Nt,Nt);
    MS3=sparse(I3,J3,S3,Nt,Nt);
    %
    Nt2=3*Nt;
    MSA=sparse(II,JJ+Nt,SS,Nt2,Nt2)+sparse(II+Nt,JJ+2*Nt,SS,Nt2,Nt2)+...
        sparse(I2+2*Nt,J2+2*Nt,-S2,Nt2,Nt2)+sparse(I1+2*Nt,J1+Nt,-S1,Nt2,Nt2)+...
        sparse(I0+2*Nt,J0,-S0,Nt2,Nt2);
    %
    MSB=sparse(II,JJ,SS,Nt2,Nt2)+sparse(II+Nt,JJ+Nt,SS,Nt2,Nt2)+...
        sparse(I3+2*Nt,J3+2*Nt,S3,Nt2,Nt2);

    if(sp==0)
        MA=full(MSA); MB=full(MSB);
        e=eig(MA,MB);
    else
%         e=eigs(MSA,MSB,5,0.5+3.0i);
        e=eigs(MSA,MSB,100,sm);
%         e=eigs(MSA,MSB,5,'li');
    end
%     e=polyeig(MS0,MS1,MS2,MS3);
%     [X,e]=polyeig(MS0,MS1,MS2,MS3);
    ind=find(imag(e)==max(imag(e)));
%     ind=find(abs(e)<=1.1*min(abs(e)));

    w=e(ind(1));
    ww=[ww,w];
    
    [X,D]=eigs(MSA,MSB,5,w);
    
    Xjt=[Xjt,X(1:Nt,1)];    
    
    ya=-2.0; yb=2.0;
    xa=-5.0; xb=1.0;
    subplot(2,2,2*jp-1);
    plot(real(e),imag(e),'+',real(w),imag(w),'rx','Linewidth',2);
    xlabel('\omega_r'); ylabel('\gamma');
%     title('(c) all solutions for \theta_k=0.0');
    xlim([xa,xb]); ylim([ya,yb]);
    text(0.98*xa+0.02*xb,0.1*ya+0.9*yb,['(',char(2*jp+95),') \epsilon_n=',num2str(epsn)],'FontSize',14);
    %
    ya=-1.0; yb=1.4;
    xa=-3.0; xb=3.0;
    subplot(2,2,2*jp);
    plot(t/pi,real(X(1:Nt,1))/max(abs(X(1:Nt,1))),t/pi,imag(X(1:Nt,1))/max(abs(X(1:Nt,1))),'--','Linewidth',2);
    xlabel('\vartheta/\pi'); %ylabel('\delta\phi');
    ylabel('$\delta\hat{\phi}$','interpreter','latex');
    xlim([xa,xb]); ylim([ya,yb]);
    text(0.98*xa+0.02*xb,0.1*ya+0.9*yb,['(',char(2*jp+96),') \omega=',num2str(num2str(D(1,1)))],'FontSize',14);
    
%     title(['(d) \omega=',num2str(D(1,1))]); axis tight;
    hlgd=legend('Re','Im',4); legend('boxoff');
    hLegendLines = findobj(hlgd, 'type', 'line', '-and', '-regexp','Tag','[^'']');
    set(hLegendLines, 'XData', [.45, 0.6]);
    
end

%
% save;
print(gcf,'-dpng',['itg_par=',num2str(par),',etai=',num2str(etai),...
    ',s=',num2str(s),',q=',num2str(q),',sp=',num2str(sp),...
    ',tmax=',num2str(tmax),',nt=',num2str(Nt),'.png']);


