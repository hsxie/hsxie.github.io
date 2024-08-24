% Hua-sheng XIE, IFTS-ZJU, huashengxie@gmail.com, 2013-01-06 15:57
% Tokamak Ballooning Mode Visualization, draft
% need extend: multi-modes coupling

close all;clear;clc;
R0=3.0;a=R0/3;kappa=1.5;delta=0.5;b=0.0;r0=0.8*a;dr=0.1*r0;

figure('Unit','normalized','Position',[0.01 0.3 0.6 0.6]);
set(gcf,'DefaultAxesFontSize',15);

h=axes('position',[0 0 1 0.92]);axis(h);
title('Tokamak Ballooning Mode Visualization (Artificial)'); hold on;
axis tight; axis off; hidden off;
for ipl=1:3
    m=5*ipl; n=4*ipl;

    % X=R.*cos(phi);
    % Y=R.*sin(phi);
    % Z=r0.*sin(theta);

    % without flux surface shift
    [theta1,phi1]=meshgrid(0:2*pi/600:2*pi,-0.5*pi:2*pi/500:1.5*pi);
    envelope1=exp(-0.2.*(abs(theta1-pi)-pi).^2);
    X1=(R0-b+(r0+b.*cos(theta1)).*cos(theta1+delta.*sin(theta1))).*cos(phi1);
    Y1=(R0-b+(r0+b.*cos(theta1)).*cos(theta1+delta.*sin(theta1))).*sin(phi1);
    Z1=(kappa*r0).*sin(theta1);
    F1=(dr.*envelope1).*cos(m.*theta1-n.*phi1);
    wid=0.45;
    h1=axes('position',[1/6-wid/2+(ipl-1)/3 0.3 wid 0.75]);axis(h1);
%     subplot(2,3,ipl);
    h3=surf(X1,Y1,Z1,F1,'EdgeColor','none');axis equal; hold on;
%     colorbar('location','South'); 
    str=['(m=',num2str(m),', n=',num2str(n),')'];
    title(str);axis equal; axis tight; axis off; hidden off;

    % with surf shift
    [theta2,phi2]=meshgrid(0:2*pi/600:2*pi,-0.5*pi:2*pi/500:1.2*pi);
    envelope2=exp(-0.2.*(abs(theta2-pi)-pi).^2);
    r=r0+(dr.*envelope2).*cos(m.*theta2-n.*phi2);
    X2=(R0-b+(r+b.*cos(theta2)).*cos(theta2+delta.*sin(theta2))).*cos(phi2);
    Y2=(R0-b+(r+b.*cos(theta2)).*cos(theta2+delta.*sin(theta2))).*sin(phi2);
    Z2=(kappa*r).*sin(theta2);
    F2=(dr.*envelope2).*cos(m.*theta2-n.*phi2);
%     subplot(2,3,ipl+3);
    wid=0.45;
    h2=axes('position',[1/6-wid/2+(ipl-1)/3 -0.2 wid 0.8]);axis(h2);
    h4=surf(X2,Y2,Z2,F2,'EdgeColor','none'); axis equal;hold on;
    title(str);axis equal; axis tight; axis off; hidden off;

    % add a poloidal cross section
    phi3=-0.46*pi;
    [r3,theta3]=meshgrid(0:0.01:a,0:2*pi/600:2*pi);
    envelope3=exp(-0.2.*(abs(theta3-pi)-pi).^2);
    r30=r3+(dr.*envelope3).*cos(m.*theta3-n*phi3);
    envr=exp(-30.*(r30-r0).^2);
    X3=(R0-b+(r30+b.*cos(theta3)).*cos(theta3+delta.*sin(theta3))).*cos(phi3);
    Y3=(R0-b+(r30+b.*cos(theta3)).*cos(theta3+delta.*sin(theta3))).*sin(phi3);
    Z3=(kappa.*r30).*sin(theta3);
    F3=envr.*(dr.*envelope3).*cos(m.*theta3-n*phi3);
    surf(X3,Y3,Z3,F3,'EdgeColor','none');
end
pstr=['ballooning_mode,kappa=',num2str(kappa),'delta=',num2str(kappa),'b=',num2str(b)];
print('-dpng',strcat(pstr,'.png'));
set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto');
print('-dpdf',strcat(pstr,'.pdf'));
print('-deps',strcat(pstr,'.eps'));
