%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright: Copyright (c) 2019
%Created on 2019-1-6 
%Author:MengDa (github:pilidili)
%Version 1.0 
%Title: MicroscopicRamanSimulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;close all;

format long;
n=2.5+.25*1i;%复折射率
NA=.9;%数值孔径
lambda=632.8e-9;%波长 m
theta2=asin(NA/real(n));%圆锥竖截面半顶角θ_2

ml=0.5e-6;%材料厚度m

x_edge=ml/2*1e6;
spx=linspace(-x_edge,-mod(x_edge,.05),floor(x_edge/.05)+1)*1e-6;%光源位置坐标x m
spy=0;%光源位置坐标y m
spz=2e-6;%光源位置坐标z m
clear x_edge;


W=zeros(size(spx,2),6);
t0=clock;
for ii=1:size(spx,2)
    disp(['progress: ',num2str(ii),'/',num2str(size(spx,2))])
    load(['Bulk_in_',num2str(ii)]);
    E1=E_1.E;
    x1=E_1.x;
    y1=E_1.y;
    z1=E_1.z;
    t=E_1.t;
    clear E_1;
    H1=H_1.H;
    clear H_1;
    
    E2=E_2.E;
    x2=E_2.x;
    y2=E_2.y;
    z2=E_2.z;
    clear E_2;
    H2=H_2.H;
    clear H_2;
    
    E3=E_3.E;
    x3=E_3.x;
    y3=E_3.y;
    z3=E_3.z;
    clear E_3;
    H3=H_3.H;
    clear H_3;
    
    E4=E_4.E;
    x4=E_4.x;
    y4=E_4.y;
    z4=E_4.z;
    clear E_4;
    H4=H_4.H;
    clear H_4;
    
    
    Ex1=reshape(E1(:,1,:),size(x1,1),size(y1,1),size(z1,1),size(t,1));
    Ey1=reshape(E1(:,2,:),size(x1,1),size(y1,1),size(z1,1),size(t,1));
    Ez1=reshape(E1(:,3,:),size(x1,1),size(y1,1),size(z1,1),size(t,1));
    clear E1;
    Ex2=reshape(E2(:,1,:),size(x2,1),size(y2,1),size(z2,1),size(t,1));
    Ey2=reshape(E2(:,2,:),size(x2,1),size(y2,1),size(z2,1),size(t,1));
    Ez2=reshape(E2(:,3,:),size(x2,1),size(y2,1),size(z2,1),size(t,1));
    clear E2;
    Ex3=reshape(E3(:,1,:),size(x3,1),size(y3,1),size(z3,1),size(t,1));
    Ey3=reshape(E3(:,2,:),size(x3,1),size(y3,1),size(z3,1),size(t,1));
    Ez3=reshape(E3(:,3,:),size(x3,1),size(y3,1),size(z3,1),size(t,1));
    clear E3;
    Ex4=reshape(E4(:,1,:),size(x4,1),size(y4,1),size(z4,1),size(t,1));
    Ey4=reshape(E4(:,2,:),size(x4,1),size(y4,1),size(z4,1),size(t,1));
    Ez4=reshape(E4(:,3,:),size(x4,1),size(y4,1),size(z4,1),size(t,1));
    clear E4;
    
    Hx1=reshape(H1(:,1,:),size(x1,1),size(y1,1),size(z1,1),size(t,1));
    Hy1=reshape(H1(:,2,:),size(x1,1),size(y1,1),size(z1,1),size(t,1));
    Hz1=reshape(H1(:,3,:),size(x1,1),size(y1,1),size(z1,1),size(t,1));
    clear H1;
    Hx2=reshape(H2(:,1,:),size(x2,1),size(y2,1),size(z2,1),size(t,1));
    Hy2=reshape(H2(:,2,:),size(x2,1),size(y2,1),size(z2,1),size(t,1));
    Hz2=reshape(H2(:,3,:),size(x2,1),size(y2,1),size(z2,1),size(t,1));
    clear H2;
    Hx3=reshape(H3(:,1,:),size(x3,1),size(y3,1),size(z3,1),size(t,1));
    Hy3=reshape(H3(:,2,:),size(x3,1),size(y3,1),size(z3,1),size(t,1));
    Hz3=reshape(H3(:,3,:),size(x3,1),size(y3,1),size(z3,1),size(t,1));
    clear H3;
    Hx4=reshape(H4(:,1,:),size(x4,1),size(y4,1),size(z4,1),size(t,1));
    Hy4=reshape(H4(:,2,:),size(x4,1),size(y4,1),size(z4,1),size(t,1));
    Hz4=reshape(H4(:,3,:),size(x4,1),size(y4,1),size(z4,1),size(t,1));
    clear H4;
    
    Pxy1=abs(Ey1.*Hz1);   Pyz1=abs(Ez1.*Hx1);   Pzx1=abs(Ex1.*Hy1);
    Pxz1=abs(-Ez1.*Hy1);  Pyx1=abs(-Ex1.*Hz1);  Pzy1=abs(-Ey1.*Hx1);
    clear Ex1 Ey1 Ez1 Hx1 Hy1 Hz1;
    
    Pxy2=abs(Ey2.*Hz2);   Pyz2=abs(Ez2.*Hx2);   Pzx2=abs(Ex2.*Hy2);
    Pxz2=abs(-Ez2.*Hy2);  Pyx2=abs(-Ex2.*Hz2);  Pzy2=abs(-Ey2.*Hx2);
    clear Ex2 Ey2 Ez2 Hx2 Hy2 Hz2;
    
    Pxy3=abs(Ey3.*Hz3);   Pyz3=abs(Ez3.*Hx3);   Pzx3=abs(Ex3.*Hy3);
    Pxz3=abs(-Ez3.*Hy3);  Pyx3=abs(-Ex3.*Hz3);  Pzy3=abs(-Ey3.*Hx3);
    clear Ex3 Ey3 Ez3 Hx3 Hy3 Hz3;
    
    Pxy4=abs(Ey4.*Hz4);   Pyz4=abs(Ez4.*Hx4);   Pzx4=abs(Ex4.*Hy4);
    Pxz4=abs(-Ez4.*Hy4);  Pyx4=abs(-Ex4.*Hz4);  Pzy4=abs(-Ey4.*Hx4);
    clear Ex4 Ey4 Ez4 Hx4 Hy4 Hz4;
    
    
    [Y1,X1,Z1]=meshgrid(y1,x1,z1);
    [Y2,X2,Z2]=meshgrid(y2,x2,z2);
    [Y3,X3,Z3]=meshgrid(y3,x3,z3);
    [Y4,X4,Z4]=meshgrid(y4,x4,z4);
    
    d2X1=(spx(ii)-X1).^2;
    d2Y1=(spy-Y1).^2;
    d2Z1=(spz-Z1).^2;
    clear Y1 Z1;
    
    d2X2=(spx(ii)-X2).^2;
    d2Y2=(spy-Y2).^2;
    d2Z2=(spz-Z2).^2;
    clear Y2 Z2;
    
    d2X3=(spx(ii)-X3).^2;
    d2Y3=(spy-Y3).^2;
    d2Z3=(spz-Z3).^2;
    clear Y3 Z3;
    
    d2X4=(spx(ii)-X4).^2;
    d2Y4=(spy-Y4).^2;
    d2Z4=(spz-Z4).^2;
    clear Y4 Z4;
    
    r1=sqrt(d2X1+d2Y1+d2Z1);
    r2=sqrt(d2X2+d2Y2+d2Z2);
    r3=sqrt(d2X3+d2Y3+d2Z3);
    r4=sqrt(d2X4+d2Y4+d2Z4);
    
    sinth1=sqrt(1-d2Z1./(d2X1+d2Y1+d2Z1));
    sinth2=sqrt(1-d2Z2./(d2X2+d2Y2+d2Z2));
    sinth3=sqrt(1-d2Z3./(d2X3+d2Y3+d2Z3));
    sinth4=sqrt(1-d2Z4./(d2X4+d2Y4+d2Z4));
    
    screen_matrix1=d2Z1;
    screen_matrix1(screen_matrix1<(d2X1+d2Y1)/(tan(theta2)^2))=0;
    screen_matrix1(screen_matrix1~=0)=1;
    screen_matrix1(X1<-ml/2)=0;
    clear X1 d2X1 d2Y1 d2Z1;
    
    screen_matrix2=d2Z2;
    screen_matrix2(screen_matrix2<(d2X2+d2Y2)/(tan(theta2)^2))=0;
    screen_matrix2(screen_matrix2~=0)=1;
    screen_matrix2(X2<-ml/2)=0;
    clear X2 d2X2 d2Y2 d2Z2;
    
    screen_matrix3=d2Z3;
    screen_matrix3(screen_matrix3<(d2X3+d2Y3)/(tan(theta2)^2))=0;
    screen_matrix3(screen_matrix3~=0)=1;
    screen_matrix3(X3<-ml/2)=0;
    clear X3 d2X3 d2Y3 d2Z3;
    
    screen_matrix4=d2Z4;
    screen_matrix4(screen_matrix4<(d2X4+d2Y4)/(tan(theta2)^2))=0;
    screen_matrix4(screen_matrix4~=0)=1;
    screen_matrix4(X4<-ml/2)=0;
    clear X4 d2X4 d2Y4 d2Z4;
    
    W(ii,1)=sum(sum(sum(sum(screen_matrix1.*Pxy1.*sinth1(:,:,:,ones(1,size(t,1))).*exp(-2*imag(n)/lambda*r1(:,:,:,ones(1,size(t,1))))))));
    W(ii,1)=W(ii,1)+sum(sum(sum(sum(screen_matrix2.*Pxy2.*sinth2(:,:,:,ones(1,size(t,1))).*exp(-2*imag(n)/lambda*r2(:,:,:,ones(1,size(t,1))))))));
    W(ii,1)=W(ii,1)+sum(sum(sum(sum(screen_matrix3.*Pxy3.*sinth3(:,:,:,ones(1,size(t,1))).*exp(-2*imag(n)/lambda*r3(:,:,:,ones(1,size(t,1))))))));
    W(ii,1)=W(ii,1)+sum(sum(sum(sum(screen_matrix4.*Pxy4.*sinth4(:,:,:,ones(1,size(t,1))).*exp(-2*imag(n)/lambda*r4(:,:,:,ones(1,size(t,1))))))));
    W(ii,2)=sum(sum(sum(sum(screen_matrix1.*Pxz1.*sinth1(:,:,:,ones(1,size(t,1))).*exp(-2*imag(n)/lambda*r1(:,:,:,ones(1,size(t,1))))))));
    W(ii,2)=W(ii,2)+sum(sum(sum(sum(screen_matrix2.*Pxz2.*sinth2(:,:,:,ones(1,size(t,1))).*exp(-2*imag(n)/lambda*r2(:,:,:,ones(1,size(t,1))))))));
    W(ii,2)=W(ii,2)+sum(sum(sum(sum(screen_matrix3.*Pxz3.*sinth3(:,:,:,ones(1,size(t,1))).*exp(-2*imag(n)/lambda*r3(:,:,:,ones(1,size(t,1))))))));
    W(ii,2)=W(ii,2)+sum(sum(sum(sum(screen_matrix4.*Pxz4.*sinth4(:,:,:,ones(1,size(t,1))).*exp(-2*imag(n)/lambda*r4(:,:,:,ones(1,size(t,1))))))));
    clear Pxy1 Pxy2 Pxy3 Pxy4 Pxz1 Pxz2 Pxz3 Pxz4;
    
    W(ii,3)=sum(sum(sum(sum(screen_matrix1.*Pyz1.*sinth1(:,:,:,ones(1,size(t,1))).*exp(-2*imag(n)/lambda*r1(:,:,:,ones(1,size(t,1))))))));
    W(ii,3)=W(ii,3)+sum(sum(sum(sum(screen_matrix2.*Pyz2.*sinth2(:,:,:,ones(1,size(t,1))).*exp(-2*imag(n)/lambda*r2(:,:,:,ones(1,size(t,1))))))));
    W(ii,3)=W(ii,3)+sum(sum(sum(sum(screen_matrix3.*Pyz3.*sinth3(:,:,:,ones(1,size(t,1))).*exp(-2*imag(n)/lambda*r3(:,:,:,ones(1,size(t,1))))))));
    W(ii,3)=W(ii,3)+sum(sum(sum(sum(screen_matrix4.*Pyz4.*sinth4(:,:,:,ones(1,size(t,1))).*exp(-2*imag(n)/lambda*r4(:,:,:,ones(1,size(t,1))))))));
    W(ii,4)=sum(sum(sum(sum(screen_matrix1.*Pyx1.*sinth1(:,:,:,ones(1,size(t,1))).*exp(-2*imag(n)/lambda*r1(:,:,:,ones(1,size(t,1))))))));
    W(ii,4)=W(ii,4)+sum(sum(sum(sum(screen_matrix2.*Pyx2.*sinth2(:,:,:,ones(1,size(t,1))).*exp(-2*imag(n)/lambda*r2(:,:,:,ones(1,size(t,1))))))));
    W(ii,4)=W(ii,4)+sum(sum(sum(sum(screen_matrix3.*Pyx3.*sinth3(:,:,:,ones(1,size(t,1))).*exp(-2*imag(n)/lambda*r3(:,:,:,ones(1,size(t,1))))))));
    W(ii,4)=W(ii,4)+sum(sum(sum(sum(screen_matrix4.*Pyx4.*sinth4(:,:,:,ones(1,size(t,1))).*exp(-2*imag(n)/lambda*r4(:,:,:,ones(1,size(t,1))))))));
    clear Pyz1 Pyz2 Pyz3 Pyz4 Pyx1 Pyx2 Pyx3 Pyx4; 
    
    W(ii,5)=sum(sum(sum(sum(screen_matrix1.*Pzx1.*sinth1(:,:,:,ones(1,size(t,1))).*exp(-2*imag(n)/lambda*r1(:,:,:,ones(1,size(t,1))))))));
    W(ii,5)=W(ii,5)+sum(sum(sum(sum(screen_matrix2.*Pzx2.*sinth2(:,:,:,ones(1,size(t,1))).*exp(-2*imag(n)/lambda*r2(:,:,:,ones(1,size(t,1))))))));
    W(ii,5)=W(ii,5)+sum(sum(sum(sum(screen_matrix3.*Pzx3.*sinth3(:,:,:,ones(1,size(t,1))).*exp(-2*imag(n)/lambda*r3(:,:,:,ones(1,size(t,1))))))));
    W(ii,5)=W(ii,5)+sum(sum(sum(sum(screen_matrix4.*Pzx4.*sinth4(:,:,:,ones(1,size(t,1))).*exp(-2*imag(n)/lambda*r4(:,:,:,ones(1,size(t,1))))))));
    W(ii,6)=sum(sum(sum(sum(screen_matrix1.*Pzy1.*sinth1(:,:,:,ones(1,size(t,1))).*exp(-2*imag(n)/lambda*r1(:,:,:,ones(1,size(t,1))))))));
    W(ii,6)=W(ii,6)+sum(sum(sum(sum(screen_matrix2.*Pzy2.*sinth2(:,:,:,ones(1,size(t,1))).*exp(-2*imag(n)/lambda*r2(:,:,:,ones(1,size(t,1))))))));
    W(ii,6)=W(ii,6)+sum(sum(sum(sum(screen_matrix3.*Pzy3.*sinth3(:,:,:,ones(1,size(t,1))).*exp(-2*imag(n)/lambda*r3(:,:,:,ones(1,size(t,1))))))));
    W(ii,6)=W(ii,6)+sum(sum(sum(sum(screen_matrix4.*Pzy4.*sinth4(:,:,:,ones(1,size(t,1))).*exp(-2*imag(n)/lambda*r4(:,:,:,ones(1,size(t,1))))))));
    clear Pzx1 Pzx2 Pzx3 Pzx4 Pzy1 Pzy2 Pzy3 Pzy4;
    clear sinth1 sinth2 sinth3 sinth4 r1 r2 r3 r4 screen_matrix1 screen_matrix2 screen_matrix3 screen_matrix4
    TimeCost=etime(clock,t0);
    clc;
    disp(['TimeCost:',num2str(floor(TimeCost/60)),' min ',num2str(mod(TimeCost,60)),' s']);
    TimeLeft=TimeCost/ii*(size(spx,2)-ii);
    disp(['TimeLeft:',num2str(floor(TimeLeft/60)),' min ',num2str(mod(TimeLeft,60)),' s']);
end
W=W*2;
for ii=[1 3 5]
    figure(ii)
    subplot(2,1,1)
    plot(spx,W(:,ii));
    subplot(2,1,2)
    plot(spx,W(:,ii+1));
end

output_var=[spx',W];
save('W.txt','output_var','-ascii')
TimeCost=etime(clock,t0);
clc;
disp(['TimeCost:',num2str(floor(TimeCost/60)),' min ',num2str(mod(TimeCost,60)),' s']);
disp('--Finished!--')
load chirp
sound(y,Fs);