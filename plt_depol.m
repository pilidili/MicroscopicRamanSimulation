%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright: Copyright (c) 2019
%Created on 2019-1-6 
%Author:MengDa (github:pilidili)
%Version 1.0 
%Title: MicroscopicRamanSimulation-plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear;clc;close all;

thickness=.5:.25:2;

rawdata_0=zeros(31,7,size(thickness,2));
rawdata_90=zeros(31,7,size(thickness,2));

for ii=1:size(thickness,2)
    name=['bulk_s_',strrep(num2str(thickness(ii)),'.','p'),'um'];
    rawdata_0(1:11,:,ii)=load([name,'_0d_out/w.txt']);
    indata=load([name,'_0d_in/w.txt']);
    rawdata_0(11:size(indata,1)+10,:,ii)=indata;
    rawdata_90(1:11,:,ii)=load([name,'_90d_out/w.txt']);
    indata=load([name,'_90d_in/w.txt']);
    rawdata_90(11:size(indata,1)+10,:,ii)=indata;
end
rawdata_0(rawdata_0==0)=nan;
rawdata_90(rawdata_90==0)=nan;

sc=['P_x';'P_y';'P_z'];
data_0=zeros(size(thickness,2),31,6);
data_90=data_0;
[X,Y]=meshgrid(0:.05:3,thickness);
drawdata=X;
drawdata(:,:,:)=nan;
for ii=1:3
    data_0(:,:,ii)=squeeze(rawdata_0(:,2*ii,:)+rawdata_0(:,2*ii+1,:))';
    data_90(:,:,ii)=squeeze(rawdata_90(:,2*ii,:)+rawdata_90(:,2*ii+1,:))';
    subplot(2,3,ii)
    for jj=1:size(thickness,2)
        drawdata(jj,1:floor((26+5*jj)/2))=data_0(jj,1:floor((26+5*jj)/2),ii);
        drawdata(jj,floor((26+5*jj)/2)+1:(26+5*jj))=fliplr(data_0(jj,1:ceil((26+5*jj)/2),ii));
    end
    imagesc(0:.05:3,thickness,drawdata);
    title(['z(x)激发 - ',sc(ii,:),'分量'])
    drawdata(:,:,:)=nan;
    view(0,90)
    colorbar;
    set(gca,'ytick',.5:.25:2);
    
    subplot(2,3,ii+3)
    for jj=1:size(thickness,2)
        drawdata(jj,1:floor((26+5*jj)/2))=data_90(jj,1:floor((26+5*jj)/2),ii);
        drawdata(jj,floor((26+5*jj)/2)+1:(26+5*jj))=fliplr(data_90(jj,1:ceil((26+5*jj)/2),ii));
    end
    imagesc(0:.05:3,thickness,drawdata);
    title(['z(y)激发 - ',sc(ii,:),'分量'])
    drawdata(:,:,:)=nan;
    view(0,90)
    colorbar;
    set(gca,'ytick',.5:.25:2);
end

% sc=['x(y)';'x(z)';'y(z)';'y(x)';'z(x)';'z(y)'];
% data_0=zeros(size(thickness,2),31,6);
% data_90=data_0;
% [X,Y]=meshgrid(0:.05:3,.5:.25:2);
% drawdata=X;
% drawdata(:,:,:)=nan;
% for ii=1:6
%     data_0(:,:,ii)=squeeze(rawdata_0(:,ii+1,:))';
%     data_90(:,:,ii)=squeeze(rawdata_90(:,ii+1,:))';
%     figure(ii)
%     subplot(2,1,1)
%     for jj=1:size(thickness,2)
%         drawdata(jj,1:floor((26+5*jj)/2))=data_0(jj,1:floor((26+5*jj)/2),ii);
%         drawdata(jj,floor((26+5*jj)/2)+1:(26+5*jj))=fliplr(data_0(jj,1:ceil((26+5*jj)/2),ii));
%     end
%     surf(X,Y,drawdata,'edgecolor','none');
%     title(['z(x)激发 - ',sc(ii,:),'散射'])
%     drawdata(:,:,:)=nan;
%     view(0,90)
%     colorbar;
%     subplot(2,1,2)
%     for jj=1:size(thickness,2)
%         drawdata(jj,1:floor((26+5*jj)/2))=data_90(jj,1:floor((26+5*jj)/2),ii);
%         drawdata(jj,floor((26+5*jj)/2)+1:(26+5*jj))=fliplr(data_90(jj,1:ceil((26+5*jj)/2),ii));
%     end
%     surf(X,Y,drawdata,'edgecolor','none');
%     title(['z(y)激发 - ',sc(ii,:),'散射'])
%     drawdata(:,:,:)=nan;
%     view(0,90)
% end
