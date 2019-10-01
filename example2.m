% this example shows how to perform color assignment based on the photon
% numbers in the conventional and salvaged fluorescence channel
%% load example data set, the same data shown in Supplementary Fig. 13D
clear
clc
load('data_example2.mat');

%% generare a mask for each channel
N=length(Photon);
P=[Photon(:,1),Photon(:,3)];
P0=log10(P);
R=P(:,2)./P(:,1);

threshold1=0.4;
threshold2=0.6;
level1=0.03;
level2=0.03;

ix=R<threshold1;
P1=P0(ix,:);
X1=(P1(:,1)-2.2)*200;
Y1=(P1(:,2))*60;
dmap=cHistRecon(300,300,single(Y1),single(X1),0);
dmap=double(imgaussfilt(dmap,2));
dmap=dmap/max(dmap(:));
BW=imbinarize(dmap,level1);
se=strel('disk',5);
BW=imclose(BW,se);
B=bwboundaries(BW);
B1=B{1};

ix=R>threshold2;
P1=P0(ix,:);
X1=(P1(:,1)-2.2)*200;
Y1=(P1(:,2))*60;
dmap=cHistRecon(300,300,single(Y1),single(X1),0);
dmap=double(imgaussfilt(dmap,2));
dmap=dmap/max(dmap(:));
BW=imbinarize(dmap,level2);
se=strel('disk',5);
BW=imclose(BW,se);
B=bwboundaries(BW);
B2=B{1};

%% display masks on histogram image
X=(P0(:,1)-2.2)*200;
Y=P0(:,2)*60;
dmap=cHistRecon(300,300,single(Y),single(X),0);
dmap=double(imgaussfilt(dmap,2));
dmap=dmap/max(dmap(:));
figure;imshow(dmap,[]); hold on; axis xy
plot(B2(:,2),B2(:,1),'g-','MarkerSize',2);
plot(B1(:,2),B1(:,1),'m-','MarkerSize',2);
xlabel('Conventional fluorescence');
ylabel('Salvaged fluorescence');
h=legend('AF647','CF660C');
set(h,'linewidth',3,'FontSize',16);
pause(0.1)

%% assign molecules to each channel
id1=inpolygon(X,Y,B1(:,2),B1(:,1));
id2=inpolygon(X,Y,B2(:,2),B2(:,1));

%% display channel assignment result
xout=VC(:,1);
yout=VC(:,2);
pixel=10; % pixel size for recontruction, nm
sz=168*128/pixel;

coords=[xout(id1),yout(id1)]/pixel;
im1=cHistRecon(sz,sz,coords(:,1),coords(:,2),0);
im1=single(im1);
im1=im1/max(im1(:))*10;

coords=[xout(id2),yout(id2)]/pixel;
im2=cHistRecon(sz,sz,coords(:,1),coords(:,2),0);
im2=single(im2);
im2=im2/max(im2(:))*10;

rgbImage=cat(3,im1,im2,im1);
figure;imshow(rgbImage);
