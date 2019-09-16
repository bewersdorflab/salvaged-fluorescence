% this example shows how to estimate the photon number in the salvaged fluorescence channel
%% load data in the conventional fluorescence channel
clear
clc
filestr='data_example1.mat';
if exist(filestr,'file')
    load(filestr);          % load data if exists    
else
    Conventional_fitting;   % otherwise, generate the data
end

%% read images of the salvaged fluorescence channel
filename='Salvaged_fluorescence.tif';
Im=single(tiffread(filename));
Im=imrotate(Im,90);
EMgain=200;                 % electron-multiplying(EM) gain of EMCCD
ADU=12;                     % ADU converstion factor, provided by the manufacture
Im=Im*ADU/EMgain;           % convert counts to photons

%% transform the positions to the salvaged fluorescence channel
load('two_channel_registration.mat'); % load the registration file between the two channels
x1=[];
x1(:,1)=VC(:,1);
x1(:,2)=VC(:,2);
x3=[];
XX=[ones(size(x1(:,1))),x1(:,1),x1(:,2),x1(:,1).*x1(:,2),x1(:,1).^2,x1(:,2).^2];
x3(:,1)=XX*betax;
x3(:,2)=XX*betay;
x3=x3*128/168;

%% estimate the background by median filter
totsz=size(Im,1);
bg=zeros(totsz,totsz);
for i=1:totsz
    for j=1:totsz
        I=single(Im(i,j,:));
        I=I(:);
        threshold=mean(I)+3*std(I);
        ix=I<threshold;
        I=I(ix);
        bg(i,j)=median(I);      % change to mean for faster calculation
    end
end
bg=repmat(bg,1,1,size(Im,3));

%% crop the corresponding sub regions
trans_x=x3(:,1)-1;
trans_y=x3(:,2)-1;
tc=VC(:,3)-1;
subsz=5;
im=single(Im)-bg;
[subims,t,l]=cMakeSubregions(round(trans_y),round(trans_x),tc,subsz,single(permute(im,[1 2 3])));
xcenter=trans_x(:)-l;
ycenter=trans_y(:)-t;

%% display example frame
f=100;
ix=VC(:,3)==f;
I=im(:,:,f);
figure;imshow(I,[0 max(I(:))*0.5]);hold on;plot(x3(ix,1),x3(ix,2),'bo');

%% generate weighted-Gaussian distributions
sigma=0.75;
sz=size(subims,1);
Model=single(zeros(sz,sz,size(subims,3)));
for ii=1:1
    yf=ycenter(:,ii);
    xf=xcenter(:,ii);
    ROI=finitegausspsf(sz,sigma,1,0,[xf,yf]);
    Model(:,:,:,ii)=single(ROI);
end

%% estimate the photons
L=length(VC);
P=zeros(L,3);
P(:,1)=PC;                      % photon number in the conventional fluorescence channel
for i=1:L
    I=subims(:,:,i);
    I(I<0)=0;
    model=Model(:,:,i);
    Inom=I.*(model);
    P(i,2)=sum(I(:));           % photon number in the salvaged fluorescence channel     
    P(i,3)=sum(Inom(:))*10;     % photon number in the salvaged fluorescence channel, weighted by Gaussian 
end
R=zeros(L,2);
R(:,1)=P(:,2)./P(:,1);
R(:,2)=P(:,3)./P(:,1);

%% display results on a logarithmic scale
figure;loglog(P(:,1),P(:,3),'r.','MarkerSize',0.5);
axis([100 10000 1 10000]);
xlabel('Conventional fluorescence (photon)');
ylabel('Salvaged fluorescence (AU)');
