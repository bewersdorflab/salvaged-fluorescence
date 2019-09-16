% this script shows how to perform single-molecule localization in the
% conventional fluorescence channel and estimate the positions and photon
% numbers
%% read images of the conventiola fluorescence channel
clear
clc
addpath('fit3Dcspline');
filename='Conventional_fluorescence.tif';
Im=single(tiffread(filename));

%% single-molecule detection
sz=size(Im,1);      % image size
subsz=7;            % sub region size
det_threshold=4;    % threshold for single-molecule detection
[sub_regions,tlz,locmaxc1]=single_molecue_detection(Im,det_threshold,subsz);

%% MLE GPU fitting
[P,CRLB,LL]=mleFit_LM(sub_regions,2,50,1.1,0,0,0); 
xf=P(:,2);
yf=P(:,1);
llr=-2*LL;
CRLB(:,1:2)=sqrt(CRLB(:,1:2));
xest=xf+tlz(:,2);
yest=yf+tlz(:,1);
VC=zeros(length(xf),3);
VC(:,1)=xest+1;
VC(:,2)=yest+1;
VC(:,3)=tlz(:,3)+1;

%% reject molecules with bad fitting results
r=(subsz-1)/2;
pixel=128;          % pixel size nm
crlb=sqrt((CRLB(:,1).^2+CRLB(:,2).^2)/2)*pixel;
mask=abs(xf-r)<1 & abs(yf-r)<1 & crlb<40 & llr<100;
VC=VC(mask,:);      % x position, y position, frame number
crlb=crlb(mask);    % cramer-rao lower bound for xy positions, nm
PC=P(mask,3);       % photon number in the conventional fluorescence channel
llr=llr(mask);      % log-likelihood ratio

%% show example localization result
f=100;
V0=VC(VC(:,3)==f,:);
figure;imshow(Im(:,:,f),[30 300]);hold on;plot(V0(:,1),V0(:,2),'bo');pause(1)

%% save result
filestr='data_example1.mat';
save(filestr,'VC','PC','crlb','llr');
