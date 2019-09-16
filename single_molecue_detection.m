function [subims,tlz,locmaxc]=single_molecue_detection(im,thresh,subsz)

imsz=size(im,1);
[filteredim1]=varunif(squeeze(im),5);
[filteredim2]=varunif(squeeze(filteredim1),9);
im_unif=filteredim1-filteredim2;
L=size(im_unif,3);
for i=1:L
    A=im_unif(:,:,i);
    deta=median(median(abs(A-median(median(A)))));
    t=thresh*deta/0.67;
    A=(A>=t).*A;
    im_unif(:,:,i)=A;
end

se=ones(5,5);
dilatedI=imdilate(im_unif,se);
im_max=(im_unif==dilatedI)&(im_unif>0);
a=find(im_max);
z=floor(a/imsz/imsz);
pnum=mod(a,imsz*imsz);
y=mod(pnum,imsz)-1;
x=floor(pnum/imsz);
locmaxc=[x y z];

r=(subsz-1)/2;
rangemin=[r-1 r-1];
rangemax=[imsz-r+1 imsz-r+1];
mask=locmaxc(:,1)<rangemax(1)&locmaxc(:,1)>rangemin(1)...
    &locmaxc(:,2)<rangemax(2)&locmaxc(:,2)>rangemin(2);
locmaxc=locmaxc(mask,:);

if isempty(locmaxc)
    subims=[];
    tlz=[];
else
    [subims,t,l]=cMakeSubregions(locmaxc(:,2),locmaxc(:,1),locmaxc(:,3),subsz,single(permute(im,[1 2 3])));
    tlz=[t l locmaxc(:,3)];
end