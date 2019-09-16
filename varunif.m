function unifim=varunif(im,sz)
if sz==5
    kernel=[1/16,1/4,3/8,1/4,1/16];
    km=kernel'*kernel;
elseif sz==9
    kernel=[1/16,0,1/4,0,3/8,0,1/4,0,1/16];
    km=kernel'*kernel;
else
    km=ones(sz,sz);
end
unifim=convn(im,km,'same');