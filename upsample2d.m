function out = upsample2d(in,n,phase_r,phase_c)
%zeros out n-1 pixels starting at pixel [phase_r,phase_c]
r = upsample(ones(size(in,1),1),n,phase_r);
c = upsample(ones(1,size(in,2)),n,phase_c);
m = r*c;
out = m.*imresize(in,n,'box');
