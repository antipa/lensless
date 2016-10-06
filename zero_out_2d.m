function out = zero_out_2d(in,comb)
%Only leave every pth entry nonzero in the input, in
out = zeros(size(in));
u = ceil(p/2):p:size(in,1);
v = ceil(p/2):p:size(in,2);
[U,V] = meshgrid(u,v);
i = sub2ind(size(in),U,V);
out(i) = in(i);
