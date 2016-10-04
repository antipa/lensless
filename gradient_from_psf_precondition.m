function [f, g] = gradient_from_psf_precondition(x,H,H_conj,W,cr,cc,Atb,obj_r,ds,p)
%W: window function (2d rect)
%r: cropping boundaries
%x: variable 
%cr: crop rows
%cc: crop columns
%H: fft2(zeropad(psfd))
%H_conj: conj(H)
%Atb: ifft2(H_conj).*fft2(x);
%ds: relative pixel size after downsampling (1/downsample number)
%p: precondition window function. This should normalize columns of A
b_u = ifftshift(ifft2(H.*fft2(x./p)));   %crop(b_u) is same as Ax
AtAx = ifftshift(ifft2(H_conj.*fft2(W.*b_u)))./p;
g = (AtAx-Atb)*ds.*(p>1e-4);
b_un = W.*b_u;
f = norm(b_un(cr,cc)-obj_r,'fro')^2;