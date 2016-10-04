function [f, g] = gradient_from_psf(x,H,H_conj,W,cr,cc,Atb,obj_r,ds)
%W: window function (2d rect)
%r: cropping boundaries
%x: variable 
%cr: crop rows
%cc: crop columns
%H: fft2(zeropad(psfd))
%H_conj: conj(H)
%Atb: ifft2(H_conj).*fft2(x);
%ds: relative pixel size after downsampling (1/downsample number)

b_u = ifftshift(ifft2(H.*fft2(x)));   %crop(b_u) is same as Ax
AtAx = ifftshift(ifft2(H_conj.*fft2(W.*b_u)));
g = (AtAx-Atb)*ds;
b_un = W.*b_u;
f = norm(b_un(cr,cc)-obj_r,'fro')^2;         