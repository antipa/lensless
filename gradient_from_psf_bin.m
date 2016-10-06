function [f, g] = gradient_from_psf_bin(x,H,H_conj,W,crop_obj,Atb,obj_r,ds,upsamp)
%W: window function (2d rect)
%r: cropping boundaries
%x: variable 
%cr: crop rows
%cc: crop columns
%H: fft2(zeropad(psfd))
%H_conj: conj(H)
%Atb: ifft2(H_conj).*fft2(x);
%ds: relative pixel size after downsampling (1/downsample number)

%xup = upsample2d(x,upsamp,1,1);
xup = imresize(x,upsamp,'nearest');
b_u = ifftshift(ifft2(H.*fft2(xup)));   %crop(b_u) is same as Ax
AtAx = imresize(ifftshift(ifft2(H_conj.*fft2(W.*b_u))),1/upsamp,'box');
g = (AtAx-Atb)*ds;
b_un = W.*b_u;
f = norm(crop_obj(b_un)-obj_r,'fro')^2;