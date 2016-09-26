function [f, g] = gradient_from_psf_bin(x,H,H_conj,W_obj,W,crop_obj,Atb,obj_r,ds,filt,downsamp)
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
AtAx = ifftshift(ifft2(H_conj.*fft2(W.*imfilter(b_u,filt))));
g = (AtAx-Atb)*ds;
b_un = W_obj.*imresize(b_u,downsamp,'box');
f = norm(crop_obj(b_un)-obj_r,'fro')^2;