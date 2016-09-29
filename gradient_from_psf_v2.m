function [f, g] = gradient_from_psf_v2(x,H,H_conj,W,inds,Atb,obj_r,ds)
%W: window function (2d rect)
%r: cropping boundaries
%x: variable 
%inds: pull out nonzero entries in forward model
%H: fft2(zeropad(psfd))
%H_conj: conj(H)
%Atb: ifft2(H_conj).*fft2(x);
%ds: relative pixel size after downsampling (1/downsample number)
%obj_r: input object. Enter as a vector!
b_u = ifft2(H.*fft2(x));  %crop(b_u) is same as Ax
AtAx = ifft2(H_conj.*fft2(W.*b_u));
g = (AtAx-Atb)*ds;
f = norm(b_u(inds)-obj_r)^2;