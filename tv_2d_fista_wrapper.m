function [out, norm_out] = tv_2d_fista_wrapper(x,tau,niters,minval,maxval,pad,use_gpu)
%[out, norm_out] = tv_2d_fista_wrapper(x,tau,niters,minval,maxval,pad,use_gpu)
a.MAXITER = niters;
a.gpu = use_gpu;
denoised = TV2DFista(x,tau,minval,maxval,a);
out = pad(denoised);
norm_out = tau*TVnorm(denoised);