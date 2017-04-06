function [out, norm_out] = tv_2d_fista_wrapper(x,tau,niters,minval,maxval,pad)
a.MAXITER = niters;
denoised = TV2DFista(x,tau,minval,maxval,a);
out = pad(denoised);
norm_out = tau*TVnorm(denoised);