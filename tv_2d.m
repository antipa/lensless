function [out, norm_out] = tv_2d(x,tau,niters,minval,maxval,pad)
denoised = min(max(tvdenoise(x,1/tau,niters),minval),maxval);
out = pad(denoised);
norm_out = tau*TVnorm(denoised);