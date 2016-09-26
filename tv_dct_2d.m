function [out, norm_out] = tv_dct_2d(x,tau_tv,tau_dct,niters,minval,maxval,pad)
%[out, norm_out] = tv_dct_2d(x,tau_tv,tau_dct,niters,minval,maxval,pad)
dct_vals = dct2(x);
d = soft(dct_vals,tau_dct);
x_dct = idct2(d);
out = pad(min(max(tvdenoise(x_dct,2/tau_tv,niters),minval),maxval));
norm_out = 1/2 * tau_dct*norm(d(:),1) + 1/2 * tau_tv * TVnorm(out);