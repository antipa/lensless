function [out, norm_out] = soft_dct2(x,tau,minval,maxval,W,pad)
dct_vals = dct2(x);
d = soft(dct_vals,tau);
out = pad(min(max(idct2(d),minval),maxval));
norm_out = tau*norm(d(:),1);

