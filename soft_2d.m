function [out,norm_out] = soft_2d(x,tau,minval,maxval)
out = min(max(soft(x,tau),minval),maxval);
norm_out = tau*norm(out(:),1);