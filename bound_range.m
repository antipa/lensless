function [out, norm_out] = bound_range(x,minval,maxval,pad)
out = pad(min(max(x,minval),maxval));
norm_out = 0;