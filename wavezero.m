function [out, norm_out] = wavezero(x,wavelev,wavetype,zvec,minval,maxval,pad)
%Do wavelet decomposition, then zero out indices in zvec, then recompose.
%This is only for debiasing and the zeros must be found in advance.

[C,S] = wavedec2(x,wavelev,wavetype);
C(zvec)=0;
out = pad(min(max(waverec2(C,S,wavetype),minval),maxval));
norm_out = 0;
%norm(C,1);
