function [out, norm_out,varargout] = adaptive_soft_wvlt_2d(x,wavelev,wavetype,pct_nz,minval,maxval,pad)
%Inputs:
%x: variable
%wavelev: wavelet level for wavedec
%wavetype: wavelet type
%tau: soft threshold amount
%minval: minimum allowed value in final output
%maxval: same but max
%W: windowing function. Optional, but needed for diffuser 2d problem. Set
%to ones if not using.
[C,S] = wavedec2(x,wavelev,wavetype);
Csort = sort(abs(C),'descend');
cutoff = floor(pct_nz*length(Csort));
tau = Csort(cutoff);

C_soft = soft(C,tau);
norm_out = tau*norm(C,1);
out = pad(min(max(waverec2(C_soft,S,wavetype),minval),maxval));
if nargout >0
    varargout{1} = tau;
end
