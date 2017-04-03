function [out, norm_out] = wavelet_detail_denoise(x,wavelev,wavetype,denoise_level,tau,minval,maxval,pad)
%Inputs:
%x: variable
%wavelev: wavelet level for wavedec
%wavetype: wavelet type
%denoise_level:  array containing list of detail coefficient levels to apply soft thresholding
%tau: soft threshold amount
%minval: minimum allowed value in final output
%maxval: same but max
%W: windowing function. Optional, but needed for diffuser 2d problem. Set
%to ones if not using.
[C,S] = wavedec2(x,wavelev,wavetype);
C_soft = C;
norm_out = 0;
for n = 1:numel(denoise_level)
    C_soft = soft_wave_detail('c',C_soft,S,denoise_level(n),tau);
    norm_out = norm_out + tau*norm(detcoef2('c',C_soft,S,denoise_level(n)),1);
end

out = pad(min(max(waverec2(C_soft,S,wavetype),minval),maxval));