function [denoised, norm_out] = bm3d_wrapper(x,snr,pad)
sc = max(max(x));
[~,denoised] = BM3D(x/sc,x/sc,snr,'lc',0);
denoised = pad(gpuArray(denoised*sc));
norm_out = 0;