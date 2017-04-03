function [denoised, norm_out] = bm3d_wrapper(x,snr)
sc = max(max(x));
[~,denoised] = BM3D(x/sc,x/sc,snr,'lc',0);
denoised = denoised*sc;
norm_out = 0;