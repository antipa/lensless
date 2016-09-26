function [out,norm_out] = sharpen_2d(x,rad,amt,minval,maxval,pad)
out = pad(min(max(imsharpen(x,'radius',rad,'amount',amt),minval),maxval));
norm_out = 0;
