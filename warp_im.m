function warped = warp_im(im,mag)


tform = affine2d([mag 0 0; 0 mag 0; 0 0 1]);
width = size(im,2);
height = size(im,1);
hwarp = imwarp(im,tform,'cubic');
ulx = size(hwarp,2)/2-width/2;
uly = size(hwarp,1)/2-height/2;
if mag > 1
    warped = imcrop(hwarp,[ulx,uly,width-1,height-1]);
elseif mag < 1
    pad_size = size(im)-size(hwarp);
    
    warped = padarray(hwarp,floor(pad_size/2),'pre');
    warped = padarray(warped,ceil(pad_size/2),'post');
end