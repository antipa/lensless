function out = bin_2d(im,n)
out = imresize(imresize(im,1/n,'box'),n,'nearest');
