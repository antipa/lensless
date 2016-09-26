%Make psf from hdr stack

folder = '/Users/nick.antipa/Documents/Diffusers/Lensless/cman_recovered/Box/psf_hdr'
files = dir(folder);

for n = 1:length(files);
    fgood(n) = ~(strcmp(files(n).name,'.')||strcmp(files(n).name,'..'));
end

files_good = files(fgood);  %Throw out '.' and '..' from list of filenames
stuff = imfinfo([folder,'/',files_good(2).name]);   %Get image metadata
nsamples = 500;  %Number of samples for linearization

ulr = 120;  %Upper left location of PSF (don't bother with zeros)
ulc = 260;
psfwidth = 900;   %Width of PSF in image;
ridx = ulr+randsample(psfwidth,nsamples);  %Randome samples over psf width
cidx = ulc+randsample(psfwidth,nsamples);
iidx = sub2ind([stuff.Height,stuff.Width],ridx,cidx);
pvals = zeros(numel(files_good),nsamples);

for n = 1:numel(files_good)  %Looks for file of form blah_exp###.tif' where ### is exposure time
    fname = files_good(n).name;   %Get filename
    exp_find = strfind(fname,'exp');   %Find string 'exp'
    dots = strfind(fname,'.tif');   %Find  file extension
    dt(n) = str2double(fname(exp_find(end)+3:dots(end)-1)); %Pull out exposure time
    im = imread([folder,'/',fname]);   %Read image
    pvals(n,:) = im(iidx);
end
[exp_sort,sort_i] = sort(dt,'ascend');
pvals = pvals(sort_i,:);


