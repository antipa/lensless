ds = 1/4;  %Amount to downsample problem. Use 1/integer.
useGpu = 1;
useDouble = 1;
%Load psf
%psf = double(imread('/Users/nick.antipa/Documents/Diffusers/Lensless/diffuser_flatcam/psf_again.tif'));
%psf = double(imread('/Users/nick.antipa/Documents/Diffusers/Lensless/diffuser_flatcam/Baffle/psf_16bit_baffle.tif'));
psf = double(imread('C:\Users\herbtipa\Dropbox\Light fields without lenslets NICK\data\diffuser_flatcam\Baffle\psf_16bit_baffle.tif'));
%psf = double(imread('/Users/nick.antipa/Documents/Diffusers/Lensless/diffuser_flatcam/Box/psf_hdr/psf_box_exp16.tif'));
%psf3 = double(imread('/Users/nick.antipa/Documents/Diffusers/Lensless/diffuser_flatcam/Box/psf_hdr/psf_box_exp24.tif'));
%psf = mean(cat(3,psf,psf2),3);
psfc = circshift(psf,[-80,0]);
psfd = imresize(psfc,ds,'box')*1000/max(psfc(:));  %Normalize and antialias downsample



pad = @(x)padarray(x,[size(psfd,1),size(psfd,2)],'both');
nopad = @(x)x;
psf_z = pad(psfd);

W = pad(ones(size(psfd)));
H = fft2(psf_z);
H_conj = conj(H);

bin = 0;
ii = find(W);
[r,c] = ind2sub(size(W),ii);
cr = min(r):max(r);
cc = min(c):max(c);
crop = @(x)x(cr,cc);
if ~bin
    crop_obj = crop;
else
    W_obj = imresize(W,1/bin,'box');
    ii = find(W_obj);
    [r,c] = ind2sub(size(W_obj),ii);
    cr = min(r):max(r);
    cc = min(c):max(c);
    crop_obj = @(x)x(cr,cc);
    psf_b = imresize(psfd,1/bin,'box');
    pad_obj = @(x)padarray(x,[size(psf_b,1),size(psf_b,2)],'both');
end
A = @(x)crop(ifftshift(ifft2(H.*fft2(x))));
data_type = 'measured';   %use  'measured' to load real image, 'cameraman' to make fake data from cameraman
im_type = 'peppers';

dct_sparse = 0;
lpf_true = 0;
window_object = 1;
add_noise = 0;
sense_compressively = 0;
cs_type = 'pepper';
switch lower(data_type)
    case('measured')
        
        obj = double(imread('C:\Users\herbtipa\Dropbox\Light fields without lenslets NICK\data\diffuser_flatcam\Baffle\baffle_hand.tif'));
        %obj2 = double(imread('/Users/nick.antipa/Documents/Diffusers/Lensless/diffuser_flatcam/Box/scene_hdr/box_scene_exp110.tif'));
        %obj3 = double(imread('/Users/nick.antipa/Documents/Diffusers/Lensless/diffuser_flatcam/Box/scene_hdr/box_scene_exp120.tif'));
        %obj4 = double(imread('/Users/nick.antipa/Documents/Diffusers/Lensless/diffuser_flatcam/Box/scene_hdr/box_scene_exp130.tif'));
        %obj5 = double(imread('/Users/nick.antipa/Documents/Diffusers/Lensless/diffuser_flatcam/Box/scene_hdr/box_scene_exp140.tif'));
        %obj = mean(cat(3,obj,obj2,obj3,obj4,obj5),3);
        %obj = double(imread('/Users/nick.antipa/Documents/Diffusers/Lensless/diffuser_flatcam/Baffle/baffle_scene.tif'));
        %obj = double(imread('/Users/nick.antipa/Documents/Diffusers/Lensless/diffuser_flatcam/cman_laptop_causticier.tif'));
        %psf = double(imread('/Users/nick.antipa/Documents/Diffusers/Lensless/diffuser_flatcam/psf_again.tif'));
        obj_r = imresize(obj,ds,'box');
    case('simulated')
        switch lower(im_type)
            case('cameraman')
                obj = imresize(double(imread('cameraman.tif')),2.5*2*ds,'bicubic')*255/255;
            case('deltas')
                %imsize = [400 400];
                %ndeltas = 10000;
                %randrows = randsample(imsize(1),ndeltas);
                %randcols = randsample(imsize(2),ndeltas);
                %randidx = sub2ind(imsize,randrows,randcols);
                %randidx = randsample(imsize(1)*imsize(2),ndeltas);
                %obj = zeros(imsize);
                %obj(randidx) = 18*rand(size(randidx));
                obj = delta([90,120],3*size(psfd))*255;
            case('circ')
                obj = imresize(double(imread('circles.png')),2);
            case('peppers')
                objin = imread('peppers.png');
                obj = double(objin(:,:,3));
                obj= imresize(obj,2,'bicubic');
                
            case('flat')
                obj = ones(3*size(psfd));
            case('usaf')
                obj_in = imresize(double(imread('/Users/nick.antipa/Downloads/USAF.tiff')),1/15,'box');
                obj = obj_in(:,:,2);
                
        end
        if mod(size(obj,1),2)
            obj = cat(1,zeros(1,size(obj,2)),obj);
        end
        if mod(size(obj,2),2)
            obj = cat(2,zeros(size(obj,1),1),obj);
        end
        x_in = padarray(obj,[size(psfd,1)*3/2-(size(obj,1)/2),size(psfd,2)*3/2-(size(obj,2)/2)],'both');
        if lpf_true==1
            lpf = (abs(H)<.01*max(abs(H(:))));
            x_lpf= (ifft2(fft2(x_in).*(1-lpf)));
        elseif lpf_true==0
            x_lpf = x_in;
            
        end
        if dct_sparse ==1
            x_lpf = idct2(hard(dct2(x_lpf),10));
        else
            x_lpf = x_lpf;
        end
        if window_object
            x_lpf = x_lpf.*W;
        end
        
        
        
        
        %x_in = padarray(obj,[size(psfd,1)*3/2-size(obj,1)/2,size(psfd,2)*3/2-size(obj,2)/2],'both');
        %x_lpf= (ifft2(fft2(x_in).*(1-lpf)));
        sig_n = .05e4;
        if add_noise
            obj_r = A(x_lpf)+sig_n*rand(size(psfd));
        else
            obj_r = A(x_lpf);
        end
        
        %sig_n = 0;
        %obj_r = double(uint8(obj_r*255/max(obj_r(:))+sig_n*rand(size(obj_r))));
    case('deltas')
end

if sense_compressively
    switch lower(cs_type)
        case('pepper')
            inds = randsample(numel(psfd),floor(.8*numel(psfd)));
            WD = ones(size(psfd));
            WD(inds) = 0;
            WD = pad(WD);
        case('tile')
            
            modr = mod(1:size(psfd,1),220/2);
            lr = 1:size(psfd,1);
            modc = mod(1:size(psfd,2),300/2);
            lc = 1:size(psfd,2);
            rsub = double(modr<160/2);
            csub = double(modc<210/2);
            WD = pad(rsub'*csub);
        case('random_lines')
            inds = randsample(size(psfd,1),floor(.7*size(psfd,1)));
            WD = W;
            WD(size(psfd,1)+sort(inds),:) = 0;
        case('regular_lines')
            inds_one = [1:10:size(psfd,1)];
            inds = 1:size(psfd,1);
            inds = inds(~(ismember(inds,inds_one)));
            WD = W;
            WD(size(psfd,1)+sort(inds),:) = 0;
    end
    
    obj_r = crop(WD.*pad(obj_r));
end

if bin
    obj_r = imresize(obj_r,1/bin,'box');
end

if useGpu
    if ~useDouble
        psfd = gpuArray(single(psfd));
        H = gpuArray(single(H));
        H_conj = gpuArray(single(H_conj));
        obj_r = gpuArray(single(obj_r));
    else
        psfd = gpuArray(double(psfd));
        H = gpuArray(double(H));
        H_conj = gpuArray(double(H_conj));
        obj_r = gpuArray(double(obj_r));
    end
    W = gpuArray(logical(W));
    cr = gpuArray(cr);
    cc = gpuArray(cc);
end
obj_r = obj_r/max(obj_r(:))*7.83e8;

%%
wavelev = 8;
wavetype = 'db9';
%[nrows,ncols] = size(im_in);
tau = 7e-7;
tau2 = tau;
lambda = tau;
options.stepsize = .0017;
options.convTol = 8.2e6;
%options.xsize = [256,256];
options.maxIter = 18000;
options.residTol = .2;
options.momentum = 'nesterov';
options.disp_figs = 1;
options.disp_fig_interval = 200;   %display image this often
options.xsize = 3*size(psfd);
options.disp_crop = @(x)gather(crop(abs(x)));
options.disp_gamma = 1/2.2;
options.known_input = 1;
options.force_real = 1;
if options.known_input
    options.crop = crop;
    filt = fspecial('gaussian',[5,5],.3);
    if useGpu
        if ~useDouble
            options.xin = single(conv2(x_lpf,filt,'same'));
        else
            options.xin = conv2(x_lpf,filt,'same');
        end
    else
        options.xin = double(conv2(x_lpf,filt,'same'));
    end
end
filt = fspecial('gaussian',[5,5],.2);

%%prox_handle = @(x)max(soft(x,.000001),0);
%prox_handle = @(x)max(W.*x,0);
%prox_handle = @(x)x;
%prox_handle = @(x)(x);
%prox_handle = @(x)max(W.*soft_dct2(x,tau),0);
%prox_handle = @(x)min(max(W.*tvdenoise(x,2/tau,6),0),255);
%prox_handle = @(x)W.*x;

niters = 8;
minval = 0;
maxval = inf;
%prox_handle = @(x)soft_dct2(crop(real(x)),tau,minval,maxval,ones(size(W)),pad);
%prox_handle = @(x)bound_range(crop(abs(x)),minval,maxval,pad);
%prox_handle = @(x)bound_range(bin_2d(x,2),minval,maxval,nopad);
%nopad = @(x)x;
prox_handle = @(x)bound_range(x,minval,maxval,nopad)
%prox_handle = @(x)soft_wavelet_2d(crop(x),wavelev,wavetype,tau,minval,maxval,pad)
%prox_handle = @(x)adaptive_soft_wvlt_2d(crop(x),wavelev,wavetype,.9,minval,maxval,pad);
 
%prox_handle = @(x)tv_2d(bin_2d(x,3),tau,niters,minval,maxval,nopad);
%prox_handle = @(x)tv_2d(crop(real(x)),tau,niters,minval,maxval,pad);
amt = .001;
rad = 1;

%prox_handle = @(x)sharpen_2d(crop(x),rad,amt,minval,maxval,pad);
%prox_handle = @(x)tv_dct_2d(crop(x),tau,tau2,6,minval,maxval,pad);
nopad = @(x)x;
%prox_handle = @(x)tv_dct_2d(crop(x),tau,tau2,6,minval,maxval,pad)
%prox_handle = @(x)soft_2d(x,tau,minval,maxval);
%prox_handle = @(x)x;
%prox_handle = @(x)tvdenoise(min(max(conv2(x,filt,'same'),0),255),2/tau,4);
%prox_handle = @(x)tvdenoise(min(max(conv2(x,filt,'same'),0),18),2/tau,4);
%prox_handle = @(x)max(W.*tvdenoise(x,2/tau,6),0);
%prox_handle = @(x)eight_bit(x);
%prox_handle = @(x)soft(x,.001);
%prox_handle = @(x)min(max(tvdenoise(x,2/tau,8),0),255);
%prox_handle =@(x)min(max(soft_wavelet_2d(x,wavelev,wavetype,tau),0),18);
%prox_handle = @(x)max(x,0);
%prox_handle = @(x)max(soft(x,.000001),0);
%prox_handle = @(x)max(x,0);
%prox_handle = @(x)nonneg(soft_wavelet_2d(permute(reshape(x,options.xsize),[4,3,2,1]),wavelev,wavetype,tau));
%prox_handle = @(x)max(tvdenoise(x,2/tau,6),0);
%prox_handle = @(x)max(soft(x,tau),0);

%prox_handle = @(x)nonneg(soft(reshape(TVdenoise4d(reshape(x,options.xsize),2/tau,8,1,0,apclip),...
%[numel(x),1]),2e5))
%x = 1+.9*cos((-floor(M/2):floor(M/2))'/10);
%x = soft_wavelet_2d(im_in,8,'db9',3);
%prox_handle = @(x)x;
nvar = 0;

%b = nonneg(A*x(:)+nvar*randn(size(A,1),1));
%b_in = load('/Users/nick.antipa/Documents/MATLAB/Output/king_leaf_b.mat');
%b = b_in.b*255/max(b_in.b(:));  %Scale to 8-bit
%b = b*255/max(b(:));
%Atb = A'*b;




%GradErrHandle = @(x) matrixError(A,A',x,b);
%norm_handle = @(x)lambda*norm(reshape(dct2(x),[numel(x),1]),1);
%norm_handle = @(x)lambda*TVnorm(x);
%norm_handle = @(x)lambda*norm(reshape(wdec2(x),wave
if sense_compressively
    Atb = ifftshift(ifft2(H_conj.*fft2(pad(obj_r))));
    GradErrHandle = @(x) gradient_from_psf(x,H,H_conj,WD,cr,cc,Atb,obj_r,1/(numel(psfd)^2));
else
    if ~bin
        Atb = ifftshift(ifft2(H_conj.*fft2(pad(obj_r))));
        
        GradErrHandle = @(x) gradient_from_psf(x,H,H_conj,W,cr,cc,Atb,obj_r,1/(numel(psfd)^2));
    else
        filt = ones(bin)/bin^2;
        Atb = ifftshift(ifft2(H_conj.*fft2(pad(imresize(obj_r,bin,'nearest')))));
        GradErrHandle = @(x) gradient_from_psf_bin(x,H,H_conj,W_obj,W,crop_obj,Atb,obj_r,1/(numel(psfd)^2),filt,1/bin);
    end
end



%norm_hanlde = @(x)TVnorm(x);


%x_init = deconvwnr(pad(obj_r),psfd,1e3).*W;
%x_init = xhat;%rand(size(Atb));
%x_init = nspace;
%x_init = zeros(size(xhat));
%x_init = zeros(size(pad(psfd)));\
if useGpu
    if ~useDouble
        x_init = gpuArray(single(zeros(size(Atb))));
    elseif useDouble
        x_init = gpuArray(zeros(size(Atb)));
        %x_init = xhat;
    end
    %x_init = xhat;
else
    x_init = zeros(size(Atb));
end
%x_init = x_lpf+randn(size(x_lpf));
%x_init = 255*rand(size(Atb));
%x_init = W.*(xhat+.07*rand(size(xhat)).*max(xhat(:)));
figure(1),clf
[xhat, funvals] = proxMin(GradErrHandle,prox_handle,x_init,obj_r,options);


