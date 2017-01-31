%Load impulse response stack, h
ds = 2;
colors = 'mono';
psf_in = imread('Y:\Diffusers''nstuff\Color_pco_2d_data\darpa_calibration.png');
%psf_in = imread('Y:\Grace\2d\psf_med.tif');
psf_demosaic = demosaic(psf_in,'bggr');
psf_in = mean(double(psf_demosaic),3);
enforce_obj_support = 0;

try   %Figure out if there's a usable GPU
    gpuDevice;
    use_gpu = 1;
catch
    use_gpu = 0;
end

object_close = 0;
if object_close
    mag = 1.00;
    %mag = 1.05;
    tform = affine2d([mag 0 0; 0 mag 0; 0 0 1]);
    width = size(psf_in,2);
    height = size(psf_in,1);
    hwarp = imwarp(psf_in,tform,'cubic');
    ulx = size(hwarp,2)/2-width/2;
    uly = size(hwarp,1)/2-height/2;
    if mag > 1
        psf_warp = imcrop(hwarp,[ulx,uly,width-1,height-1]);
    elseif mag < 1
        pad_size = size(psf_in)-size(hwarp);
        
        psf_warp = padarray(hwarp,floor(pad_size/2),'pre');
        psf_warp = padarray(psf_warp,ceil(pad_size/2),'post');
    end
else
    psf_warp = psf_in;
end

h_in = imresize(double(psf_warp),1/ds,'box');
    


switch colors
    case('mono')
        h = mean(h_in,3);
    case('all')
    case('red')
        h = h_in(:,:,1);
    case('green')
        h = h_in(:,:,2);
    case('blue')
        h = h_in(:,:,3);
end
h = h-100;
h = h/norm(h,'fro');



%%
%define problem size
NX = size(h,2);
NY = size(h,2);

%define crop and pad operators to handle 2D fft convolution
pad = @(x)padarray(x,[size(h,1)/2,size(h,2)/2],0,'both');
cc = (size(h,2)/2+1):(3*size(h,2)/2);
rc = (size(h,1)/2+1):(3*size(h,1)/2);
crop = @(x)x(rc,cc);

H = fft2(pad(h));
H_conj = conj(H);

% Define function handle for forward A(x)
A2d = @(x)real(crop(ifftshift(ifft2(H.*fft2(x)))));
Aadj_2d = @(x)real(ifftshift(ifft2(H_conj.*fft2(pad(x)))));

% Make or load sensor measurement
meas_type = 'pre_loaded';
switch lower(meas_type)
    case 'simulated'
        obj = zeros(size(h));
        obj(250,320,10) = 1;
        obj(250,400,20) = 1;
        b = A3d(obj);
    case 'measured'
        bin = double(imread('C:\Users\herbtipa\Documents\MATLAB\robin_close.png'));
        b = imresize(bin,1/4,'box');
    case 'pre_loaded'
        b = imresize(bin,1/ds,'box')-100;
end

% Define gradient handle
GradErrHandle = @(x) linear_gradient(x,A2d,Aadj_2d,b);

% Prox handle
tau = .001;
niters = 8;
minval = 0;
maxval = inf;
nopad = @(x)x;
%prox_handle = @(x)tv_2d(crop(x),tau,niters,minval,maxval,pad);
if ~enforce_obj_support
    prox_handle = @(x)tv_2d(x,tau,niters,minval,maxval,nopad);
    %prox_handle = @(x)bound_range(soft(x,.01),minval,maxval,nopad);
else
    prox_handle = @(x)tv_2d(crop(x),tau,niters,minval,maxval,pad);
    %prox_handle = @(x)bound_range(crop(x),minval,maxval,pad);
end
h1 = figure(1)
clf
drawnow
options.fighandle = h1;
if ds == 8
    options.stepsize = .0002;
elseif ds == 4
    options.stepsize = .000078*.9;
elseif ds == 2
    options.stepsize = .000014;
end
maxeig = (max(max(abs(fft2(pad(h))))))^2;
options.stepsize = .85*2/maxeig;
options.convTol = .005;
%options.xsize = [256,256];
options.maxIter =600;
options.residTol = 50;
options.momentum = 'nesterov';
options.disp_figs = 1;
options.disp_fig_interval = 5;   %display image this often
options.xsize = size(h);
nocrop = @(x)x;

options.disp_gamma = 1/1.5;
cm = round(.1*size(pad(h)));
cmx = round(.9*size(pad(h)));
options.disp_crop = @(x)(max(abs(x(cm(1):cmx(1),cm(2):cmx(2)))/max(x(:)),0)).^options.disp_gamma;
options.known_input = 0;
options.force_real = 1;
options.color_map = 'gray';
init_style = 'zero';

switch lower(init_style)
    case('zero')     
        xinit = zeros(size(h)*2);
        %[xhat, ~] = proxMin(GradErrHandle,prox_handle,zeros(size(h)*2),b,options);
    case('xhat')
        xinit = xhat;
        %[xhat, ~] = proxMin(GradErrHandle,prox_handle,xhat,b,options);
    case('wiener')
        xinit = deconvwnr(pad(b-min(b(:))),pad(h-min(min(h))),1);
        
end
if use_gpu
    H = gpuArray(H);
    H_conj = gpuArray(H_conj);
    b = gpuArray(b);
    xinit = gpuArray(xinit);
end
[xhat, ~] = proxMin(GradErrHandle,prox_handle,xinit,b,options);