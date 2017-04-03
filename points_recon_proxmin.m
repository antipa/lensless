%function [im_out, mse] = points_recon_proxmin.m(A, A_adj, sensor_meas, ground_truth, options)
%% Setup
ds = 1/40;
mcount = 0;
ncount = 0;
h2 = figure(2),clf
h1 = figure(1),clf
erasure_vec = linspace(0,.99,20);
%erasure_vec = [0, .225, .45, .675, .9];
%erasure_vec = 0;
%erasure_vec = .8;
%erasure_vec = 0;
%point_vec = linspace(.01,.5,10);
point_vec = [    0.0100
    0.0358
    0.0616
    0.0874
    0.1132
    0.1389
    0.1647
    0.1905
    0.2163
    0.2421
    0.2679]';
snrvec = Inf;
%snrvec = Inf
mse_out = zeros(length(point_vec),length(erasure_vec),length(snrvec));
psf = double(imread('/Users/nick.antipa/Dropbox/Light fields without lenslets NICK/data/diffuser_flatcam/Box/psf_hdr/psf_box_exp8.tif'));
psfd = imresize(psf,ds,'box');
psfd = psfd/norm(psfd);
pad = @(x)padarray(x,[size(psfd,1)/2,size(psfd,2)/2],'both');
W = logical(pad(ones(size(psfd))));
[rr,cc] = find(W);
crop = @(x)x(min(rr):max(rr),min(cc):max(cc));

H = fft2(pad(psfd));
H_conj = conj(H);

for m = erasure_vec  %Erasure density
    mcount = mcount+1;
    ncount = 0;
    t_outer = tic;
    for n = point_vec  %Point source density
        ncount = ncount+1;
        
        erasure_density = m;
        point_density = n;
        
        Erasures = ones(size(W)/2);
        if erasure_density ~= 0
            
            %Erasures = full(sprand(size(W,1)/2,size(W,2)/2,erasure_density)>0);
            Erasure_idx = randperm(numel(W)/4,round(erasure_density*numel(W)/4));
            Erasures(Erasure_idx) = 0;
        end
        A = @(x)Erasures.*real(crop(ifftshift(ifft2(H.*fft2(x)))));
        Ah = @(x)real(ifftshift(ifft2(H_conj.*fft2(pad(Erasures.*x)))));
        rng default
        
        %im_in = pad(sprand(size(psfd,1),size(psfd,2),point_density));
        im_in = pad(full(sprand(30,40,point_vec(ncount))));
        kcount = 0;
        for k = 1:length(snrvec)
            kcount = kcount+1;
            sensor_meas = Erasures.*abs(awgn(A(im_in),snrvec(k)));
            
            
            tau = 1e-5;
            tau2 = 5e-5;
            lambda = tau;
            options.stepsize = 3e-3;
            options.convTol = 2.1e-12;
            %options.xsize = [256,256];
            options.maxIter = 10000;
            options.residTol = .00002;
            options.momentum = 'nesterov';
            options.disp_figs = 1;
            options.disp_fig_interval = 1000;   %display image this often
            options.xsize = 2*size(sensor_meas);
            %options.disp_crop = @(x)(crop(abs(x)));
            options.disp_crop = crop;
            options.disp_gamma = 1/2.2;
            options.known_input = 0;
            options.restarting = 1;
            options.fighandle = h1;
            % Function call
            
            gradErrHandle = @(x)linear_gradient(x,A,Ah,sensor_meas);
            prox_handle = @(x)soft_2d(pad(crop(x)),tau,0,Inf);
            %prox_handle = @(x)bound_range(crop(x),0,Inf,pad);
            
            x_init = zeros(size(sensor_meas)*2);
            [im_out, ~] = proxMin(gradErrHandle,prox_handle,x_init,sensor_meas,options);
            mse_current = sum(sum(abs(im_out-im_in).^2))/norm(im_in,'fro')^2;
            mse_out(ncount,mcount,kcount) = mse_current;
            set(0,'CurrentFigure',h2);
            subplot(1,3,1)
            imagesc(im_in), axis image, title('ground truth')
            subplot(1,3,2)
            imagesc(im_out), axis image, title(['reconstruction. MSE: ', num2str(mse_current)])
            subplot(1,3,3)
            imagesc(sensor_meas), axis image, title('measurement')
            drawnow
        end
    end
    
    toc(t_outer);
end

figure(3),clf
imagesc(fliplr(mse_out),'YData',point_vec,'XData',fliplr(1-erasure_vec));
xlabel('Fraction of pixels kept')
ylabel('point source density')

%%
a1 = load('/Users/nick.antipa/Dropbox/EE290S Project/FISTA results/SNR_sweep.mat');
figure(3),clf
subplot(2,1,1)
plot(snrvec,fliplr(squeeze(a1.mse_out)'))
legendCell = cellstr(num2str(100*flip(a1.erasure_vec)', '%.1f'));
lh = legend(legendCell)
 v = get(lh,'title');
set(v,'string','% Pixel Erasure')
xlabel('SNR')
ylabel('MSE')
ylim([0 20])
title('FISTA MSE with Pixel Erasure')

subplot(2,1,2)
plot(snrvec,fliplr(squeeze(mse_out)'))
legendCell = cellstr(num2str(100*(point_vec)', '%.0f'));
lh = legend(legendCell)
 v = get(lh,'title');
set(v,'string','Point Source Density')
xlabel('SNR')
ylabel('MSE')
ylim([0 12])
title('FISTA MSE with Image Sparsity')
