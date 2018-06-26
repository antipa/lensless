%% Load test image
count = 0;
psnr_list = [];
st = 4;
xorig = double(imread('cameraman.tif'));
ds = 1/2;
x = imresize(xorig,1/ds,'box')/255;
[Ny, Nx] = size(x);
Nmax = 64;
im_list = zeros(Ny,Nx,Nmax/st);

for n = 1:st:64
    pcount = 0;
    count = count+1;
    for Np = 25
        pcount = pcount + 1;
        rng(1);
        xorig = double(imread('cameraman.tif'));

        x = imresize(xorig,1/ds,'box')/255;
        [Ny, Nx] = size(x);

        %% Loop over number of lenslets 
        %% Generate PSF
        upsamp = 8;
        Nyu = upsamp*Ny;
        Nxu = upsamp*Nx;
        % Generate random locations





        Nu = n;

        AperArea = Nu/Nmax*Nyu*Nxu;
        AperLen = floor(sqrt(AperArea)/4);
        subAper = zeros(min(AperLen,Nyu), min(AperLen,Nxu));

        if Nu>1
            ridx_random = randsample(numel(subAper),Nu-1);
        else
            ridx = []
        end
        ridx = sub2ind([AperLen,AperLen],floor(AperLen/2)+1,floor(AperLen/2)+1);
        psf = zeros(Nyu,Nxu);
        psf_gt = psf;
        subAper(ridx) = 1;
        psf_gt(ceil(Nyu/2)-floor(AperLen/2)+1:ceil(Nyu/2)+ceil(AperLen/2),ceil(Nxu/2)-floor(AperLen/2)+1:ceil(Nxu/2)+ceil(AperLen/2)) = subAper;
        if Nu>1
            ridx = cat(1,ridx,ridx_random);
        end
        subAper(ridx(2:end)) = 1;
        psf(ceil(Nyu/2)-floor(AperLen/2)+1:ceil(Nyu/2)+ceil(AperLen/2),ceil(Nxu/2)-floor(AperLen/2)+1:ceil(Nxu/2)+ceil(AperLen/2)) = subAper;
        psf_gt = psf_gt/sum(psf_gt(:));

        imagesc(psf)





        % Blur
        blur = .5;
        kernel = fspecial('gaussian',[10,10]*upsamp+1,blur*2/ds*upsamp);
        psf = conv2(psf,kernel,'same');
        psf_gt = conv2(psf_gt,kernel,'same');
        psf = imresize(psf,1/upsamp/2,'box');
        psf_gt = imresize(psf_gt,1/upsamp/2,'box');
        imagesc(psf)
        kernel2 = fspecial('gaussian',[21,21],blur*1/ds);

        x_gt = conv2(x,psf_gt/sum(psf_gt(:)),'same');

        %% Convolve and crop
        Ny = Ny/2;
        Nx = Nx/2;
        pad_len = @(x)padarray(x,[floor(Ny/2), floor(Nx/2)],'both');
        crop_len = @(x)x(floor(Ny/2)+1:end-floor(Ny/2),floor(Nx/2)+1:end-floor(Ny/2));
        psf = psf/sum(psf(:));
        H = fft2(ifftshift(pad_len(psf)));
        Hfor = @(x)crop_len(real(fftshift(ifft2(H.*fft2(ifftshift(x))))));
        b = Hfor(x);
        photon_count = Nu*Np;
        b = b*photon_count;
        b_noise = imnoise(b*1e-12,'poisson')*1e12;
        QE = .8;
        gain = .46;
        gnoise = 5;
        b_noise = round(b_noise*QE/gain)+gnoise*randn(size(b));
        imagesc(b_noise)


        %% Calculate and add noise

        %% Reconstruct (loop Tau add later)
        [xhat, f] = admm2d_solver(b_noise,[],psf,[]);
        im_list(:,:,count) = xhat;
        %% Plot
        psnr_list(count,pcount) = psnr(crop_len(xhat)/photon_count*gain/QE,crop_len(x_gt))
    end

end
%% Graduate