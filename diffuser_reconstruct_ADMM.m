%psf = imread('/Users/nick.antipa/Documents/Diffusers/Lensless/pco_2d_color/psf_med.tif');
%psfin = load('/Volumes/PhantomData/Diffusers''nstuff/colorPCO_collimated_calibration_far/pinhole_far_green_coherent.mat');

%psf = psfin.imave;
ds = 4;
demos = 0;
random_erasure = .1;   %If enabled, randomly delete pixels before recon

if demos
    psf = demosaic(imread('/Users/nick.antipa/Documents/pinhole_far_green_00013.png'),'rggb');
    psf_in = psf(:,:,2)-100;
else
    psf = imread('/Users/nick.antipa/Documents/pinhole_far_green_00013.png');
    psf_in = double(psf(2:2:end,1:2:end))/2 + double(psf(1:2:end,2:2:end))/2-100;
    ds = ds/2;
end
%%
color_proc = 'rgb'

switch lower(color_proc)
    case('red')
        color_vec = 3;
    case('green')
        color_vec = 2;
    case('blue')
        color_vec = 1;
    case('rgb')
        color_vec = 1:3;
end
      
ccount = 0;

for color = color_vec
    object_close = 1;
    ccount = ccount + 1;

    if object_close
        %mag = .9744;
        %mag = 1.02;
        mag = 1.04;
        %mag = 1.005;
        %mag = .98;%1.01;
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

    if ds ~= 1
        psf_ = imresize(double(psf),ds,'box');
    else
        psf_in = psf;
    end
    psfd = imresize(double(psf_warp),1/ds,'box');


    h = psfd/norm(psfd(:),2);
    Nx = size(h,2);
    Ny = size(h,1);
    if ccount == 1
        im_out = zeros(2*Ny,2*Nx,length(color_vec));
        y_out = zeros(Ny,Nx,length(color_vec));
    end
    %%
    imin = imread('cameraman.tif');
    %imin = imread('/Users/nick.antipa/Documents/Diffusers/Lensless/tapeycam/tape_40cm_shrey1.png');
    %imin = imread('/Users/nick.antipa/Documents/Diffusers/Lensless/tapeycam/tape_40cm_hand1.png');

    if size(imin,1) ~= size(psfd,1)
        im = (double(imresize(imin,2*[Ny,Nx],'box')));
    else
        im = double(imin)/255;
    end
    h1 = figure(1),clf;

    %define crop and pad operators to handle 2D fft convolution
    rng(0)
    if random_erasure
        erasure_window = full(double(sprand(Ny,Nx,random_erasure)>0));
    else
        erasure_window = ones(Ny,Nx);
    end
    pad = @(x)padarray(x,[Ny/2,Nx/2],0,'both');
    cc = (Nx/2+1):(3*Nx/2);
    rc = (Ny/2+1):(3*Ny/2);
    crop = @(x)x(rc,cc);

    H = fft2(ifftshift(pad(h)));
    H_conj = conj(H);

    % Define forward and adjoints that would be used in gradient descent
    A = @(x)erasure_window.*crop(ifftshift(real(ifft2(H.*fft2(ifftshift(x))))));
    At = @(b)real(ifftshift(ifft2(H_conj.*fft2(ifftshift(pad(erasure_window.*b))))));


    %y = A(im);
    %yin = imread('/Volumes/PhantomData/Diffusers''nstuff/Processing_results/20170419_152324_AVG_Video_rgb/AVG_Video.png');
    %yin = demosaic(yin,'bggr');
    %yin = double(imread('/Users/nick.antipa/Documents/Diffusers/Lensless/tapeycam/tape_40cm_shrey1.png'));
    %yin = double(imread('/Users/nick.antipa/Documents/Diffusers/Lensless/tapeycam/tape_40cm_hand1.png'));
    %%
    if demos
        yinb = demosaic(imread('/Users/nick.antipa/Documents/diffuser_cam_results/grace_dog_close_scale_102.png'),'rggb');
        yin = double(yinb(:,:,color))-100;
    else
        yinb = imread('/Users/nick.antipa/Documents/diffuser_cam_results/grace_dog_close_scale_102.png');
        yinred = double(yinb(1:2:end,1:2:end))-100;
        yingreen = double(yinb(2:2:end,1:2:end))/2 + double(yinb(1:2:end,2:2:end))/2-100;
        yinblue = double(yinb(2:2:end,2:2:end))-100;
        if color == 1
            yin = yinred;
        elseif color == 2
            yin = yingreen;
        elseif color ==3
            yin = yinblue;
        end
    end
    y = imresize(yin,1/ds,'box');

    y = erasure_window.*y/(2^16);
    y_out(:,:,ccount) = y;
    %if size(y)~=size(h)
     %   y = imresize(y,[Ny,Nx],'box');
    %end
    imagesc(y), axis image

    %%


    s0 = zeros;

    % Adjoint of TV finite difference
    Ltv = @(P1,P2)cat(1,P1(1,:),diff(P1,1,1),-P1(end,:)) + cat(2,P2(:,1),diff(P2,1,2),-P2(:,end));

    % Forward finite difference (returns two images!) [P1, P2] = L(x)
    L = @(D)deal(-diff(D,1,1),-diff(D,1,2));

    % Uncropped forward
    Hfor = @(x)real(ifftshift(ifft2(H.*fft2(ifftshift(x)))));
    Hadj = @(x)real(ifftshift(ifft2(H_conj.*fft2(ifftshift(x)))));

    % Precompute spectrum of L'L
    d = zeros(Ny,Nx);
    d(floor(Ny/2),floor(Nx/2)) = 1;
    % [P1,P2] = L(d);
    % lapl = Ltv(P1,P2);
    % LtL = real(psf2otf(lapl,2*[Ny Nx]));
    lapl = zeros(2*Ny,2*Nx);
    lapl(1) = 4;
    lapl(1,2) = -1;
    lapl(2,1) = -1;
    lapl(1,end) = -1;
    lapl(end,1) = -1;
    LtL = abs(fft2(lapl));
    imagesc(abs(LtL));

    %Precompute H'H
    %psfacorr = Hadj(Hfor(pad(d)));
    HtH = H.*H_conj;
    imagesc(abs(HtH));

    mu1 = 0.020972;
    mu2 = 6.815744;
    mu3 = 1.363149;
    tau = .0003;

    % mu1 = 1e-1;
    % mu2 = 1;
    % mu3 = 4;
    % tau = .1;
    tauinc = 2;
    taudec = tauinc;
    muinc = 10;

    I(1,1) = 1;
    Smult = 1./(mu1*HtH + mu2*LtL + mu3);  %This is the frequency space division fo S update
    if find(isinf(Smult))
        error('Dividing by zero!')
    end
    imagesc(abs(Smult))
    colorbar
    % Compute V denominator
    CtC = pad(erasure_window);
    Vmult = 1./(CtC + mu1);


    % Initialize variables
    %sk = im;
    sk = zeros(Ny*2,Nx*2);
    alpha2k_1 = zeros(2*Ny-1,2*Nx);
    alpha2k_1p = alpha2k_1;
    alpha2k_2 = zeros(2*Ny,2*Nx-1);
    alpha2k_2p = alpha2k_2;
    alpha1k = zeros(2*Ny,2*Nx);
    alpha1kp = alpha1k;
    alpha3k = zeros(2*Ny,2*Nx);
    alpha3kp = alpha3k;
    Cty = pad(y);
    tk = 1;
    tv_iso = 1;

    n = 0;
    niter = 200;
    while n<niter
        n = n+1;
        [Lsk1,Lsk2] = L(sk);
        uk_1 = Lsk1+alpha2k_1/mu2;
        uk_2 = Lsk2+alpha2k_2/mu2;
        if tv_iso
            mag = sqrt(cat(1,uk_1,zeros(1,size(uk_1,2))).^2 + ...
            cat(2,uk_2,zeros(size(uk_2,1),1)).^2);
            magt = soft(mag,tau/mu2);
            mag(mag==0) = eps;
            mmult = magt./mag; 
            mmult(mag==0) = 0;
            ukp_1 = uk_1.*mmult(1:end-1,:);
            ukp_2 = uk_2.*mmult(:,1:end-1);
        else
            ukp_1 = soft(uk_1,tau/mu2);
            ukp_2 = soft(uk_2,tau/mu2);
        end

        vkp = Vmult.*(mu1*(alpha1k/mu1 + Hfor(sk)) + Cty);
        wkp = max(alpha3k/mu3 + sk,0);

        skp = real(ifftshift(ifft2(Smult .* fft2(ifftshift(mu3*(wkp-alpha3k/mu3) + mu2*Ltv(ukp_1 - alpha2k_1/mu2,ukp_2 - alpha2k_2/mu2) + mu1*Hadj(vkp - alpha1k/mu1))))));

        t_kp1 = (1+sqrt(1+4*tk^2))/2;
        beta_kp1 = (tk-1)/t_kp1;

        alpha1kp = alpha1k + mu1*(Hfor(skp) - vkp);
        alpha1kp = alpha1kp + beta_kp1*(alpha1kp-alpha1k);
        ea1 = norm(alpha1kp-alpha1k,'fro');
        alpha1k = alpha1kp;
        [Lskp1, Lskp2] = L(skp);
        alpha2k_1p = alpha2k_1 + mu2*(Lskp1 - ukp_1) + beta_kp1*(alpha2k_1p - alpha2k_1);
        alpha2k_2p = alpha2k_2 + mu2*(Lskp2 - ukp_2) + beta_kp1*(alpha2k_2p - alpha2k_2);
        alpha2k_1p = alpha2k_1p + beta_kp1*(alpha2k_1p-alpha2k_1);
        alpha2k_2p = alpha2k_2p + beta_kp1*(alpha2k_2p-alpha2k_2);
        ea2a = norm(alpha2k_1p-alpha2k_1,'fro');
        ea2b = norm(alpha2k_2p-alpha2k_2,'fro');
        alpha2k_1 = alpha2k_1p;
        alpha2k_2 = alpha2k_2p;
        alpha3kp = alpha3k + mu3*(skp - wkp) + beta_kp1*(alpha3kp-alpha3k);
        alpha3kp = alpha3kp + beta_kp1*(alpha3kp-alpha3k);
        ea3 = norm(alpha3kp-alpha3k,'fro');
        alpha3k = alpha3kp;
        e_primal = norm(sk-skp,'fro');
        sk = skp;
        tk = t_kp1;

        if mod(n,10)==0
            set(0,'CurrentFigure',h1);
            f = norm(A(sk)-y,'fro') + tau*TVnorm(sk);
            fprintf('%.4f\t %.4f\t %.4f\t %.4f\t %.4f\t %.4f\n',ea1,ea2a,ea2b,ea3,e_primal,f)
            %nrm = TVnorm(sk) + 
            imagesc(skp);
            axis image

            caxis([prctile(skp(:),.001) prctile(skp(:),99.99)])
           % axis(ds/.1*[128.5166  356.6414  107.6861  300.1663]);
            colorbar
            colormap gray
            drawnow
        end

    end
    im_out(:,:,ccount) = sk;
end
%%
imagesc(im_out/max(im_out(:)))
