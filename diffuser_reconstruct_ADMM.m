%psf = imread('/Users/nick.antipa/Documents/Diffusers/Lensless/pco_2d_color/psf_med.tif');
%psfin = load('/Volumes/PhantomData/Diffusers''nstuff/colorPCO_collimated_calibration_far/pinhole_far_green_coherent.mat');

%psf = psfin.imave;
ds = 4;
tic
demos = 0;
hsc=50;
regularizer_type = 'tv';
camera_type = 'pco';
erasure_type = 'none' %Choices: random, tile, none
object_close = 0;   %Magnify PSF to account for object distance changes
color_correct = 0;
psf_color_proc = 'rgb';
mag_vec = .95;  %Adjust 1PSF distance
color_proc = 'rgb';
random_erasure = 3/16;   %If enabled, randomly delete pixels before recon. Use a number between 0 and 1 to do erasure, and 0 for full res
y_tiles = 4;
x_tiles = 4;
tile_gap = .5;  %Fraction of individual sensor width to leave as a gap
% ADMM stuff
%tau_vec = [8e-4 4e-4 10e-4]/hsc;  %TV params for BSI
tau_vec = [1e-3 1e-3 5e-3]/.01;  %TV params for PCO
taumax = 10;
taumin = 2;

niter = 500;
tauiter = logspace(log10(taumax),log10(taumin),niter);
tau2 = 20e-4;  %DCT


quantize_error = 1;
B = 8;
noise_std = 0/(2^B-1);  %fraction of max, only used when simulating
mu1 = .09;
mu2 = 5.9;
mu3 = 5;
mu4 = .7; %DCT

mu_inc = 1.2;
mu_dec = 1.2;
resid_tol = 1.5;
if demos
    
    switch lower(camera_type)
        case('pco')
            psf = demosaic(imread('/Users/nick.antipa/Documents/pinhole_far_green_00013.png'),'bggr');
            psf_in = psf(:,:,2)-100;
        case('ids_bsi')
            psf_in = double(demosaic(imread('/Users/nick.antipa/Documents/Diffusers/Lensless/BSI_color/2p5ft_LED_12bit_raw.png'),'rggb'));
            [Ny_in, Nx_in,~] = size(psf_in);
            c_crop = 1:(ds*2*floor(Nx_in/ds/2));
            r_crop = 1:(ds*2*floor(Ny_in/ds/2));
            psf_in = psf_in(r_crop,c_crop,:)-128;

    end
else
    imbias = 100;
    %psf = imread('/Users/nick.antipa/Documents/pinhole_far_green_00013.png');
    switch lower(camera_type)
        case('pco')
            psf_in = load('/Users/nick.antipa/Documents/Diffusers/Lensless/pco_2d_color/pinhole_far_green_coherent.mat');
            %psf = double(imread('/Volumes/PhantomData/Diffusers''nstuff/Big_DiffuserMicroscope/171201/789_00001.tif'));
            psf = psf_in.imave;
            %psf_in = load('/Users/nick.antipa/Documents/developer/lensless/poisson_lenslets_best.mat');
            %psf_in = load('/Users/nick.antipa/Documents/developer/lensless/poisson_lenslets_best_sphere_k2000.mat');
            %psf_in = load('./uniform_lenslets_best_moresphere_k300.mat');
            %psf_in = load('/Users/nick.antipa/Documents/developer/lensless/poisson_lenslets_worst_k110e6.mat');
            %psf_in = load('/Users/nick.antipa/Documents/developer/lensless/poisson_lenslets_best_k2600.mat');
            %psf_in = load('./uniform_miniscope_k5k_f36.mat');
            %psf_in = load('/Users/nick.antipa/Documents/developer/lensless/poisson_lenslets_best_k2600.mat');
            %psf_in = load('/Users/nick.antipa/Documents/Diffusers/Miniscope/Masks/uniform_quantized_H_498e3_correct_f=10.mat');
            %psf_in = load('/Users/nick.antipa/Documents/Diffusers/Miniscope/Masks/uniform_quantized_H_330e3_225lenslets_f=10.mat');
            %psf_in = load('/Users/nick.antipa/Documents/Diffusers/Miniscope/Masks/regular_quantized_H_111e6_225lenslets_f=10.mat');
            %psf_in = load('/Users/nick.antipa/Documents/Diffusers/Miniscope/Masks/uniform_quantized_H_6e9_225lenslets_f=10.mat');
            %psf_in = load('/Users/nick.antipa/Documents/Diffusers/Miniscope/Masks/poisson_quantized_sidelobes_4_225lenslets_f=10.mat');
            %psf = imread('/Users/nick.antipa/Documents/Diffusers/Miniscope/Masks/Zstacks_lenslets_longer/lenslets_longer_ 40.png');
            %psf_in = load('/Users/nick.antipa/Documents/Diffusers/Miniscope/Masks/gaussian_diffuser_f=10.mat');
            %psf_in = load('/Users/nick.antipa/Documents/Diffusers/Miniscope/Masks/uniform_quantized_H_7e9_correct_f=10.mat');
            %psf = imread('/Users/nick.antipa/Documents/Diffusers/Miniscope/Masks/zstack_diffuser/diffuser_  0.png');
            %psf_in = load('/Users/nick.antipa/Documents/Diffusers/Miniscope/Masks/regular_quantized_sidelobes_1p15_K_2239e6.mat')
            %psf_in = load('/Users/nick.antipa/Documents/Diffusers/Miniscope/Masks/uniform_quantized_H_3899.mat');
            %psf_in = load('./uniform_miniscope_k136k_f36_bad_fillfactor.mat');
            %psf_in = load('poisson_miniscope_k20k_f36.mat');
            %psf_in = load('./uniform_miniscope_kbad_f36.mat');
            %psf = padarray(psf_in.Id,[200,200],'both');
            %psf = psf_in.Id;
           % psf = padarray(psf_in.Ibest,[00,00],'both');
           %psf = psf_in.Id(115:366,255:506);
            %psf = padarray(psf,[206,206],'both');
            %psf = psf_in.IworstF;
            psf_in = double(psf(2:2:end,1:2:end))/2 + double(psf(1:2:end,2:2:end))/2-imbias;
            %psf_in = double(psf_in(1:end-1,1:end-1));
            %ds = ds/2;
        case('ids_bsi')
            psf_in = double(imread('/Users/nick.antipa/Documents/Diffusers/Lensless/BSI_color/LED_2ft_hotpix_corrected.png'));
            [Ny_in, Nx_in,~] = size(psf_in);
            c_crop = 1:(ds*2*floor(Nx_in/ds/2));
            r_crop = 1:(ds*2*floor(Ny_in/ds/2));
            psf_in = psf_in(r_crop,c_crop,:);
        case('grasshopper')
            psf_in = double(imread('/Users/nick.antipa/Documents/psfasdf.png'));
            [Ny_in, Nx_in,~] = size(psf_in);
            c_crop = 1:(ds*2*floor(Nx_in/ds/2));
            r_crop = 1:(ds*2*floor(Ny_in/ds/2));
            psf_in = psf_in(r_crop,c_crop,:);
            
    end
end
%%


switch lower(color_proc)
    case('red')
        color_vec = 3;
    case('green')
        color_vec = 2;
    case('blue')
        color_vec = 1;
    case('rgb')
        color_vec = 1:3;
    case('mono')
        color_vec = 1;
end


%mag_vec = 1.04;

mag_count = 0;
bigT = tic
for mag_set = mag_vec
    mag_count = mag_count + 1;
    if mag_count == 1
        im_stack = cell(length(mag_vec),1);
    end
    ccount = 0;
    for color = color_vec
        tau = tau_vec(color);
        ccount = ccount + 1;
        switch lower(camera_type)
            case('pco')
                if object_close
                    psf_warp = warp_im(psf_in,mag_set);
                else
                    psf_warp = psf_in;
                end
            case{'ids_bsi','grasshopper'}
                if strcmpi(psf_color_proc,'mono')
                    psf_warp = warp_im(mean(double(psf_in),3),mag_set);
                else
                    if object_close
                        psf_warp = warp_im(psf_in(:,:,color),mag_set);
                    else
                        psf_warp = psf_in(:,:,color);
                    end
                end
         
        end
        
        psfd = imresize(double(psf_warp),1/ds,'box');
        
        
        %h = psfd/norm(psfd(:),2)*2^8;
        h = psfd/max(psfd(:)); 
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
        h1 = figure(2),clf;
        
        %define crop and pad operators to handle 2D fft convolution
        rng(0)
        switch lower(erasure_type)
            case('random')
                erasure_window = full(double(sprand(Ny,Nx,random_erasure)>0));
            case('none')
                erasure_window = ones(Ny,Nx);
            case('tile')
                
                x_tile_size = floor(Nx/(1+(x_tiles-1)*(1+tile_gap)));
                y_tile_size = floor(Ny/(1+(y_tiles-1)*(1+tile_gap)));
                erasure_window_x = zeros(1,Nx);
                erasure_window_y = zeros(Ny,1);
                for ewx = 1:x_tiles
                    xstart = floor(1+(ewx-1)*x_tile_size*(1+tile_gap));
                    ystart = floor(1+(ewx-1)*y_tile_size*(1+tile_gap));
                    erasure_window_x(xstart:xstart+x_tile_size-1) = 1;
                    erasure_window_y(ystart:ystart+y_tile_size-1) = 1;
                end
                erasure_window = erasure_window_y*erasure_window_x;
        end
        pad = @(x)padarray(x,[Ny/2,Nx/2],0,'both');
        cc = (Nx/2+1):(3*Nx/2);
        rc = (Ny/2+1):(3*Ny/2);
        crop = @(x)x(rc,cc);
        
        H = (fft2(ifftshift(pad(h))));
        H_conj = conj(H);
        
        % Define forward and adjoints that would be used in gradient descent
        A = @(x)erasure_window.*crop(real(ifft2(H.*fft2(x))));
        At = @(b)real(ifft2(H_conj.*fft2(pad(erasure_window.*b))));
        
        
        %y = A(im);
        %yin = imread('/Volumes/PhantomData/Diffusers''nstuff/Processing_results/20170419_152324_AVG_Video_rgb/AVG_Video.png');
        %yin = demosaic(yin,'bggr');
        %yin = double(imread('/Users/nick.antipa/Documents/Diffusers/Lensless/tapeycam/tape_40cm_shrey1.png'));
        %yin = double(imread('/Users/nick.antipa/Documents/Diffusers/Lensless/tapeycam/tape_40cm_hand1.png'));
        %%
        if demos
            switch lower(camera_type)
                case('pco')
                    yinb = demosaic(imread('/Users/nick.antipa/Documents/diffuser_cam_results/grace_dog_close_scale_102.png'),'rggb');
                    yin = double(yinb(:,:,color))-100;
                case('ids_bsi')
                    yinb = double(demosaic(imread('/Users/nick.antipa/Documents/Diffusers/Lensless/BSI_color/nick3_8bit_raw_nowb.png'),'rggb'));
                    yin = yinb(r_crop,c_crop,:)-4;
                    yinred = yin(r_crop,c_crop,1);
                    yingreen = yin(r_crop,c_crop,2);
                    yinblue = yin(r_crop,c_crop,3);
            end
        else
            
            color_calibration = load('/Users/nick.antipa/Documents/Diffusers/Lensless/pco_2d_color/calibration/bayer_filter_coefficients.mat');
            %yinb = double(imread('/Users/nick.antipa/Documents/Diffusers/Lensless/pco_2d_color/laura_moreTV_540.png'));
            %yinb = double(imread('/Users/nick.antipa/Documents/diffuser_cam_results/gra'ce_dog_close_scale_102.png'));
            %yinb = double(imread('/Users/nick.antipa/Google Drive/dcam/test_012.png'));
            %yinb = double(imread('/Volumes/PhantomData/Diffusers''nstuff/Big_DiffuserMicroscope/171201/ghfh_00011.tif'));
            %yinb = double(imread('/Users/nick.antipa/Documents/Diffusers/Lensless/BSI_color/nick2_test_no_awb.png'));
            %yinb = double(imread('/Users/nick.antipa/Documents/Diffusers/Lensless/pco_2d_color/lab_scene_ave-raw.png'))-imbias;
            yinb = double(imread('/Users/nick.antipa/Documents/Diffusers/Diffuser data/emrah_pco_rgb.png'));
            %yinb = A(im);
           
            if numel(size(yinb)) == 2
                if demos
                    switch lower(camera_type)
                        case('pco')
                            yinred = double(yinb(1:2:end,1:2:end));
                            yingreen = double(yinb(2:2:end,1:2:end))/2 + double(yinb(1:2:end,2:2:end))/2;
                            yinblue = double(yinb(2:2:end,2:2:end));
                        case{'ids_bsi','grasshopper'}
                            yinred = double(yinb(1:2:end,1:2:end));
                            yingreen = double(yinb(2:2:end,1:2:end))/2 + double(yinb(1:2:end,2:2:end))/2;
                            yinblue = double(yinb(2:2:end,2:2:end));
                            
                    end
                else
                    yinred = double(yinb);
                    yingreen = yinred;
                    yinblue = yinred;
                end
            else
                switch lower(camera_type)
                    case('pco')
                        yinred = imresize(double(yinb(:,:,1)),.5,'box')-100;
                        yingreen = imresize(double(yinb(:,:,2)),.5,'box')-100;
                        yinblue = imresize(double(yinb(:,:,3)),.5,'box')-100;
                    case{'ids_bsi','grasshopper'}
                        %yinred = yinb(r_crop,c_crop,1);
                        %yingreen = yinb(r_crop,c_crop,2);
                        %yinblue = yinb(r_crop,c_crop,3);
                        yinred = yinb(:,:,1);
                        yingreen = yinb(:,:,2);
                        yinblue = yinb(:,:,3);
                end
            end
            if color_correct
                Ideconv = PCO_demosaic(yinb,color_calibration.couplingmatrix);
                %Ideconv = Ideconv/max(Ideconv(:)) * max(yinb(:) - 100);
                yinred = Ideconv(:,:,1)-100;
                yingreen = Ideconv(:,:,2)-100;
                yinblue = Ideconv(:,:,3)-100;
            end

            

            
          
        end
        if color == 1
            yin = yinred;
        elseif color == 2
            yin = yingreen;
        elseif color ==3
            yin = yinblue;
        end
        if strcmpi(color_proc,'mono')
            yin = (yinred+yingreen+yinblue)/3;
        end
        y = imresize(yin,[Ny,Nx],'box');
        
        y = erasure_window.*y/max(yin(:))*255;
        y = y+randn(size(y))*noise_std;
        if quantize_error
            y = double(floor(y*2^B))*2^(-B);
        end
        y_out(:,:,ccount) = y;
        %if size(y)~=size(h)
        %   y = imresize(y,[Ny,Nx],'box');
        %end
        imagesc(y), axis image
        
        %%
        
        
        s0 = zeros;
        
        % Adjoint of TV finite difference
        % Precompute spectrum of L'L
        d = zeros(Ny,Nx);
        d(floor(Ny/2),floor(Nx/2)) = 1;
        % Uncropped forward
        Hfor = @(x)real(ifft2(H.*fft2(x)));
        Hadj = @(x)real(ifft2(H_conj.*fft2(x)));
        % Initialize variables
        sk = zeros(Ny*2,Nx*2);
        alpha1k = zeros(2*Ny,2*Nx);
        alpha3k = zeros(2*Ny,2*Nx);
        alpha3kp = alpha3k;
        Cty = pad(y);
        tk = 1;
        tv_iso = 1;
        
        % TV stuff
        Ltv = @(P1,P2)cat(1,P1(1,:),diff(P1,1,1),-P1(end,:)) + cat(2,P2(:,1),diff(P2,1,2),-P2(:,end));
        
        % Forward finite difference (returns two images!) [P1, P2] = L(x)
        L = @(D)deal(-diff(D,1,1),-diff(D,1,2));
        
        lapl = zeros(2*Ny,2*Nx);
        lapl(1) = 4;
        lapl(1,2) = -1;
        lapl(2,1) = -1;
        lapl(1,end) = -1;
        lapl(end,1) = -1;
        LtL = abs(fft2(lapl));
        alpha2k_1 = zeros(2*Ny-1,2*Nx);
        alpha2k_2 = zeros(2*Ny,2*Nx-1);
        
        
        
        % DCT stuff
        DCT_adj = @(Q)idct2(Q);
        DCT_for = @(D)dct2(D);
        
        %LtL = ones(2*Ny,2*Nx);%abs(fft2(dct2(pad(d))));
        alpha4k = zeros(2*Ny,2*Nx);
        alpha4kp = alpha4k;
        
        %Precompute H'H
        %psfacorr = Hadj(Hfor(pad(d)));
        HtH = H.*H_conj;
        imagesc(abs(HtH));
        
        
        I(1,1) = 1;
        switch lower(regularizer_type)
            case('tv')
                Smult = 1./(mu1*HtH + mu2*LtL + mu3);  %This is the frequency space division fo S update
            case('dct')
                Smult = 1./(mu1*HtH + mu3 + mu4);
            case('both')
                Smult = 1./(mu1*HtH + mu2*LtL + mu3 + mu4);
        end
        if find(isinf(Smult))
            error('Dividing by zero!')
        end
        imagesc(abs(Smult))
        colorbar
        % Compute V denominator
        CtC = pad(erasure_window);
        Vmult = 1./(CtC + mu1);
        Hskp = zeros(size(Smult));
        
        
        
        n = 0;

        dual_resid_s = zeros(1,niter)./0;
        primal_resid_s = zeros(1,niter)./0;
        dual_resid_u = dual_resid_s;
        primal_resid_u = dual_resid_u;
        dual_resid_w = dual_resid_s;
        primal_resid_w = dual_resid_s;
        
        f = primal_resid_u;
        [ukp_1, ukp_2] = L(zeros(2*Ny, 2*Nx));
        Lsk1 = ukp_1;
        Lsk2 = ukp_2;
        %----------------- Main Loop 0----------------------
        
        while n<niter
            n = n+1;
            tau = tauiter(n);
            Lskm1 = Lsk1;
            Lskm2 = Lsk2;
            [Lsk1,Lsk2] = L(sk);
            
            uk_1 = ukp_1;
            uk_2 = ukp_2;
            [ukp_1, ukp_2] = soft_2d_gradient(Lsk1+alpha2k_1/mu2, Lsk2+alpha2k_2/mu2, tau/mu2);
            Hsk = Hskp;
            vkp = Vmult.*(mu1*(alpha1k/mu1 + Hsk) + Cty);
            wkp = max(alpha3k/mu3 + sk,0);
            switch lower(regularizer_type)
                case('tv')
                    skp_numerator = mu3*(wkp-alpha3k/mu3) + ...
                        mu2*Ltv(ukp_1 - alpha2k_1/mu2,ukp_2 - alpha2k_2/mu2) + ...
                        mu1*Hadj(vkp - alpha1k/mu1);
                case('dct')
                    zkp = soft(DCT_for(sk) + alpha4k/mu4,tau2/mu4);
                    skp_numerator = mu3*(wkp-alpha3k/mu3) + ...
                        DCT_adj(mu4*zkp - alpha4k) + ...
                        mu1*Hadj(vkp - alpha1k/mu1);
                case('both')
                    zkp = soft(DCT_for(sk) + alpha4k/mu4,tau2/mu4);
                    skp_numerator = mu3*(wkp-alpha3k/mu3) + ...
                        mu2*Ltv(ukp_1 - alpha2k_1/mu2,ukp_2 - alpha2k_2/mu2) + ...
                        mu1*Hadj(vkp - alpha1k/mu1) + ...
                        DCT_adj(mu4*zkp - alpha4k);
            end
            
            skp = real(ifft2(Smult .* fft2(skp_numerator)));
            
            %Update dual and parameter for Hs=v constraint
            Hskp = Hfor(skp);
            r_sv = Hskp-vkp;
            alpha1k = alpha1k + mu1*r_sv;
            %dual_resid_s(n) = mu1*norm(Hfor(sk-skp),'fro');
            dual_resid_s(n) = mu1*norm(Hskp - Hsk,'fro');
            primal_resid_s(n) = norm(r_sv);
            [mu1, mu1_update] = update_param(mu1,resid_tol,mu_inc,mu_dec,primal_resid_s(n),dual_resid_s(n));
            
            % Update dual and parameter for Ls=v
            [Lskp1, Lskp2] = L(skp);
            r_su_1 = Lskp1 - ukp_1;
            r_su_2 = Lskp2 - ukp_2;
            alpha2k_1 = alpha2k_1 + mu2*r_su_1;
            alpha2k_2 = alpha2k_2 + mu2*r_su_2;
            dual_resid_u(n) = mu2*sqrt(norm(Lsk1 - Lskm1,'fro')^2+norm(Lsk2 - Lskm2,'fro')^2);
            primal_resid_u(n) = sqrt(norm(r_su_1,'fro')^2 + norm(r_su_2,'fro')^2);
            [mu2, mu2_update] = update_param(mu2,resid_tol,mu_inc,mu_dec,primal_resid_u(n),dual_resid_u(n));
            
            % Update nonnegativity dual and parameter (s=w)
            r_sw = skp-wkp;
            alpha3k = alpha3k + mu3*r_sw;
            dual_resid_w(n) = mu3*norm(sk - skp,'fro');
            primal_resid_w(n) = norm(r_sw,'fro');
            [mu3, mu3_update] = update_param(mu3,resid_tol,mu_inc,mu_dec,primal_resid_w(n),dual_resid_w(n));
            
            if strcmpi(regularizer_type,'dct') || strcmpi(regularizer_type,'both')
                
                alpha4k = alpha4k + mu4*(DCT_for(skp) - zkp);
            end
            
            %Update filters
            if mu1_update || mu2_update || mu3_update
                %fprintf('Mu updates: %i \t %i \t %i\n',mu1_update, mu2_update, mu3_update);
                mu_update = 1;
            end
            if mu_update
                switch lower(regularizer_type)
                    case('tv')
                        Smult = 1./(mu1*HtH + mu2*LtL + mu3);  %This is the frequency space division fo S update
                    case('dct')
                        Smult = 1./(mu1*HtH + mu3 + mu4);
                    case('both')
                        Smult = 1./(mu1*HtH + mu2*LtL + mu3 + mu4);
                end
                Vmult = 1./(CtC + mu1);
            end
            
            
            %alpha4kp = alpha4kp + beta_kp1*(alpha4kp - alpha4k);
            
            sk = skp;
            f(n) = norm(crop(Hskp)-y,'fro')^2 + tau*TVnorm(skp);
            if exist('im')
                sig_to_noise(n) = psnr(crop(skp),crop(im)/max(yinb(:)),255/max(yinb(:)));
            end
            if mod(n,20)==0
                toc
                fprintf('iter: %i \t cost: %.4f \t Primal v: %.4f \t Dual v: %.4f \t Primal u: %.4f \t Dual u: %.4f \t Primal w: %.4f \t Dual w: %.4f \t mu1: %.4f \t mu2: %.4f \t mu3: %.4f \n',...
                    n,f(n),primal_resid_s(n), dual_resid_s(n),primal_resid_u(n), dual_resid_u(n),primal_resid_w(n), dual_resid_w(n),mu1,mu2,mu3)
                
                set(0,'CurrentFigure',h1);
                
                %fprintf('%.4f\t %.4f\t %.4f\t %.4f\t %.4f\t %.4f\n',ea1,ea2a,ea2b,ea3,e_primal,f)
                %fprintf('%.4f\n',norm(Hfor(skp) - vkp,'fro'))
                %nrm = TVnorm(sk) +
                %subplot(2,2,1)
                imagesc(crop(skp));
                axis image
                colormap gray
                colorbar
                %caxis([prctile(skp(:),.001) prctile(skp(:),99)])
                drawnow
                
                
%                 % axis(ds/.1*[128.5166  356.6414  107.6861  300.1663]);
%                 colorbar
%                 colormap gray
%                 drawnow
%                 
%                 subplot(2,2,2)
%                 plot(f)
%                 title('cost')
%                 drawnow
%                 
%                 subplot(2,2,3)
%                 plot(dual_resid_s)
%                 hold on
%                 plot(primal_resid_s)
%                 hold off
%                 legend('s dual','s primal')
%                 title('residuals')
%                 drawnow
                
                 tic
            end
        end
        im_out(:,:,ccount) = sk;
    end
end

fprintf('total time: %.2f seconds \n',toc(bigT))
%%
gamma = .8;
%white_bal= [1 .7 1 ];
%white_bal = [1.1 .8 .5];
white_bal = [1.1 .9 .8];
im_out = im_out/prctile(im_out(:),99.9);
switch lower(color_proc)
    case('rgb')
        im_out_wb = cat(3,white_bal(1).*im_out(:,:,1),white_bal(2).*im_out(:,:,2),white_bal(3).*im_out(:,:,3));
    otherwise
        im_out_wb = im_out;
end
im_to_write = uint8(255*max(im_out_wb,0).^(gamma));
imagesc(im_to_write)

axis image
%axis([Nx/2 3*Nx/2 Ny/2 3*Ny/2])
im_stack{mag_count} = im_to_write;


