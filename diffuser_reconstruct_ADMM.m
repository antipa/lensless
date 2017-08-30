%psf = imread('/Users/nick.antipa/Documents/Diffusers/Lensless/pco_2d_color/psf_med.tif');
%psfin = load('/Volumes/PhantomData/Diffusers''nstuff/colorPCO_collimated_calibration_far/pinhole_far_green_coherent.mat');

%psf = psfin.imave;
ds = 8;
demos = 0;
regularizer_type = 'tv';
erasure_type = 'none' %Choices: random, tile, none
%random_erasure = .05;   %If enabled, randomly delete pixels before recon. Use a number between 0 and 1 to do erasure, and 0 for full res
y_tiles = 4;
x_tiles = 4;
tile_gap = .5;  %Fraction of individual sensor width to leave as a gap
% ADMM stuff
tau = 1e-4;  %TV
tau2 = 8e-4;  %DCT
color_correct = 1;

mu1 = .2;
mu2 = .5;
mu3 = 1.6;
mu4 = 1.6; %DCT

mu_inc = 1.3;
mu_dec = 1.3;
resid_tol = 2;
if demos
    psf = demosaic(imread('/Users/nick.antipa/Documents/pinhole_far_green_00013.png'),'bggr');
    psf_in = psf(:,:,2)-100;
else
    psf = imread('/Users/nick.antipa/Documents/pinhole_far_green_00013.png');
    psf_in = load('/Users/nick.antipa/Documents/Diffusers/Lensless/pco_2d_color/pinhole_far_green_coherent.mat');
    psf = psf_in.imave;
    psf_in = double(psf(2:2:end,1:2:end))/2 + double(psf(1:2:end,2:2:end))/2-100;
    ds = ds/2;
end
%%
color_proc = 'rgb';

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


%mag_vec = 1.04;
mag_vec = 1.02;
mag_count = 0;

for mag_set = mag_vec
    mag_count = mag_count + 1;
    if mag_count == 1
        im_stack = cell(length(mag_vec),1);
    end
    ccount = 0;
    for color = color_vec
        object_close = 1;
        ccount = ccount + 1;

        if object_close
            %mag = .9744;
            mag = mag_set;
            %mag = 1.04;  %This is good for gracedog
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
            color_calibration = load('/Users/nick.antipa/Documents/Diffusers/Lensless/pco_2d_color/calibration/bayer_filter_coefficients.mat');
            %yinb = double(imread('/Users/nick.antipa/Documents/diffuser_cam_results/grace_dog_close_scale_102.png'));
            yinb = double(imread('/Users/nick.antipa/Documents/Diffusers/Lensless/pco_2d_color/lab_scene_ave-raw.png'));
            if numel(size(yinb)) == 2
                yinred = double(yinb(1:2:end,1:2:end))-100;
                yingreen = double(yinb(2:2:end,1:2:end))/2 + double(yinb(1:2:end,2:2:end))/2-100;
                yinblue = double(yinb(2:2:end,2:2:end))-100;
            else
                yinred = imresize(double(yinb(:,:,1)),.5,'box')-100;
                yingreen = imresize(double(yinb(:,:,2)),.5,'box')-100;
                yinblue = imresize(double(yinb(:,:,3)),.5,'box')-100;
            end
            if color_correct
                Ideconv = PCO_demosaic(yinb-100,color_calibration.couplingmatrix);
                %Ideconv = Ideconv/max(Ideconv(:)) * max(yinb(:) - 100);
                yinred = Ideconv(:,:,1);
                yingreen = Ideconv(:,:,2);
                yinblue = Ideconv(:,:,3);
            end
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
        % Precompute spectrum of L'L
        d = zeros(Ny,Nx);
        d(floor(Ny/2),floor(Nx/2)) = 1;
        % Uncropped forward
        Hfor = @(x)real(ifftshift(ifft2(H.*fft2(ifftshift(x)))));
        Hadj = @(x)real(ifftshift(ifft2(H_conj.*fft2(ifftshift(x)))));
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




        n = 0;
        niter = 70;
        
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
            Lskm1 = Lsk1;
            Lskm2 = Lsk2;
            [Lsk1,Lsk2] = L(sk);
            
            uk_1 = ukp_1;
            uk_2 = ukp_2;
            [ukp_1, ukp_2] = soft_2d_gradient(Lsk1+alpha2k_1/mu2, Lsk2+alpha2k_2/mu2, tau/mu2);
            vkp = Vmult.*(mu1*(alpha1k/mu1 + Hfor(sk)) + Cty);
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
            
            skp = real(ifftshift(ifft2(Smult .* fft2(ifftshift(skp_numerator)))));
            
            %Update dual and parameter for Hs=v constraint
            Hskp = Hfor(skp);
            r_sv = Hskp-vkp;
            alpha1k = alpha1k + mu1*r_sv;
            dual_resid_s(n) = mu1*norm(Hfor(sk-skp),'fro');
            primal_resid_s(n) = norm(r_sv);
            [mu1, mu1_update] = update_param(mu1,resid_tol,mu_inc,mu_dec,primal_resid_s(n),dual_resid_s(n));
            
            % Update dual and parameter for Ls=v
            [Lskp1, Lskp2] = L(skp);
            r_su_1 = Lskp1 - ukp_1;
            r_su_2 = Lskp2 - ukp_2;
            alpha2k_1 = alpha2k_1 + mu2*r_su_1;
            alpha2k_2 = alpha2k_2 + mu2*r_su_2;            
            dual_resid_u(n) = mu2*sqrt(norm(Lskm1 - Lsk1,'fro')^2+norm(Lskm2 - Lsk2,'fro')^2);
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
                fprintf('Mu updates: %i \t %i \t %i\n',mu1_update, mu2_update, mu3_update);
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
            fprintf('iter: %i \t cost: %.4f \t Primal v: %.4f \t Dual v: %.4f \t Primal u: %.4f \t Dual u: %.4f \t Primal w: %.4f \t Dual w: %.4f \t mu1: %.4f \t mu2: %.4f \t mu3: %.4f \n',...
        n,f(n),primal_resid_s(n), dual_resid_s(n),primal_resid_u(n), dual_resid_u(n),primal_resid_w(n), dual_resid_w(n),mu1,mu2,mu3)
            if mod(n,5)==0
                set(0,'CurrentFigure',h1);
                
                %fprintf('%.4f\t %.4f\t %.4f\t %.4f\t %.4f\t %.4f\n',ea1,ea2a,ea2b,ea3,e_primal,f)
                %fprintf('%.4f\n',norm(Hfor(skp) - vkp,'fro'))
                %nrm = TVnorm(sk) +
                subplot(2,2,1)
                imagesc(crop(skp));
                axis image

                caxis([prctile(skp(:),.001) prctile(skp(:),99.99)])
                % axis(ds/.1*[128.5166  356.6414  107.6861  300.1663]);
                colorbar
                colormap gray
                drawnow
                
                subplot(2,2,2)
                plot(f)
                title('cost')
                drawnow
                
                subplot(2,2,3)
                plot(dual_resid_s)
                hold on
                plot(primal_resid_s)
                hold off
                legend('s dual','s primal')
                title('residuals')
                drawnow
                
                
            end
        end
    im_out(:,:,ccount) = sk;
    end
end


    %%
gamma = .7;
white_bal= [1 1 1];
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

im_stack{mag_count} = im_to_write;

    
