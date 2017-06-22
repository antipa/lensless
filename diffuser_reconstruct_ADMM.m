%psf = imread('/Users/nick.antipa/Documents/Diffusers/Lensless/pco_2d_color/psf_med.tif');
%psfin = load('/Volumes/PhantomData/Diffusers''nstuff/colorPCO_collimated_calibration_far/pinhole_far_green_coherent.mat');

%psf = psfin.imave;
psf = double(imread('/Users/nick.antipa/Documents/Diffusers/Lensless/tapeycam/tape_40cm_psf.png'));
ds = 1;
if ds ~= 1
    psfd = imresize(double(psf(:,:,1)),ds,'box');
else
    psfd = psf;
end
h = psfd/norm(psfd(:),2);
Nx = size(h,2);
Ny = size(h,1);

%%
%imin = imread('cameraman.tif');
imin = imread('/Users/nick.antipa/Documents/Diffusers/Lensless/tapeycam/tape_40cm_shrey1.png');
if size(imin) ~= size(psfd)
    im = (double(imresize(imin,2*[Ny,Nx])));
else
    im = double(imin)/255;
end
h1 = figure(1),clf;

%define crop and pad operators to handle 2D fft convolution
pad = @(x)padarray(x,[size(h,1)/2,size(h,2)/2],0,'both');
cc = (size(h,2)/2+1):(3*size(h,2)/2);
rc = (size(h,1)/2+1):(3*size(h,1)/2);
crop = @(x)x(rc,cc);

H = fft2(ifftshift(pad(h)));
H_conj = conj(H);

% Define forward and adjoints that would be used in gradient descent
A = @(x)crop(ifftshift(real(ifft2(H.*fft2(ifftshift(x))))));
At = @(b)real(ifftshift(ifft2(H_conj.*fft2(ifftshift(pad(b))))));


%y = A(im);
%yin = imread('/Volumes/PhantomData/Diffusers''nstuff/Processing_results/20170419_152324_AVG_Video_rgb/AVG_Video.png');
%yin = demosaic(yin,'bggr');
yin = double(imread('/Users/nick.antipa/Documents/Diffusers/Lensless/tapeycam/tape_40cm_shrey1.png'));
if ds ~= 1
    y = imresize(double(yin(:,:,2)),ds,'box');
else
    y = double(yin)/255;
end

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
tau = .001;

% mu1 = 1e-1;
% mu2 = 1;
% mu3 = 4;
% tau = .1;
tauinc = 2;
taudec = tauinc;
muinc = 10;

I(1,1) = 1;
Smult = 1./(mu1*HtH + mu2*LtL+mu3);  %This is the frequency space division fo S update
if find(isinf(Smult))
    error('Dividing by zero!')
end
imagesc(abs(Smult))
colorbar
% Compute V denominator
CtC = pad(ones(size(im)));
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
%%
n = 0;
niter = 500;
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
        
        caxis([0 prctile(skp(:),99.9)])
       % axis(ds/.1*[128.5166  356.6414  107.6861  300.1663]);
        colorbar
        colormap gray
        drawnow
    end
    
end