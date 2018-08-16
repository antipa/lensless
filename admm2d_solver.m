function [xhat, f] = admm2d_solver(y,h,sett_file)
% Solves xhat using ADMM from measurement b
% Inputs: y - array measurement
% x_init: initialization for variable. If empty, use zeros.
% h: psf (as double);
% sett_file : m-file that generates all settings. Use admm2d_settings.m as
% default (or template to make your own)
h1 = figure(1);
clf
tic

h = h/norm(h,'fro');
if isempty(sett_file)
    admm2d_settings
else
    run(sett_file);
end
tau = cast(tau,'like',y);
[Ny, Nx] = size(h);

pad = @(x)padarray(x,[floor(Ny/2), floor(Nx/2)],'both');
%cc = (Nx/2+1):(3*Nx/2);
%rc = (Ny/2+1):(3*Ny/2);
crop = @(x)x(floor(Ny/2)+1:end-floor(Ny/2),floor(Nx/2)+1:end-floor(Nx/2));

H = fft2(ifftshift(pad(h)));
H_conj = conj(H);

% Adjoint of TV finite difference
% Precompute spectrum of L'L

% Uncropped forward
Hfor = @(x)real(ifft2(H.*fft2(x)));
Hadj = @(x)real(ifft2(H_conj.*fft2(x)));
% Initialize variables
Cty = pad(y);
sk = 0*Cty;
alpha1k = sk;
alpha3k = sk;
alpha3kp = sk;

tv_iso = 1;

% TV stuff
Ltv = @(P1,P2)cat(1,P1(1,:),diff(P1,1,1),-P1(end,:)) + cat(2,P2(:,1),diff(P2,1,2),-P2(:,end));

% Forward finite difference (returns two images!) [P1, P2] = L(x)
L = @(D)deal(-diff(D,1,1),-diff(D,1,2));

lapl = sk;
lapl(1) = 4;
lapl(1,2) = -1;
lapl(2,1) = -1;
lapl(1,end) = -1;
lapl(end,1) = -1;
LtL = abs(fft2(lapl));
%alpha2k_1 = zeros(2*Ny-1,2*Nx);
%alpha2k_2 = zeros(2*Ny,2*Nx-1);
alpha2k_1 = sk(1:end-1,:);
alpha2k_2 = sk(:,1:end-1);



% DCT stuff
DCT_adj = @(Q)idct2(Q);
DCT_for = @(D)dct2(D);

%LtL = ones(2*Ny,2*Nx);%abs(fft2(dct2(pad(d))));
alpha4k = sk;
alpha4kp = sk;

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
if find(isinf(gather(Smult)))
    error('Dividing by zero!')
end


% Compute V denominator
CtC = pad(ones(Ny, Nx,'like',y));
Vmult = 1./(CtC + mu1);
Hskp = zeros(2*Ny, 2*Nx, 'like',Vmult);



n = 0;

f.dual_resid_s = zeros(1,niter,'like',y)./0;
f.primal_resid_s = f.dual_resid_s;
f.dual_resid_u = f.dual_resid_s;
f.primal_resid_u = f.dual_resid_u;
f.dual_resid_w = f.dual_resid_s;
f.primal_resid_w = f.dual_resid_s;

f.cost = f.dual_resid_s;
[ukp_1, ukp_2] = L(zeros(2*Ny, 2*Nx,'like',y));
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
    f.dual_resid_s(n) = mu1*norm(Hsk - Hskp,'fro');
    f.primal_resid_s(n) = norm(r_sv,'fro');
    [mu1, mu1_update] = update_param(mu1,resid_tol,mu_inc,mu_dec,f.primal_resid_s(n),f.dual_resid_s(n));
    alpha1k = alpha1k + mu1*r_sv;
    %[mu1, mu1_update] = update_param(mu1,resid_tol,mu_inc,mu_dec,f.primal_resid_s(n),f.dual_resid_s(n));
    
    % Update dual and parameter for Ls=v
    [Lskp1, Lskp2] = L(skp);
    r_su_1 = Lskp1 - ukp_1;
    r_su_2 = Lskp2 - ukp_2;
    f.dual_resid_u(n) = mu2*sqrt(norm(Lskm1 - Lsk1,'fro')^2+norm(Lskm2 - Lsk2,'fro')^2);
    f.primal_resid_u(n) = sqrt(norm(r_su_1,'fro')^2 + norm(r_su_2,'fro')^2);
    [mu2, mu2_update] = update_param(mu2,resid_tol,mu_inc,mu_dec,f.primal_resid_u(n),f.dual_resid_u(n));
    alpha2k_1 = alpha2k_1 + mu2*r_su_1;
    alpha2k_2 = alpha2k_2 + mu2*r_su_2;
    %[mu2, mu2_update] = update_param(mu2,resid_tol,mu_inc,mu_dec,f.primal_resid_u(n),f.dual_resid_u(n));
    
    % Update nonnegativity dual and parameter (s=w)
    r_sw = skp-wkp;
    f.dual_resid_w(n) = mu3*norm(sk - skp,'fro');
    f.primal_resid_w(n) = norm(r_sw,'fro');
    [mu3, mu3_update] = update_param(mu3,resid_tol,mu_inc,mu_dec,f.primal_resid_w(n),f.dual_resid_w(n));
    alpha3k = alpha3k + mu3*r_sw;
    %[mu3, mu3_update] = update_param(mu3,resid_tol,mu_inc,mu_dec,f.primal_resid_w(n),f.dual_resid_w(n));
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
    f.cost(n) = norm(crop(Hskp)-y,'fro')^2 + tau*TVnorm(skp);
    if mod(n,disp_interval)==0
        toc
        fprintf('iter: %i \t cost: %.4f \t Primal v: %.4f \t Dual v: %.4f \t Primal u: %.4f \t Dual u: %.4f \t Primal w: %.4f \t Dual w: %.4f \t mu1: %.4f \t mu2: %.4f \t mu3: %.4f \n',...
            n,f.cost(n),f.primal_resid_s(n), f.dual_resid_s(n),f.primal_resid_u(n), f.dual_resid_u(n),f.primal_resid_w(n), f.dual_resid_w(n),mu1,mu2,mu3)
        
        set(0,'CurrentFigure',h1);
        
        %fprintf('%.4f\t %.4f\t %.4f\t %.4f\t %.4f\t %.4f\n',ea1,ea2a,ea2b,ea3,e_primal,f)
        %fprintf('%.4f\n',norm(Hfor(skp) - vkp,'fro'))
        %nrm = TVnorm(sk) +
        %subplot(2,2,1)
        imagesc((sk));
        axis image
        colormap gray
        drawnow
        tic
    end
end

xhat = skp;


