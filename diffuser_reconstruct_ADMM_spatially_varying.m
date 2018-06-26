% Load data

psf_in = load('/Users/nick.antipa/Documents/Diffusers/Lensless/pco_2d_color/calibration/low_rank/PCA_v1_10pc.mat');
meas_in = load('/Users/nick.antipa/Documents/Diffusers/Lensless/pco_2d_color/calibration/low_rank/rawData_flower_PCA100_v1.mat');
%%
autotune = 0;
save_metrics = 1;
ds = 4;
N_comps = 1;
pc = double(imresize(psf_in.pc(:,:,1:N_comps),1/ds,'box'));
weights = double(imresize(psf_in.weights(:,:,1:N_comps),1/ds,'box'));
 %pc_n = norm(pc(:,:,1),'fro');
 %weights = weights/pc_n;
 %pc = pc/pc_n;
[Ny, Nx] = size(pc(:,:,1));
Nx = Nx/2;
Ny = Ny/2;
s_orig = imresize(meas_in.gndtruth,1/ds,'box');
%b = imresize(meas_in.btruth,1/ds,'box');


mu_inc = 1.3;
mu_dec = 1.3;
resid_tol = 4;
% Tuning params
sc = 1;
mu1 = .004;
mu2 = 2.7;
mu3 = 6.5;
mu4 = .5;
mu5 = .8;
mu6 = 12.6;
tau = 5e1;



% Setup diagonal operators for x update
Lambda = @(x)x.*weights;   %x must be 3D stack of ims
LambdaT = @(x)sum(x.*weights,3);
lamb_weight = sum(weights.^2,3);
Lambda_mult = 1./(mu2 + mu1*lamb_weight);
%%
% PSF-spectrum
H_freq = fft2(pc);
H_conj_freq = conj(H_freq);
HtH_freq = 1./(mu4*abs(H_freq.*H_conj_freq) + mu1);

H = @(x)real(ifftshift2(ifft2(H_freq.*fft2(x))));
H_adj = @(x)real(ifftshift2(ifft2(H_conj_freq.*fft2(x))));
%Crop 
pad = @(x)padarray(x,[Ny/2,Nx/2],0,'both');
cc = (Nx/2+1):(3*Nx/2);
rc = (Ny/2+1):(3*Ny/2);
crop = @(x)x(rc,cc);
CtC = pad(ones(Ny,Nx));
CtC_mult = 1./(mu5+CtC);
b = crop(sum(H(Lambda(s_orig)),3));  %Create fake data
Ctb = pad(b);

Ibart = @(x)repmat(x,[1 1 N_comps]);
Ibar = @(x)sum(x,3);
% Setup TV operators for s and u updates
Ltv = @(P1,P2)cat(1,P1(1,:),diff(P1,1,1),-P1(end,:)) + cat(2,P2(:,1),diff(P2,1,2),-P2(:,end));

% Forward finite difference (returns two images!) [P1, P2] = L(x)
L = @(D)deal(-diff(D,1,1),-diff(D,1,2));

% Get spectrum of laplacian
lapl = zeros(2*Ny,2*Nx);
lapl(1) = 4;
lapl(1,2) = -1;
lapl(2,1) = -1;
lapl(1,end) = -1;
lapl(end,1) = -1;
LtL = abs(fft2(lapl));

L_mult_freq = 1./(mu2 + mu6 + mu3*LtL);

% Initialize variables

s = zeros(2*Ny,2*Nx);
x = s;
v = zeros(size(weights));
w = v;
y = s;
p = s;

alpha1 = w;
alpha2 = s;
alpha3_v = zeros(2*Ny-1,2*Nx);
alpha3_h = zeros(2*Ny,2*Nx-1);
alpha4 = v;
alpha5 = y;
alpha6 = s;
niter = 1000;
s = zeros(2*Ny,2*Nx);
for k = 1:niter
    % Primals
    
    
    [uv, uh] = L(s);
    [uv, uh] = soft_2d_gradient(uv + alpha3_v/mu3, uh + alpha3_h, tau/mu3);
    
    %y
    y_km1 = y;
    y = CtC_mult .* (mu5*sum(v,3) + alpha5 + Ctb);
    
    %x
    x = Lambda_mult .* (LambdaT(mu1*w - alpha1) + mu2*s - alpha2);
    
    %v
    v_km1 = v;
    v = ItI_inv(mu4*H(w)+alpha4 + Ibart(mu5*y - alpha5),N_comps,mu4/mu5)/mu5;
    
    %w
    wnum = (mu1*Lambda(x) + alpha1 + H_adj(mu4*v - alpha4));
    w_km1 = w;
   %HtH_inv = real(ifftshift2(ifft2(HtH_freq.*fft2(ifftshift2(HtH_in)))));
    w = real(ifftshift2(ifft2(HtH_freq.*fft2(ifftshift2(wnum)))));
    
    % p
    p_km1 = p;
    p = max(s+alpha6/mu6,0);
    
    %s
    snum = mu2*x + alpha2 + Ltv(mu3*uv-alpha3_v,mu3*uh-alpha3_h) + mu6*p-alpha6;
    s_km1 = s;
    s = real(ifftshift(ifft2(L_mult_freq.*fft2(ifftshift(snum)))));
    
    % Duals
    %a1
    r_xw = Lambda(x)-w;
    alpha1 = alpha1 + mu1*r_xw;
    %a2
    r_xs = x-s;
    alpha2 = alpha2 + mu2*r_xs;
    %a3

    [Lsv, Lsh] = L(s);
    r_su_v = Lsv - uv;
    r_su_h = Lsh - uh;
    alpha3_v = alpha3_v + mu3*(r_su_v);% + beta_kp1*(alpha2k_1p - alpha2k_1);
    alpha3_h = alpha3_h + mu3*(r_su_h);% + beta_kp1*(alpha2k_2p - alpha2k_2);
    %a4
    r_wv = H(w)-v;
    alpha4 = alpha4 + mu4*r_wv;
    %a5
    r_vy = sum(v,3)-y;
    alpha5 = alpha5 + mu5*(r_vy);
    %a6
    r_sp = s-p;
    alpha6 = alpha6 + mu6*r_sp;
    grad_mag = sqrt(cat(1,Lsv,zeros(1,size(Lsv,2))).^2+cat(2,Lsh,zeros(size(Lsh,1),1)).^2);
    subplot(1,2,1)
    imagesc(s)
    colorbar
    axis image
    subplot(1,2,2)
    imagesc(uv)
    axis image
    
    drawnow
    if mod(k-1,10)==0
        fprintf(['Total cost\t Data fidelity   ||Ls||_1',...
        '    norm(Lambda(x)-w)   Lambda(X)-w dual    norm(L(s)-u)   norm(x-s)',...
        'norm(H(w)-v)   norm(I_bar(v)-y)  sum(negs)',...
        '\n'])
    end
    %npix = 4*Ny*Nx;
    data_fidelity = norm(crop(sum(H(Lambda(s)),3))-b)^2;
    TV_penalty = sum(grad_mag(:));
    
    x_resid = norm(vec(r_xw));
    x_resid_dual = norm(mu1*LambdaT(w-w_km1),'fro');
    
    u_resid = sqrt(norm(r_su_v,'fro')^2+norm(r_su_h,'fro')^2);
    [u_dual_v, u_dual_h] = L(s-s_km1);
    u_resid_dual = mu3*sqrt(norm(u_dual_v,'fro')^2 + norm(u_dual_h,'fro'));
    
    s_resid = norm(r_xs,'fro');
    s_resid_dual = mu2*norm(s-s_km1,'fro');
    
    w_resid = norm(vec(r_wv));
    w_resid_dual = mu4*norm(vec(H(w - w_km1)));
    
    v_resid = norm(sum(v,3)-y);
    v_resid_dual = mu4*norm(vec(Ibart(y - y_km1)));
    
    y_resid = norm(r_vy,'fro');
    y_resid_dual = mu5*norm(sum(v - v_km1,3));
    
    p_resid = norm(r_sp,'fro');
    p_resid_dual = norm(mu6*(p - p_km1));
    
    if autotune
        [mu1, mu1_update] = update_param(mu1,resid_tol,mu_inc,mu_dec,x_resid,x_resid_dual);
        [mu2, mu2_update] = update_param(mu2,resid_tol,mu_inc,mu_dec,s_resid,s_resid_dual);
        [mu3, mu3_update] = update_param(mu3,resid_tol,mu_inc,mu_dec,u_resid,u_resid_dual);
        [mu4, mu4_update] = update_param(mu4,resid_tol,mu_inc,mu_dec,w_resid,w_resid_dual);
        [mu5, mu5_update] = update_param(mu5,resid_tol,mu_inc,mu_dec,y_resid,y_resid_dual);
        [mu6, mu6_update] = update_param(mu6,resid_tol,mu_inc,mu_dec,p_resid,p_resid_dual);
        %If any mu updated, update cached factorizations.
        
        if mu1_update || mu2_update || mu3_update || mu4_update || mu5_update || mu6_update
            Lambda_mult = 1./(mu2 + mu1*lamb_weight);
            HtH_freq = 1./(mu4*H_freq.*H_conj_freq + mu1);
            L_mult_freq = 1./(mu2 + mu6 + mu3*LtL);
            CtC_mult = 1./(mu5+CtC);
            fprintf('Updating parameters:\n mu1\t mu2 \t mu3 \t mu4 \t mu5 \t mu6\n')
            fprintf([repmat('%i \t',[1,6]),'\n'],mu1_update, mu2_update, mu3_update,...
                mu4_update, mu5_update, mu6_update);
        end

    end
    


    fprintf('%.4e\t',data_fidelity + tau*TV_penalty);
    fprintf('%.4e\t',data_fidelity);
    fprintf('%.4e\t',TV_penalty);
    fprintf('%.4e\t',x_resid);
    fprintf('%.4e\t',x_resid_dual);
    fprintf('%.4e\t',u_resid);
    fprintf('%.4e\t',s_resid);
    fprintf('%.4e\t',w_resid);
    fprintf('%.4e\t',v_resid);
    fprintf('%.4e\t',y_resid);
    fprintf('%.4e\n',sum(vec(s.*double(s<0))))
    if save_metrics
        if k == 1
            date_of_run = datestr(now,'yyyymmdd_HHMMSS');
            fname = ['/Users/nick.antipa/Documents/developer/logs/admm_psf_decomp',date_of_run,'_log.txt'];
            metric_file = fopen(fname,'w');
            fprintf(metric_file,['Total_cost\tdata_fidelity\tTV_norm\t ',...
                '||Lambda(x)-w||\t||mu1*LambdaT(w-w_km1)||\t',...
                '||Ls-u||\t||mu3*L(s-s_km1)||\t',...
                '||s-x||\t||mu2*(s-s_km1)||\t',...
                '||H(w)-v||\t||mu5*H(w-w_km1)||\t',...
                '||I_bar(v)-y||\t||mu4*Ibart(y-y_km1)||\t',...
                '||mu5*Ibar(v-v_km1)||\t',...
                'sum_negs\tnumber_negs\n',...
                ]);
        else
            fopen(fname,'a+');
        end
       fprintf(metric_file,[repmat('%f\t',[1,15]),'%f\n'],data_fidelity+tau*TV_penalty,...
           data_fidelity, TV_penalty,...
           x_resid, x_resid_dual, u_resid, u_resid_dual, s_resid, s_resid_dual,...
           w_resid, w_resid_dual, v_resid, v_resid_dual, y_resid_dual,...
           sum(vec(s.*double(s<0))), abs(sum(sum(s<0))));
       fclose(metric_file);
    end
    %fprintf('norm(Lambda(x)-w): %.2f',norm(Lambda(x)(:)-w(:),'fro'));
    
end





