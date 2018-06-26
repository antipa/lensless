N = 809;
A = generate_2D_MURA(N);
up = 3;
pxdown = 3;

zvec = linspace(0,200e-6,100);
pad = @(x)padarray(x,[ceil(N/4) ceil(N/4)],'both');
Aup = imresize(pad(A),up,'nearest');
crop = @(x)x(ceil(N/2)+1:end-ceil(N/2),ceil(N/2)+1:end-ceil(N/2));
fftconv = @(x)(ifftshift(ifft2(fft2(abs(x)).^2)));
ps = 3e-6;
fx = linspace(-1/2/ps*up,1/2/ps*up,size((Aup),1));
[Fx, Fy] = meshgrid(fx,fx);
lambda= 550e-9;
dint_lam = zeros(size(Fx));
lamvec = linspace(450e-9,500e-9,10);
%lamvec = 510e-9;
for z0 = zvec(end)
    for lambda = lamvec
        d = propagate2(Aup,lambda,z0,Fx,Fy);
        dint_lam = dint_lam + abs(d).^2/length(lamvec);
        imagesc(dint_lam)
        title(sprintf('lambda = %.1f nm',lambda*1e9));
        drawnow
    end
    %d = abs(propagate_plane_decomp((1:N*up)*(3e-6)/up,(1:N*up)*(3e-6)/up,Aup,z0,550e-9,1)).^2;

    dint = imresize(abs(d).^2,1/pxdown,'box');
    dint_lam_ds = imresize(dint_lam,1/pxdown,'box');
    dac = fftconv(dint);
    dac_lam = fftconv(dint_lam_ds);
    [m,i] = max(dac);
    x = (1:size(dac,2))*ps*pxdown/up;
    x = x-mean(x);
    dac_crop = dac(i(m==max(m)),:);
    dac_crop_lam = dac_lam(i(m==max(m)),:);
    plot(x*1e6,dac_crop/max(dac_crop));
    hold on
    plot(x*1e6,dac_crop_lam/max(dac_crop_lam))
    %imagesc(dint)
    xlabel('\mu m')
    legend('coherent','incoherent')
    hold off
    

    drawnow
    
end

