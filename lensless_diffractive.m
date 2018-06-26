%%

h2 = figure(2);
clf


%%

xmin = -1.8;
xmax = 1.8;
npx = 1023.75;
x = linspace(xmin,xmax,npx); %mm
y = x;
dx = mean(diff(x))*1e-3;  %m
dx_mm = dx*1e3;
[X,Y] = meshgrid(x,y);

fno = 65;
gpitch = 5;
n_lenslets = 6;
D = gpitch*dx_mm*floor((xmax-xmin)/n_lenslets/(dx_mm*gpitch));
f = D*fno;
n= 1.56;
delta_n = 1-n;
R = f*delta_n;

sph = @(cx,cy,r)real(sqrt(r.^2-(X-cx).^2-(Y-cy).^2));

centers_y = (xmin+D/2):D:(xmax-D/2);
centers_x = centers_y;
diff_surf = zeros(size(X));
for n = 1:length(centers_x)
    for m = 1:length(centers_y)
         diff_surf = max(diff_surf,sph(centers_x(n),centers_y(m),R));
         %imagesc(diff_surf)
         %axis image
         %drawnow 
    end
end

zmin = 0;
zmax = f*1e-3;
nz = 10;
z_prop = linspace(zmin,zmax,nz);
lambda = 550e-9;
diffuser_phase = 2*pi*delta_n/lambda*diff_surf*1e-3;

fx = linspace(-1/2/dx,1/2/dx,numel(x));
fy = fx;
[Fx, Fy] = meshgrid(fx,fy);




    %gpitch = 8;
grating = zeros(size(X));
    for gp = 1:gpitch
        grating(gp:2*gpitch:end,:) = 1;
    end
    for gp = 1:gpitch
        grating(:,gp:2*gpitch:end) = 1;
    end

if gpitch == 0
    grating = ones(size(grating));
end
    ui = exp(-1i*diffuser_phase);
    for z0 = z_prop(end)
        
        up = propagate2(grating.*ui,lambda,z0,Fx,Fy);
        up_nograting = propagate2(ui,lambda,z0,Fx,Fy);
        psf = abs(up).^2;
        psf_nograting = abs(up_nograting).^2;
        set(0,'CurrentFigure',h2)
        imagesc(psf_nograting)
        axis image
        axis([408.2962  450.0717  243.5062  274.8523])
        
        drawnow
    end

psf_spect = ordfilt2(abs(fftshift(fft2(ifftshift(psf)))),n_lenslets^2,ones(n_lenslets,n_lenslets));
psf_nograting_spect = ordfilt2(abs(fftshift(fft2(ifftshift(psf_nograting)))),n_lenslets^2,ones(n_lenslets,n_lenslets));

%axis([  226.2963  288.4196  234.4934  281.1075])
figure(3),clf
imagesc(psf_spect(1:n_lenslets:end,1:n_lenslets:end))
axis image
%%

figure(1)
clf
plot(psf_nograting_spect(1:n_lenslets:end,middle(x))/max(psf_nograting_spect(:)))
hold on
plot(psf_spect(1:n_lenslets:end,middle(x))/max(psf_spect(:)))
legend('no grating','grating')






