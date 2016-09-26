%%
ny = 256;
nx = ny;
lambda = .55;
crop = @(x)x(501:500+ny,501:500+nx);
%in = load('/Users/nick.antipa/Documents/MATLAB/Output/1deg_flea_4f_gptie_dz381um_surface_1552.mat');
%diffuser_in = crop(in.filtered)*1e-6;   %Scale height
%x = in.x*1e-6;   %get physical units grid from file
filt = fspecial('gaussian',[150,150],10);
diffuser_in = conv2(rand(ny,nx),filt,'same');
diffuser_in = diffuser_in./max(diffuser_in(:));
x = (1:ny)*.25;
px = mean(diff(x));
W = -.68*diffuser_in;
[Fx,Fy] = gradient(W);
Fx = Fx/px;
Fy = Fy/px;
[Fxx, Fxy] = gradient(Fx);
Fxx = Fxx/px;
Fxy = Fxy/px;
[Fyx,Fyy] = gradient(Fy);
Fyy = Fyy/px;
Fyx = Fyx/px;

zv = 0:10:7000;
for z = zv;
    Iz = ones(size(diffuser_in)).*(1-(Fxx+Fyy)*z);
    g = (1:size(diffuser_in))*px;
    [X,Y] = meshgrid(g,g);

    Gx = Fx*z;
    Gy = Fy*z;
    Iz_interp = max(interp2(X,Y,Iz,X-Gx,Y-Gy),0);
    imagesc(Iz_interp);

    %colorbar
    drawnow
    
    
end