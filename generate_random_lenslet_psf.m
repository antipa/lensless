function psf = generate_random_lenslet_array(x0,y0,Xs, Ys, R_lenslet, varargin)
% Generates a lenslet array with curvature R_lenslet, with each lenslet
% centered at xi, yi in physical units. Xs and Ys define the grid
% and must be in physical units matching R.
% Optional 6th argument
% adds curvature (i.e. a main lens).


sph = @(x0,y0,R)real(sqrt(R^2-(Xs-x0).^2 - (Ys-y0).^2));

for n = 