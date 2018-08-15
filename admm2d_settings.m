ds = 4;  %Downsampling
demos = 1;
niter = 800;
regularizer_type = 'tv';
camera_type = 'ids_bsi';
erasure_type = 'none' %Choices: random, tile, none
object_close = 1;   %Magnify PSF to account for object distance changes
color_correct = 1;
psf_color_proc = 'mono';
mag_vec = 1.0;  %Adjust PSF distance
color_proc = 'rgb';
random_erasure = .25;   %If enabled, randomly delete pixels before recon. Use a number between 0 and 1 to do erasure, and 0 for full res
y_tiles = 4;
x_tiles = 4;
tile_gap = .5;  %Fraction of individual sensor width to leave as a gap
% ADMM stuff
tau_vec = [8e-4 4e-4 10e-4];  %TV
tau = 50;
tau2 = .5e-7;  %DCT
disp_interval = 10;

mu1 = .14;
mu2 = 3;
mu3 = 3.2;
mu4 = .7; %DCT

mu_inc = 1.2;
mu_dec = 1.2;
resid_tol = 1.5;