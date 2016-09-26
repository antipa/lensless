% Script to trace one ray per diffuser pixel, then intersect each traced location
% with the pixel grid
filt = fspecial('gaussian',[200,200],3);
diffuser_surf = conv2(randn(200,200),filt,'same');
diffuser_surf = diffuser_surf/(max(diffuser_surf(:)));
imagesc(diffuser_surf);
px = 2.5;  %Micron

[Fx, Fy] = gradient(diffuser);   %2d gradient

%Convert gradient to correct units
Fx = Fx/px;
Fy = Fy/px;
%Trace




%Bin