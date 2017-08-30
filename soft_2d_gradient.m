function [vsoft, hsoft] =  soft_2d_gradient(v,h,tau)

mag = sqrt(cat(1,v,zeros(1,size(v,2))).^2 + ...
    cat(2,h,zeros(size(h,1),1)).^2);
magt = soft(mag,tau);
mmult = magt./mag;
mmult(mag==0) = 0;
vsoft = v.*mmult(1:end-1,:);
hsoft = h.*mmult(:,1:end-1);