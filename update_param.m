function [mu_out, mu_update] = update_param(mu,resid_tol,mu_inc,mu_dec,r,s)
    if r > resid_tol*s
        mu_out = mu*mu_inc;
        mu_update = 1;
    elseif r*resid_tol < s
        mu_out = mu/mu_dec;
        mu_update = -1;
    else
        mu_out = mu;
        mu_update = 0;
    end
end