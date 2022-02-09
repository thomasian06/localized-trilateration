function [bias, r, rel_sig] = get_noisy_measurements(...
    t, trajectory, pseudolites, cpoly, sig_r0, dfac)
%% Function get_noisy_measurements
%
% Purpose:
% from a trajectory and pseudolite locations, generate noisy measurements
%
% Inputs:
% t - time vector (row vector)
% trajectory - traj (n rows, T columns)
% pseudolites - list of positions of pseudolites
% cpoly - clock polynomial coefficients ()
% sig_r0 - standard range deviation of highest power signal (20 m)
% dfac - coefficient of noise for rel_power (0.05)
%
% 

    [npseudolites, dims] =  size(pseudolites);
    T = length(t);
    r = reshape(sqrt(sum((trajectory - ...
        reshape(pseudolites', [dims, 1, npseudolites])).^2, 1)),...
        [T, npseudolites])';
    
    rel_sig = (r ./ min(r)); % figure out the SNR ratios
    r = r + sig_r0*rel_sig.*randn(size(r)); % apply adjusted noise
%     figure(10)
%     plot(t', rel_power')
    rel_sig = rel_sig + dfac*rel_sig.*randn(size(rel_sig));
    
    bias = polyval(cpoly, t);
    r = r + bias*3e8;
    
end