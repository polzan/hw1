function [a_tf, sigma2w_estim] = ar_model(x, N, autocorr_type, theor_noise_pow)
if nargin < 3
    autocorr_type = 'unbiased';    
end

if nargin < 4 && strcmp(autocorr_type, 'theoretical')
    theor_noise_pow = NaN;
end

switch autocorr_type
    case 'unbiased'
        r = autocorrelation_unbiased(x, N);
    case 'biased'
        r = autocorrelation_biased(x, N);
    case 'theoretical'
        r = theoretical_autocorr(theor_noise_pow, N);
    otherwise
        error('unknown autocorrelation type');
end
R = transpose(toeplitz(r(1:N)));
r_vec = r(2:N+1);
%a = -inv(R) * r_vec;
linsolve_opts.SYM = true;
a = linsolve(R, -r_vec, linsolve_opts);
a_tf = [1; a];

sigma2w_estim = r(1) + (r_vec'*a);
if ~isreal(sigma2w_estim) || sigma2w_estim < 0
    warning(['The AR model estimates a noise power of ' num2str(sigma2w_estim)]);
    sigma2w_estim = abs(real(sigma2w_estim));
end
end
