function [c_opt, a] = predictor(x, N, theor_noise_pow)
a = ar_model(x, N, 'theoretical',theor_noise_pow);   %use the theoretical one
c_opt = -a(2:N+1);
end