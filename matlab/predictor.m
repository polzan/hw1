function [c_opt, f_n, a] = predictor(x, N)
K = length(x);
r_whole=autocorrelation_unbiased(x);
r_n = r_whole(2:N+1);
R_n = toeplitz(r_whole(1:N)); 
c_opt = inv(R_n)*r_n;
a = [1;-c_opt];
f_n = conv(x,a);
end