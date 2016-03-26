function [f_n,c_opt] = predictor(x, N)
K = length(x);
r_whole=autocorrelation_unbiased(x);
r_n = r_whole(2:N+1);
R_n = toeplitz(r_whole(1:N));
f_n = x;
c_opt = inv(R_n)*r_n;
for i = N+1:K-1
f_n(i) = x(i+1) - fliplr(c_opt.')*x(i-N:i-1);
end
end