function r = autocorrelation_unbiased(x, N)
K = length(x);
if nargin < 2
    N = K-1;
end

s = size(x);
if s(2) ~= 1
    x = transpose(x);
end
if s(2) ~= 1
    error('x is neither a row nor a column vector');
end

r = zeros(N+1, 1);
for n = 0:N
    x1 = x(n+1:K);
    x2 = x(1:K-n);
    r(n+1) = 1/(K-n) * (x2'*x1);
end
end
