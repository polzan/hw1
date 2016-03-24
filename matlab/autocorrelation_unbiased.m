function r = autocorrelation_unbiased(x)
K = length(x);
r = zeros(K, 1);
s = size(x);
if s(2) ~= 1
    x = transpose(x);
end
if s(2) ~= 1
    error('x is neither a row nor a column vector');
end
for n = 0:K-1
    x1 = x(n+1:K);
    x2 = x(1:K-n);
    r(n+1) = 1/(K-n) * (x2'*x1);
end
end
