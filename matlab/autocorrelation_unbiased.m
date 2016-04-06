function r = autocorrelation_unbiased(x, N, whole)
K = length(x);
if nargin < 2
    N = K-1;
end
if nargin < 3
    whole = false;
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

if whole
    r = [flip(conj(r(2:N+1))); r];
end
end
