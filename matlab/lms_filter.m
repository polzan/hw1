function [c, e, y, k_c, k_ey] = lms_filter(d, x, N, mu)
K = length(d);
assert(length(x) >= K-1);
% Skip the first N-1 samples: we need a full vector x_k
c = zeros(N, K-N+1);
y = zeros(K-N, 1);
e = zeros(K-N, 1);
for k=1:K-N % -1 because we need the last d(k)
    x_k = flip(x(k : k + N - 1)); % The first should be the most recent
    y(k) = x_k.'*c(:,k);
    e(k) = d(k) - y(k); % should predict the next sample
    c(:,k+1) = c(:,k) + mu*e(k) .* conj(x_k);
end
k_c = N-1:K-1;
k_ey = N-1:K-2;
end
