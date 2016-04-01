function [c_opt, f_n, a] = predictor(x, N)
a = ar_model(x, N);
c_opt = -a(2:N+1);
f_n = conv(x,a);
end