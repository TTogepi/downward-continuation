function [error] = rms(in,true)

[m,n] = size(in);
delta = true - in;
sum_error = sum(sum(delta .^ 2));
error = (sum_error / (m * n)) ^ (1/2);

end
