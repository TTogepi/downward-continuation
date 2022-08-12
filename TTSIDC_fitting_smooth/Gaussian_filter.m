function [expectation, Gxx] = Gaussian_filter(f,flag,cell,sigma,k)

%% Laplacian of gaussain operator

[m,n] = size(f);
expectation = zeros(m,n);
U = wextend(2,'symh',f,[k,k]);
Gxx = zeros(2 * k + 1,2 * k + 1);
Gyy = zeros(2 * k + 1,2 * k + 1);

for i = 1 : 2 * k + 1
for j = 1 : 2 * k + 1

x = i - k - 1;
% x = x * cell;
y = j - k - 1;
% y = y * cell;

Gxx(i,j) = (x^2 - sigma^2) / sigma^4 * exp(- (x^2 + y^2) / (2 * sigma^2));
Gyy(i,j) = (y^2 - sigma^2) / sigma^4 * exp(- (x^2 + y^2) / (2 * sigma^2));

end
end

miu_xx = sum(sum(Gxx))/ (2 * k + 1) ^ (2);
miu_yy = sum(sum(Gyy))/ (2 * k + 1) ^ (2);

for i = 1 : 2 * k + 1
for j = 1 : 2 * k + 1
    
    Gxx(i,j) = (Gxx(i,j) - miu_xx) / cell^2;
    Gyy(i,j) = (Gyy(i,j) - miu_yy) / cell^2;
    
end
end

if flag == 1
    kernel = Gxx;    
else    
    kernel = Gyy;    
end

for i = 1 : m 
for j = 1 : n
    
    for ii = 1 : 2 * k + 1
    for jj = 1 : 2 * k + 1
        
        expectation(i,j) = expectation(i,j) + U(i+k+(ii-k-1),j+k+(jj-k-1)) * kernel(ii,jj);
        
    end     
    end
    
end
end

end