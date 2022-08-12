function [grad,fit,e_rms]=Guo_fs_rmse(f,f_check,h,Cell,iter,omiga,flag)
 
[m,n] = size(f);
M = m * n;
kernel = zeros(M,M);
b = zeros(M,1);
delta_z = h;
dx = Cell;
dy = Cell;

table_wait1 = waitbar(0,'The program is running');

for i = 1 : m
for j = 1 : n
    
    b((i-1) * n + j) = f(i,j);
    
end
end

parfor i = 1 : M
for j = 1 : M
    
    tab_yi = mod(j,n);
    tab_xi = 1 + floor(j / n);    
    tab_yii = mod(i,n);
    tab_xii = 1 + floor(i / n);
    
    xi =  (tab_xi - 1) * dx;
    yj =  (tab_yi - 1) * dy;
    xii = (tab_xii - 1) * dx;
    yjj = (tab_yii - 1) * dy;
    
    kernel(i,j) = new_kernel(xi,yj,delta_z,xii,yjj,dx,dy) / (pi * 2);
    
end
end

%%µü´úÇó½â
% %*********************************************************
% iter = 140;
% %**********************************************************
x = b; 
grad = zeros(1,iter);
fit = zeros(1,iter);
e_rms = zeros(1,iter);

for i = 1 : iter
    
    temp =  kernel * x;
    R = b - temp;
    x = x + omiga * kernel * R;
  
    temp2 = kernel * x;  %up_continuation
    
    temp_U1 = zeros(m,n);
    temp_U2 = zeros(m,n);
    for ii = 1 : m
    for jj = 1 : n 
        
    temp_U1(ii,jj) = x((ii-1) * n + jj);
    temp_U2(ii,jj) = temp2((ii-1) * n + jj);
    
    end
    end
    

    grad(i) =  grad_U_l2(temp_U1,dx);
    fit(i) = sum(sum((temp2 - b).^2)) / sum(sum(b.^2));
    fit(i) = fit(i) ^ (1/2);
    
    if flag == 1
    e_rms(i) = rms(temp_U1,f_check);
    end
    
end

close(table_wait1);
table_wait2 = waitbar(1,'The program completed');

end