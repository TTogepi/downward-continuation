function [xl2,fit,e_rms] =  intergral_iter_L_rmse(f, f_check,Cell, h, iter,flag)

[m,n] = size(f);
xl2 = zeros(1,iter);
fit = zeros(1,iter);
e_rms = zeros(1,iter);

M = m * n;
kernel = zeros(M,M);
ub = zeros(M,1);
ua = zeros(M,1);

dx = Cell;
dy = Cell;

for i = 1 : m
for j = 1 : n
    
    ub((i-1) * n + j) = f(i,j);
    
end
end

table_wait1 = waitbar(0,'The program is running');

parfor i = 1 : M
for j = 1 : M
    
    tab_yi = 1 + mod(j-1,n);
    tab_xi = 1 + floor((j-1) / n);    
    tab_yii = 1 + mod(i-1,n);
    tab_xii = 1 + floor((i-1) / n);
    
    xi =  (tab_xi - 1) * dx;
    yj =  (tab_yi - 1) * dy;
    xii = (tab_xii - 1) * dx;
    yjj = (tab_yii - 1) * dy;
    
    cosrn = h / ((xi -xii) ^ 2 + (yj - yjj) ^ 2 + h^2) ^ (1/2) ;
    kernel(i,j) = cosrn * dx * dy / (2 * pi * ((xi -xii) ^ 2 + (yj - yjj) ^ 2 + h^2));
   
end
end

for i = 1 : iter
 
     ua = ua + ub - kernel * ua;
     
     
     temp2 = kernel * ua;  %up_continuation

     temp_U1 = zeros(m,n);
     temp_U2 = zeros(m,n);
     for ii = 1 : m
     for jj = 1 : n 

     temp_U1(ii,jj) = ua((ii-1) * n + jj);
     temp_U2(ii,jj) = temp2((ii-1) * n + jj);
     
     end
     end
    
    xl2(i) = sum(sum(ua.^2));
    xl2(i) = xl2(i) ^ (1/2);
    fit(i) = sum(sum((temp2 - ub).^2)) / sum(sum(ub.^2));
    fit(i) = fit(i) ^ (1/2);
    
    if flag == 1
    e_rms(i) = rms(temp_U1,f_check);
    end
   
end

close(table_wait1);
table_wait2 = waitbar(1,'The program completed');

end
