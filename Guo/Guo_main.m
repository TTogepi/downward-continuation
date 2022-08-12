function [expectation]=Guo_main(f,h,Cell,iter,omiga)
 
[m,n] = size(f);
expectation = zeros(m,n);
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

% %**********************************************************
% iter = 140;
% %**********************************************************

x = b; 

for i = 1 : iter
    
    temp =  kernel * x;
    R = b - temp;
    x = x + omiga * kernel * R;
    
end

for i = 1 : m
for j = 1 : n 
    
    expectation(i,j) = x((i-1) * n + j);
    
end
end

close(table_wait1);
table_wait2 = waitbar(1,'The program completed');

end