function [expectation] =  intergral_iter(f, Cell, h, iter)

[m,n] = size(f);
expectation = zeros(m,n);
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
    
    tab_yi = mod(j,n);
    tab_xi = 1 + floor(j / n);    
    tab_yii = mod(i,n);
    tab_xii = 1 + floor(i / n);
     
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
     
end

for i = 1 : m
for j = 1 : n 
    
    expectation(i,j) = ua((i-1) * n + j);
    
end
end

close(table_wait1);
table_wait2 = waitbar(1,'The program completed');

end