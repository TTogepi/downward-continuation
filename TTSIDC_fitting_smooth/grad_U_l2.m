function [out,Ux,Uy] = grad_U_l2(in,Cell)

[m,n] = size(in);
sum_Ux = 0;
sum_Uy = 0;
Ux  = zeros(m-1,n-1);
Uy  = zeros(m-1,n-1);

for i = 2 : m 
for j = 2 : n 
    
    Ux(i-1,j-1) = (in(i,j) - in(i-1,j)) / Cell;
    Uy(i-1,j-1) = (in(i,j) - in(i,j-1)) / Cell;
    
    sum_Ux = sum_Ux +Ux(i-1,j-1)^2;
    sum_Uy = sum_Uy +Uy(i-1,j-1)^2;
    
end
end

out =  (sum_Ux + sum_Uy) ^ (1/2);

end

