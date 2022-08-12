function [temp] = sum_loop(b,x,L,h)

if L ~=x
    temp = b(x,:,:) + h ./  sum_loop(b,x+1,L,h);
else
    temp = b(x,:,:);
end
    
end