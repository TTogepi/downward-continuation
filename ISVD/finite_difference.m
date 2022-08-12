function [out] = finite_difference(U,flag,cell,k)

[m,n] = size(U);

out = zeros(m,n);
f = wextend(2,'symh',U,[k,k]);

if k == 1
    coff = [1,-2,1];
elseif k == 2
    coff = [-1/12, 4/3, -5/2, 4/3, -1/12];
elseif k == 3
    coff = [ 1/90, -3/20, 3/2, -245/90, 3/2, -3/20, 1/90];
elseif k == 4
    coff = [-1/560, 8/315, -1/5, 8/5, -1453/504, 8/5, -1/5, 8/315, -1/560];
end

%%
for i = 1 : m
for j = 1 : n
    
    if flag == 1

        for ii = 1 : 2 * k + 1
            out(i,j) = out(i,j) + coff(ii) * f(i+ii-1,j+k);
        end
        out(i,j) =  out(i,j) *  (1 / cell ^ 2); 
        
    else
        
        for ii = 1 : 2 * k + 1
            out(i,j) = out(i,j) + coff(ii) * f(i+k,j+ii-1);
        end
        out(i,j) =  out(i,j) *  (1 / cell ^ 2); 
    end
    
end
end

%%

%  k = 1      1, -2, 1
%  k = 2      -1/12, 4/3, -5/2, 4/3, -1/12
%  k = 3      1/90, -3/20, 3/2, -245/90, 3/2, -3/20, 1/90
%  k = 4      -1/560, 8/315, -1/5, 8/5, -1453/504, 8/5, -1/5, 8/315, -1/560

end