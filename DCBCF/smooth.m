function [out] = smooth(in,cell)

[m,n] = size(in);

mid = floor(cell/2) + 1;
k = floor(cell/2);
in_extend = wextend(2,'symh',in,[k,k]);
out = zeros(m,n);

for i = 1 + k : m + k
for j = 1 + k : n + k
        sum = 0;
        for i_count  = 1 : cell
        for j_count = 1 : cell
             sum = sum + in_extend(i + i_count - mid,j+j_count - mid);
        end
        end
        out(i - k , j - k) = sum / cell ^ 2;
end
end

end