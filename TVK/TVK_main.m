function [gdown] = TVK_main(f, expand,Cell, h,n_Taylor)

nn = n_Taylor;
dh = h / nn;
table_wait1 = waitbar(0,'The program is running');
[m,n] = size(f);
field_i = zeros(n_Taylor,m,n);
field_0 = zeros(m,n);
gdown = zeros(m,n);
if nn ==3
    cof = [4,-6,4,-1];
elseif nn ==5
    cof = [6,-15,20,-15,6,-1];
elseif nn ==8
    cof = [9,-36,84,-126,126,-84,36,-9,1];
elseif nn == 10
    cof = [11,-55,165,-330,462,-462,330,-165,55,-11,1];
end

for i = 1 : n_Taylor 
    
    z = i * dh;
    field_i(i,:,:)= up_continue(f,expand,Cell,'symh',z);
    
end 

field_0 = f;

for k = 1 : nn
    
    for i = 1 : m
    for j = 1 : n

    gdown(i,j) =  cof(1) * field_0(i,j);
    for s = 1 : nn
    gdown(i,j) = gdown(i,j) + cof(s+1) * field_i(s,i,j);
    end

    for s = 1 : nn-1
    field_i(nn-s+1,i,j) = field_i(nn-s,i,j);
%     field_i(8,i,j) = field_i(7,i,j);
%     field_i(7,i,j) = field_i(6,i,j);
%     field_i(6,i,j) = field_i(5,i,j);       
%     field_i(5,i,j) = field_i(4,i,j);    
%     field_i(4,i,j) = field_i(3,i,j);
%     field_i(3,i,j) = field_i(2,i,j);
%     field_i(2,i,j) = field_i(1,i,j);
    end
    field_i(1,i,j) = field_0(i,j);
    field_0(i,j) = gdown(i,j);
    
    end
    end
    
end

close(table_wait1);
table_wait2 = waitbar(1,'The program completed');

end

%%coefficient of n_Taylor
% N = 3    4    -6    4    -1
% N = 5    6    -15    20    -15    6    -1
% N = 8    9    -36    84    -126    126    -84    36    -9    1
% N = 10  11    -55    165    -330    462    -462    330    -165    55    -11    1

