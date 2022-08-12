function [C] = loop(i,L,b,C1,C2)

if L == 1
    C = C1(i,:,:);
elseif L ==2
    C = C2(i,:,:);
else
    C = loop(i+1,L-2,b,C1,C2) - b(L,:,:) .* loop(i+1,L-1,b,C1,C2);
end

end
