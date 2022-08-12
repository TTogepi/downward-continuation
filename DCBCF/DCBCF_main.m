function [expectation] = DCBCF_main(f, expand,Cell, h,L,window1,window2)

[m,n] = size(f);
expectation = zeros(m,n);
i_Taylor = L;
b = zeros(L,m,n);
c = zeros(L,m,n);
mode = 'symh';

Ci0 = cal_Tay(f,expand,Cell,mode,i_Taylor);
C1i = zeros(L-1,m,n);
Ci1 = zeros(L-1,m,n);

table_wait1 = waitbar(0,'The program is running');

if mod(window1,2) || mod(window2,2)== 0
    errordlg('The input windows must be odd','Error') ;
    close(table_wait1);
end

for i = 2 : L - 1
    Ci1(i,:,:)= -Ci0(i+1,:,:) ./ Ci0(2,:,:);
end

for i = 1 : L
   
    if i == 1       
        b(i,:,:) = Ci0(1,:,:);        
    end
    
    if i == 2
        b(i,:,:) = 1 ./ Ci0(2,:,:);
    end 
        
    if i >= 3
        for j = i - 2  : i - 1       
                 C1i(j,:,:) = loop(2,j,b,Ci0,Ci1);
        end
        b(i,:,:) = C1i(i-2,:,:) ./ C1i(i-1,:,:); 
    end
    
end 

 for i = 1 : L
         for k = 1 : m
         for l = 1 : n
             temp1(k,l) = b(i,k,l);
         end
         end
         c(i,:,:) = smooth(temp1,window1);
 end

temp = sum_loop(c,1,L,h);

for i = 1 : m
for j = 1 : n
    
    expectation(i,j) = temp(1,i,j);
    
end
end

expectation = medfilt2(expectation,[window2,window2]);

close(table_wait1);
table_wait2 = waitbar(1,'The program completed');

end