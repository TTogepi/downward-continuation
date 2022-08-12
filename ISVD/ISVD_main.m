function [expectation] = ISVD_main(f,Cell,h,expand,i_taylor,k)

mode = 'symh';
table_wait1 = waitbar(0,'The program is running');
[m,n] = size(f);
expectation = zeros(m,n);

if k ~= 1 || k ~= 2  || k ~= 3 || k ~= 4
    errordlg('The input k should not exceed 4 and must be a minimum of 1. It must be an integer','Error') ;
    close(table_wait1);
end

for i = 1 : i_taylor
    
    ind = i - 1;
    derivative = zeros(m,n);
    
    if ind == 0
        
        derivative = f;
        
    elseif ind == 1
        
        integ_up = frequency_integration(f,expand,Cell,mode);
        dxx =  finite_difference(integ_up,1,Cell,k);
        dyy =  finite_difference(integ_up,2,Cell,k);
        derivative = -(dxx + dyy);
        
        cache_odd_z = derivative; 
     
    elseif ind == 2
        
        dxx =  finite_difference(f,1,Cell,k);
        dyy =  finite_difference(f,2,Cell,k);       
        derivative = -(dxx + dyy);
        
        cache_even_z = derivative;  
        
    elseif mod(ind,2) == 0
        
        dxx =  finite_difference(cache_even_z,1,Cell,k);
        dyy =  finite_difference(cache_even_z,2,Cell,k);       
        derivative = -(dxx + dyy);
        
        cache_even_z = derivative; 
        
    elseif mod(ind,2) == 1
        
        dxx =  finite_difference(cache_odd_z,1,Cell,k);
        dyy =  finite_difference(cache_odd_z,2,Cell,k);       
        derivative = -(dxx + dyy);
        
        cache_odd_z = derivative; 
        
    end
    
    expectation = expectation + h^ind / factorial(ind) * derivative;
    
end

close(table_wait1);
table_wait2 = waitbar(1,'The program completed');

end