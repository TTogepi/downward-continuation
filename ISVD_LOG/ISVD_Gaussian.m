function [expectation] = ISVD_Gaussian(f,Cell,h,expand,i_taylor,sigma,window)

mode = 'symh';
table_wait1 = waitbar(0,'The program is running');
[m,n] = size(f);
expectation = zeros(m,n);

if mod(window,2) == 0
    errordlg('The input window must be odd','Error') ;
    close(table_wait1);
end

window = floor(window/2);

for i = 1 : i_taylor
    
    ind = i - 1;
    
    if ind == 0
        
        derivative = f;
        
    elseif ind == 1
        
        integ_up = frequency_integration(f,expand,Cell,mode);
        dxx = Gaussian_filter(integ_up,1,Cell,sigma,window);
        dyy = Gaussian_filter(integ_up,2,Cell,sigma,window);
%         dxx =   finite_difference(filter_x,1,Cell);
%         dyy =   finite_difference(filter_y,2,Cell);
        derivative = -(dxx + dyy);
        
        cache_odd_z = derivative; 
     
    elseif ind == 2
        
       dxx = Gaussian_filter(f,1,Cell,sigma,window);
       dyy = Gaussian_filter(f,2,Cell,sigma,window);
%         dxx =   finite_difference(filter_x,1,Cell);
%         dyy =   finite_difference(filter_y,2,Cell);       
        derivative = -(dxx + dyy);
        
        cache_even_z = derivative;  
        
    elseif mod(ind,2) == 0
        
        dxx = Gaussian_filter(cache_even_z,1,Cell,sigma,window);
        dyy = Gaussian_filter(cache_even_z,2,Cell,sigma,window);
%         dxx =  finite_difference(filter_x,1,Cell);
%         dyy =  finite_difference(filter_y,2,Cell);       
        derivative = -(dxx + dyy);
        
        cache_even_z = derivative; 
        
    elseif mod(ind,2) == 1
        
        dxx = Gaussian_filter(cache_odd_z,1,Cell,sigma,window);
        dyy = Gaussian_filter(cache_odd_z,2,Cell,sigma,window);
%         dxx =  finite_difference(filter_x,1,Cell);
%         dyy =  finite_difference(filter_y,2,Cell);       
        derivative = -(dxx + dyy);
        
        cache_odd_z = derivative; 
        
    end
    
    expectation = expectation + h^ind / factorial(ind) * derivative;
    
end

close(table_wait1);
table_wait2 = waitbar(1,'The program completed');

end