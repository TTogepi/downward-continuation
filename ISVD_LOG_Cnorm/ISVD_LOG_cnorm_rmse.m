function [x,cnorm,error] = ISVD_LOG_cnorm_rmse(f,f_check,Cell,h,expand,i_taylor,window,sigma_star,number_sigma,ratio,flag)

x = zeros(1,number_sigma);
cnorm = zeros(1,number_sigma);
error = zeros(1,number_sigma);
table_wait1 = waitbar(0,'The program is running');

if mod(window,2) == 0
    errordlg('The input window must be odd','Error') ;
    close(table_wait1);
end

for i = 1 : number_sigma
    
    x(i) = sigma_star;
    sigma1 = sigma_star + ratio;
    [U_con1] = ISVD_Gaussian1(f,Cell,h,expand,i_taylor,sigma_star,window);
    [U_con2] = ISVD_Gaussian1(f,Cell,h,expand,i_taylor,sigma1,window);
    U_con = U_con2 - U_con1;
    cnorm(i) = max(max(abs(U_con)));
    
    %
    if flag == 1
    error(i) = rms(U_con1,f_check);
    end
   
    sigma_star = sigma_star + ratio;

end

close(table_wait1);
table_wait2 = waitbar(1,'The program completed');
end