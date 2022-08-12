function [x,cnorm,error] = TRDC_cnorm_rmse(f,f_check,expand,Cell,z,alpha_star,number_alpha,ratio,flag)

cnorm = zeros(1,number_alpha);
x = zeros(1,number_alpha);
error = zeros(1,number_alpha);
table_wait1 = waitbar(0,'The program is running');
for i = 1 : number_alpha
 
    x(i) = alpha_star;
    alpha1 = alpha_star * ratio;
    [U_con1] = TRDC_main1(f, expand,Cell, z,alpha_star);
    [U_con2] = TRDC_main1(f, expand,Cell, z,alpha1);
    U_con = U_con2 - U_con1;
    cnorm(i) = max(max(abs(U_con)));   
    
    if flag == 1
    error(i) = rms(U_con1,f_check);
    end
    
    alpha_star = alpha_star * ratio;
    
end

close(table_wait1);
table_wait2 = waitbar(1,'The program completed');
end
