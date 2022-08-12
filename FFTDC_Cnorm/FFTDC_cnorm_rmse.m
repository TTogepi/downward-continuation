function [x,cnorm,error] = FFTDC_cnorm_rmse(f,f_check,expand,Cell,z,i_taylor,alpha_star,number_alpha,ratio,flag)

cnorm = zeros(1,number_alpha);
x = zeros(1,number_alpha);
error = zeros(1,number_alpha);
table_wait1 = waitbar(0,'The program is running');

for i = 1 : number_alpha
   
    x(i) = alpha_star;
    alpha1 = alpha_star +ratio;
    [U_con1] = FFTDC_main1(f, expand, Cell, z, alpha_star,i_taylor);
    [U_con2] = FFTDC_main1(f, expand, Cell, z, alpha1,i_taylor);
    
    U_con = U_con2 - U_con1;
    cnorm(i) = max(max(abs(U_con)));

    if flag == 1
    error(i) = rms(U_con1,f_check);
    end
    
    alpha_star = alpha_star +ratio;

end

close(table_wait1);
table_wait2 = waitbar(1,'The program completed');
end
