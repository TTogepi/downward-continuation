 function [grad,fit,e_rms]=TTSIDC_fs_rmse(f,f_check,Cell,z,expand,sigma,window,taylor_down,taylor_up,iter,flag)
 
 % truncated Taylor series iterative downward continuation (TTSIDC)
 
%  %********************************************************
%  iter = 50;
%  %********************************************************
mode = 'symh';
grad = zeros(1,iter);
fit = zeros(1,iter);
e_rms = zeros(1,iter);

table_wait1 = waitbar(0,'The program is running');
if mod(window,2) == 0
    errordlg('The input window must be odd','Error') ;
    close(table_wait1);
end
%**************************************************************************************
expectation = ISVD_Gaussian_for_TTSIDC(f,Cell,z,expand,taylor_down,sigma,window);
%**************************************************************************************
 for i = 1 : iter
     
     up = up_continue(expectation,expand,Cell,mode,z);
     residue = f - up;
     dt =  ISVD_Gaussian_for_TTSIDC(residue,Cell,z,expand,taylor_up,sigma,window);
     expectation = expectation + dt;
     
     up2 = up_continue(expectation,expand,Cell,mode,z);
     grad(i) = grad_U_l2(expectation,Cell);
     fit(i) = sum(sum((up2 - f).^2)) / sum(sum(f.^2));
     fit(i) = fit(i) ^ (1/2);
     
     if flag == 1
     e_rms(i) = rms(expectation,f_check);
     end
     
 end
close(table_wait1);
table_wait2 = waitbar(1,'The program completed');
 end