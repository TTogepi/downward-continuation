 function [xl2,fit,e_rms]=Guspi_L_rmse(f,f_check,expand,Cell,z,iter,itaylor,flag)
 
%  %********************************************************
%  iter = 50;
%  %********************************************************
 [m,n] = size(f);
 mode = 'symh';
 expectation = zeros(m,n);
 stack = zeros(m,n);
 xl2 = zeros(1,iter);
 fit = zeros(1,iter);
 e_rms = zeros(1,iter);
 
 table_wait1 = waitbar(0,'The program is running');
%*******************************************************
 expectation = f;
%********************************************************
 
 for i = 1 : iter
     
     stack = cal_Tay_guspi(expectation,expand,Cell,mode,z,itaylor);
     expectation = f - stack;
     
     up = up_continue(expectation,expand,Cell,mode,z);
     
     xl2(i) = sum(sum(expectation.^2));
     xl2(i) = xl2(i) ^ (1/2);
     fit(i) = sum(sum((up - f).^2)) / sum(sum(f.^2));
     fit(i) = fit(i) ^ (1/2);
     
     if flag == 1
     e_rms(i) = rms(expectation,f_check);
     end
     
 end
  
close(table_wait1);
table_wait2 = waitbar(1,'The program completed');

 end