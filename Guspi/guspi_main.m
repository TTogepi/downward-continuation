 function [expectation]=guspi_main(f,expand,Cell,z,iter,itaylor)
 
 [m,n] = size(f);
 mode = 'symh';
 expectation = zeros(m,n);
 stack = zeros(m,n);
 
%********************************************************
 expectation = f;
%********************************************************

 table_wait1 = waitbar(0,'The program is running');
 for i = 1 : iter
     
     stack = cal_Tay_guspi(expectation,expand,Cell,mode,z,itaylor);
     expectation = f - stack;
     
     
 end
 
close(table_wait1);
table_wait2 = waitbar(1,'The program completed');

 end