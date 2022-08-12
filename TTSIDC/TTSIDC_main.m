 function [expectation]=TTSIDC_main(f,Cell,z,expand,sigma,window,taylor_down,taylor_up,iter)
 
 % truncated Taylor series iterative downward continuation (TTSIDC)
 
%  %********************************************************
%  iter = 50;
%  %********************************************************
mode = 'symh';
table_wait1 = waitbar(0,'The program is running');
if mod(window,2) == 0
    errordlg('The input window must be odd','Error') ;
    close(table_wait1);
end
%***********************³õÖµ*****************************
expectation = ISVD_Gaussian_for_TTSIDC(f,Cell,z,expand,taylor_down,sigma,window);
%********************************************************
 for i = 1 : iter
     
     up = up_continue(expectation,expand,Cell,mode,z);
     residue = f - up;
     dt =  ISVD_Gaussian_for_TTSIDC(residue,Cell,z,expand,taylor_up,sigma,window);
     expectation = expectation + dt;
    
 end
close(table_wait1);
table_wait2 = waitbar(1,'The program completed');
 end