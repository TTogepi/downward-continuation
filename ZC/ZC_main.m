function [expectation] = ZC_main(f, expand,Cell, h,k);

table_wait1 = waitbar(0,'The program is running');
if k == 1
    C = [1];
elseif k == 2
    C = [3/2,-1/2];
elseif k == 3
    C = [23/12,-4/3,5/12];
elseif k == 4
    C = [55/24,-59/24,37/24,-3/8];
elseif k == 5
    C = [1901/720,-1387/360,109/30,-637/360,251/720];
elseif k == 6
    C = [4277/1440,-2641/480,4991/720,-3649/720,959/480,-95/288];
elseif k == 7
    C = [198721/60480,-18637/2520,235183/20160,-10754/945,135713/20160,-5603/2520,19087/60480];
elseif k == 8
    C = [16083/4480,-1152169/120960,242653/13440,-296053/13440,2102243/120960,-115747/13440, ...
             32863/13440,-5257/17280];
else
    C = errordlg('The inpurt order must be an integer from 1 to 8','Error') ;
    close(table_wait1);
end

[m,n] = size(f);
expectation = f;
z0 = 0;
f_dz = dz(f,expand,Cell,'symh');

for i = 1 : k
    
    z = z0 - (i-1)*h;
    expectation =  expectation + h*C(i) * up_continue(f_dz,expand,Cell,'symh',-z);
    
end

close(table_wait1);
table_wait2 = waitbar(1,'The program completed');

end