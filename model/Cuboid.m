function [Hax,Hay,Za,delta_T,Ta]=Cuboid(m_M,m_I,m_A,m_x0,m_y0,m_z0,m_a,m_b,m_c,x,y,z)
% Forward modeling of Cuboid model
% Liu, S., Hu, X., Zhang, H., Geng, M., Zuo, B., 2017. 3D magnetization vector inversion of magnetic data: improving and comparing methods, Pure and Applied Geophysics, 174, 4421-4444.
% m_M   Total magnetization intensity (A/m)
% m_I   Inclination  (Degree)
% m_A   90-Declination
% m_x0  X coordinate of Cuboid center
% m_y0  Y coordinate of Cuboid center
% m_z0  Z coordinate of Cuboid center
% m_a   Length of Cuboid in X direction
% m_b   Length of Cuboid in Y direction
% m_c   Length of Cuboid in Z direction

% Hax   X-component magnetic anomaly
% Hay   Y-component magnetic anomaly
% Za    Z-component magnetic anomaly
% delta_T   Total-field anomaly
% Ta        Amplitude anomaly



I=m_I*pi/180;
A=m_A*pi/180;
If = pi / 2;
Df = 0;
Mx=m_M*cos(I)*cos(A);
My=m_M*cos(I)*sin(A);
Mz=m_M*sin(I);

Hax=0;Hay=0;Za=0;
mu0=4*pi*1.0e-7;

%compute Hax
for i=1:2
    for j=1:2
        for k=1:2
           part1=2*power(-1,k+j+i+1)*...
               atan(((m_z0+power(-1,i+1)*m_c/2-z)+(m_y0+power(-1,j+1)*m_b/2-y)+sqrt(power(m_x0+power(-1,k+1)*m_a/2-x,2)+power(m_y0+power(-1,j+1)*m_b/2-y,2)+power(m_z0+power(-1,i+1)*m_c/2-z,2)))/...
               (m_x0+power(-1,k+1)*m_a/2-x));
           part2=power(-1,k+j+i+1)*...
               log(sqrt(power(m_x0+power(-1,k+1)*m_a/2-x,2)+power(m_y0+power(-1,j+1)*m_b/2-y,2)+power(m_z0+power(-1,i+1)*m_c/2-z,2))+(m_z0+power(-1,i+1)*m_c/2-z));
           part3=power(-1,k+j+i+1)*...
               log(sqrt(power(m_x0+power(-1,k+1)*m_a/2-x,2)+power(m_y0+power(-1,j+1)*m_b/2-y,2)+power(m_z0+power(-1,i+1)*m_c/2-z,2))+(m_y0+power(-1,j+1)*m_b/2-y));
           Hax=Hax+mu0/(4*pi)*(Mx*part1+My*part2+Mz*part3); 
        end
    end
end
Hax=Hax*1.0e9;

%compute Hay
for i=1:2
    for j=1:2
        for k=1:2
            part1=power(-1,k+j+i+1)*...
                log(sqrt(power(m_x0+power(-1,k+1)*m_a/2-x,2)+power(m_y0+power(-1,j+1)*m_b/2-y,2)+power(m_z0+power(-1,i+1)*m_c/2-z,2))+(m_z0+power(-1,i+1)*m_c/2-z));
            part2=2*power(-1,k+j+i+1)*...
                atan(((m_z0+power(-1,i+1)*m_c/2-z)+(m_x0+power(-1,k+1)*m_a/2-x)+sqrt(power(m_x0+power(-1,k+1)*m_a/2-x,2)+power(m_y0+power(-1,j+1)*m_b/2-y,2)+power(m_z0+power(-1,i+1)*m_c/2-z,2)))/...
                (m_y0+power(-1,j+1)*m_b/2-y));
            part3=power(-1,k+j+i+1)*...
                log(sqrt(power(m_x0+power(-1,k+1)*m_a/2-x,2)+power(m_y0+power(-1,j+1)*m_b/2-y,2)+power(m_z0+power(-1,i+1)*m_c/2-z,2))+(m_x0+power(-1,k+1)*m_a/2-x));
            Hay=Hay+mu0/(4*pi)*(Mx*part1+My*part2+Mz*part3);
        end
    end
end
Hay=Hay*1.0e9;

%compute Za
for i=1:2
    for j=1:2
        for k=1:2
            part1=power(-1,k+j+i+1)*...
                log(sqrt(power(m_x0+power(-1,k+1)*m_a/2-x,2)+power(m_y0+power(-1,j+1)*m_b/2-y,2)+power(m_z0+power(-1,i+1)*m_c/2-z,2))+(m_y0+power(-1,j+1)*m_b/2-y));
            part2=power(-1,k+j+i+1)*...
                log(sqrt(power(m_x0+power(-1,k+1)*m_a/2-x,2)+power(m_y0+power(-1,j+1)*m_b/2-y,2)+power(m_z0+power(-1,i+1)*m_c/2-z,2))+(m_x0+power(-1,k+1)*m_a/2-x));
            part3=2*power(-1,k+j+i+1)*...
                atan(((m_y0+power(-1,j+1)*m_b/2-y)+(m_x0+power(-1,k+1)*m_a/2-x)+sqrt(power(m_x0+power(-1,k+1)*m_a/2-x,2)+power(m_y0+power(-1,j+1)*m_b/2-y,2)+power(m_z0+power(-1,i+1)*m_c/2-z,2)))/...
                (m_z0+power(-1,i+1)*m_c/2-z));
            Za=Za+mu0/(4*pi)*(Mx*part1+My*part2+Mz*part3);
        end
    end
end
Za=Za*1.0e9;

%compute delta_T
delta_T=Hax*cos(If)*cos(Df)+Hay*cos(If)*sin(Df)+Za*sin(If);
 
%compute Ta
Ta=sqrt(Hax*Hax+Hay*Hay+Za*Za);