function[delta_all] = combine_cub()

m_M=5;     %total magnetization intensity (A/m)
m_I= 90;     % Inclination  (Degree)
m_A=90;      % 90-Declination

m_x0=1400;   %X coordinate of Cuboid center
m_y0=1200;   %Y coordinate of Cuboid center
m_z0=500;   %Z coordinate of Cuboid center
m_a=800;    %length of Cuboid in X direction
m_b=150;    %length of Cuboid in Y direction
m_c=50;    %length of Cuboid in Z direction

m_x1=1300;   %X coordinate of Cuboid center
m_y1=1800;   %Y coordinate of Cuboid center
m_z1=500;   %Z coordinate of Cuboid center
m_a1=1000;    %length of Cuboid in X direction
m_b1=300;    %length of Cuboid in Y direction
m_c1=50;    %length of Cuboid in Z direction

m_x2=1300;   %X coordinate of Cuboid center
m_y2=2300;   %Y coordinate of Cuboid center
m_z2=500;   %Z coordinate of Cuboid center
m_a2=1000;    %length of Cuboid in X direction
m_b2=300;    %length of Cuboid in Y direction
m_c2=50;    %length of Cuboid in Z direction

m_x3=2400;   %X coordinate of Cuboid center
m_y3=2000;   %Y coordinate of Cuboid center
m_z3=600;   %Z coordinate of Cuboid center
m_a3=150;    %length of Cuboid in X direction
m_b3=1500;    %length of Cuboid in Y direction
m_c3=100;    %length of Cuboid in Z direction

C1=100;
C2=100;
cell = 40;


delta0 = zeros(C1,C2);
delta1 = zeros(C1,C2);
delta2 = zeros(C1,C2);

for i=0:C1
    for j=0:C2
        [Hax,Hay,Za,delta_T,Ta]=...
            Cuboid(m_M,m_I,m_A,m_x0,m_y0,m_z0,m_a,m_b,m_c,cell*j,cell*i,0);
        delta0(i+1,j+1) = delta_T;
    end
end

for i=0:C1
    for j=0:C2
        [Hax,Hay,Za,delta_T,Ta]=...
            Cuboid(m_M,m_I,m_A,m_x1,m_y1,m_z1,m_a1,m_b1,m_c1,cell*j,cell*i,0);
        delta1(i+1,j+1) = delta_T;
    end
end

for i=0:C1
    for j=0:C2
        [Hax,Hay,Za,delta_T,Ta]=...
            Cuboid(m_M,m_I,m_A,m_x2,m_y2,m_z2,m_a2,m_b2,m_c2,cell*j,cell*i,0);
        delta2(i+1,j+1) = delta_T;
    end
end

for i=0:C1
    for j=0:C2
        [Hax,Hay,Za,delta_T,Ta]=...
            Cuboid(m_M,m_I,m_A,m_x3,m_y3,m_z3,m_a3,m_b3,m_c3,cell*j,cell*i,0);
        delta3(i+1,j+1) = delta_T;
    end
end


% subplot(1,2,1)
plot_cuboid(m_x0,m_y0,-m_z0,m_a,m_b,m_c);
plot_cuboid(m_x1,m_y1,-m_z1,m_a1,m_b1,m_c1);
plot_cuboid(m_x2,m_y2,-m_z2,m_a2,m_b2,m_c2);
plot_cuboid(m_x3,m_y3,-m_z3,m_a3,m_b3,m_c3);
view(0,0);
xlabel('Easting (m)');
zlabel('Depth (m)');

delta_all = delta0 + delta1 + delta2 +delta3 ;
% subplot(1,2,2)
% contourf(delta_all,10);
colorbar;
end