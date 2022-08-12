function [out] = new_kernel(alpha,beta,delta_z,xm,yn,dx,dy)

part1 = (alpha + dx / 2 - xm) * (beta + dy / 2 -yn) / (delta_z * ((alpha + dx / 2 - xm) ^2 + ...   
            (beta + dy / 2 -yn)^2 + (delta_z)^2) ^ (1/2));
part2 = (alpha + dx / 2 - xm) * (beta - dy / 2 -yn) / (delta_z * ((alpha + dx / 2 - xm) ^2 + ...   
            (beta - dy / 2 -yn)^2 + (delta_z)^2) ^ (1/2));      
part3 = (alpha - dx / 2 - xm) * (beta + dy / 2 -yn) / (delta_z * ((alpha - dx / 2 - xm) ^2 + ...   
            (beta + dy / 2 -yn)^2 + (delta_z)^2) ^ (1/2));
part4 = (alpha - dx / 2 - xm) * (beta - dy / 2 -yn) / (delta_z * ((alpha - dx / 2 - xm) ^2 + ...   
            (beta - dy / 2 -yn)^2 + (delta_z)^2) ^ (1/2));

out = atan(part1) - atan(part2) - atan(part3) + atan(part4);

end