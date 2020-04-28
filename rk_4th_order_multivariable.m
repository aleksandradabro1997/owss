function [y, x_vec] = rk_4th_order_multivariable(msh_dens, y_0, t_min, t_max, fun_vec, T_zewn, T_ref, T_refco, T_zco, u)
% x1' = f1(x1, x2, ..., xm)
% x2' = f2(x1, x2, ..., xm)
% ...
% xm' = fm(x1, x2, ..., xm)
%
% x'=f(x)

n = ceil(msh_dens*abs(t_max-t_min));
h = abs(t_max-t_min)/n;
x_vec = t_min:h:t_max;
y = zeros(length(fun_vec),length(x_vec));
y(:,1) = y_0;

for i=1:length(x_vec)-1
    for j=1:length(fun_vec)
        a(j) = h*fun_vec{j}(y(1,i), y(2, i), T_zewn(i), T_zco(i), u(i), x_vec(i+1));
        b(j) = h*fun_vec{j}(y(1,i) + (1/2)*a(j), y(2, i)+ (1/2)*a(j),T_zewn(i), T_zco(i), u(i), x_vec(i+1));
        c(j) = h*fun_vec{j}(y(1,i) + (1/2)*b(j), y(2, i)+ (1/2)*a(j),T_zewn(i), T_zco(i), u(i), x_vec(i+1));
        d(j) = h*fun_vec{j}(y(1,i) + c(j), y(2, i) + c(j),T_zewn(i), T_zco(i), u(i), x_vec(i+1));
    end
    
    y(:,i+1) = y(:,i)+transpose((1/6)*(a + 2*b + 2*c + d));
end
end