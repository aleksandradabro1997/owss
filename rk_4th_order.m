function [J] = rk_4th_order(msh_dens, x_0, T_ref, T_zewn, T_zco,t_min, t_max, tau)
%% rk_4th_order 
% Inputs
% 1. msh_dens - gêtoœæ siatki
% 2. x0 - punkt startowy
% 3. T_ref - temperatura zadana
% 4. T_zewn - temperatura zewnêtrzna
% 5. T_zco - temperatura zasilania z CO
% 6. t_min - czas startowy
% 7. t_max - czas koñcowy
% 8. tau - wektor czasów prze³¹czeñ sterowania

% dx1 = h*f(x,y)
% dx2 = h*f(x+h/2, y+(1/2)dx1)
% dx3 = h*f(x+h/2, y+(1/2)dx2)
% dx4 = h*f(x+h/2, y+dx3)
%% Inicjalizacja
u = zeros(1, length(T_zewn));
u(1) = 1;
tau = round(tau*msh_dens);
%tau = tau*msh_dens;
for i=2:length(T_zewn)
    if ismember(i, tau)
         u(i) = ~u(i-1);
    else
         u(i) = u(i-1);
    end
end
n = ceil(msh_dens*abs(t_max-t_min));
h = abs(t_max-t_min)/n;
t = t_min;
var_nb = length(x_0);
x = zeros(n, var_nb);
x(1,:) = x_0';
dx1 = zeros(var_nb, 1);
dx2 = dx1;dx3 = dx1;dx4 = dx1;
tmp = zeros(n,1);xtmp=x_0;
J = 0;
%% Rozwi¹zywanie równañ
for i=1:n
    dx1(1)=calculate_T(t,xtmp,u(i),T_zewn(i), T_zco(i));
    dx1(2)=calculate_Ts(t, xtmp,T_zewn(i));
    tmp=xtmp+(h/2)*dx1; 
    %tmp=xtmp+(1/2)*dx1; 
    t=t+(h/2);
    dx2(1)=calculate_T(t,tmp,u(i),T_zewn(i), T_zco(i));
    dx2(2)=calculate_Ts(t,tmp,T_zewn(i));
    tmp=xtmp+(h/2)*dx2;
    %tmp=xtmp+(1/2)*dx2;
    dx3(1)=calculate_T(t,tmp,u(i),T_zewn(i), T_zco(i));
    dx3(2)=calculate_Ts(t,tmp,T_zewn(i));
    tmp=xtmp+h*dx3;
    %tmp=xtmp+dx3;
    t=t+(h/2);
    dx4(1)=calculate_T(t,tmp,u(i),T_zewn(i), T_zco(i));
    dx4(2)=calculate_Ts(t,tmp,T_zewn(i));
    xtmp=xtmp+(h/6)*(dx1+2*dx2+2*dx3+dx4);
    x(i,:)=xtmp';
    J = J + (T_ref(i)-x(i,1))^2;
    if ismember(i, tau)
        %pause()
    end
end
t = linspace(t_min, t_max, n);
J = J/2;
end


