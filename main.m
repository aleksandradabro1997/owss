% Inicjalizacja danych dla modelu ogrzewania 
%{
Tzco - temperatura zasilania
T - temperatura w pomieszczeniu
Ts - temperatura œciany
Tzewn - temperatura zewnêtrzna
u - stopieñ otwarcia zaworu w termostacie <0,1>
Tref - temperatura zadana na termostacie
Trefco - temperatura referencyjna dla wêz³a CO
Tsampl - okres próbkowania
%}
%%
%{
    mshdens = 135 kroków/j.czasu
    delta_t -> ceil(mshdens * delta_t) = n -> h = delta_t/n
    Pozniej trzeba podawac wektro czasów i prze³¹czen
%}
%% Generowanie danych
clear all; close all;
ap = 25.96;
as = 6.47;
kf = 4.62;
taus = 153.54;
taup = 53.7;

[T_zewn ,name] = xlsread('Temperatura');
T_ref = zeros(1, length(T_zewn));
T_refco = zeros(1, length(T_zewn));
T_zco = zeros(1, length(T_zewn));
time = zeros(1, length(T_zewn));
u = zeros(1, length(T_zewn));
T_zewn = T_zewn+10;
for i=2:length(T_zewn) 
    time(i) = i-1;
    if i<500 || i>1000 
        T_ref(i) = 12;
        %u(i) = 1;
    else
        T_ref(i) = 8;
        %u(i) = 0;
    end
    T_refco(i) = 30;
    T_zco(i) = 80;
    
end

msh_dens = 10;
%% Ograniczenia funkcji celu
% Funkcja celu zale¿y od czasów prze³¹czeñ
% tau - wektor czasów prze³¹czeñ - zmiana sterowania
% tau(0) <tau(1) < tau(2) < ....
%tau = [80];
tau = [60 90];
F_celu = zeros(1,length(tau));
t_min = 0.1;
t_max = 200;
x0 = [0;0];
%% Zadawanie wektora sterowañ
n = ceil(msh_dens*abs(t_max-t_min));
u(1) = 1;

zmiana_sterowania = tau;
zmiana_sterowania = sort(zmiana_sterowania)*msh_dens;
% Generowanie wektora sterowania
for i=2:length(T_zewn)
    if ismember(i, zmiana_sterowania)
        u(i) = ~u(i-1);
    else
        u(i) = u(i-1);
    end
end
    %u = rand(length(T_zewn),1);
    %u = round(u);
    %[t, x, J] = rk_4th_order(msh_dens, [0;0], T_ref, T_zewn, T_zco, t_min, t_max, tau);
    %F_celu(j) = J;


%% Optymalizacja - minimalizacja wskaŸnika jakoœci
% Ograniczenia
% 0 <= tau(1) <= tau(2) <= tau(3) <= tau(4) ... <= tau(end) = t_max
cost_fun = @(tau) rk_4th_order(msh_dens, x0, T_ref, T_zewn, T_zco, t_min, t_max, tau);
n_tau = length(tau);
A = generate_A(tau);
b = zeros(n_tau, 1);
%options = optimset('display', 'iter');
%options = optimoptions('fmincon');
options = optimset('PlotFcns','optimplotfval','TolX',1e-6, 'TolFun', 1e-6,...
                   'display', 'iter', 'FinDiffRelStep', 15);
tauopt = fmincon(cost_fun, tau, A, b, [], [], [], [], [], options);
tauopt2 = fmincon(cost_fun, tauopt, A, b, [], [], [], [], [], options);
%% Wykres dla znalezionych optymalnych wartoœci
[t, x, J] = rk_4th_order2(msh_dens, x0, T_ref, T_zewn, T_zco, t_min, t_max, tau);
figure();hold on;grid on;
plot(t, x(:,1));
plot(t, x(:,2));
plot(t, T_ref(1:length(t)))
for j=1:length(tau)
    xline(tau(j));
end
legend('T', 'Ts', 'T ref') 
title('Wykres temperatur - zadane tau');
xlabel('Czas [s]');

[t, x, J] = rk_4th_order2(msh_dens, x0, T_ref, T_zewn, T_zco, t_min, t_max, tauopt);
figure();hold on;grid on;
plot(t, x(:,1));
plot(t, x(:,2));
plot(t, T_ref(1:length(t)))
for j=1:length(tauopt)
    xline(tauopt(j));
end
legend('T', 'Ts', 'T ref') 
title('Wykres temperatur - tau - fmincon 1');
xlabel('Czas [s]');

[t, x, J] = rk_4th_order2(msh_dens, x0, T_ref, T_zewn, T_zco, t_min, t_max, tauopt2);
figure();hold on;grid on;
plot(t, x(:,1));
plot(t, x(:,2));
plot(t, T_ref(1:length(t)))
for j=1:length(tauopt2)
    xline(tauopt2(j));
end
legend('T', 'Ts', 'T ref') 
title('Wykres temperatur - zadane tau - fmincon 2');
xlabel('Czas [s]');

figure();grid on; hold on;
plot(t, u(1:t_max*msh_dens-1));
title('Sterowanie optymalne');
xlabel('Czas [s]');
ylabel('Sterowanie 0, 1');
%% Sprawdzenie czy sterowanie jest optymalne
% J - jakobian
%J = [(-1-ap-kf*u)/taup ap/taup; as/taus (-1-as)/taus];
% psi' = -J'psi <- psi = [psi(1) ;psi(2)]
% x = [T Ts psi1 psi2]
% rownania_T = @(t, x, ap, as, taup, taus, kf, T_zewn, T_zco) [(T_zewn(t)-x(1))/taup + ap*(x(2)-x(1))/taup + kf*(T_zco(t) - x(1))*u(t)/taup;...
%                                                           (T_zewn(t)-x(2))/taus + as*(x(1)-x(2))/taus];
% rownania_sprzezone = @(t, x, ap, as, taup, taus, kf) [((1+ap+kf*u)/taup)*x(3) - (as/taus)*x(4); ...
%                                                       ((1+as)/taus)*x(4) - (ap/taup)*x(3)];
% rownania_sprzezone_2 = @(t, x, ap, as, taup, taus, kf, T_ref) [((1+ap+kf*u)/taup)*x(3) - (as/taus)*x(4) + 2*x(1)-2*T_ref(t); ...
%                                                       ((1+as)/taus)*x(4) - (ap/taup)*x(3)];
% 
% rownania = @(t,x,ap,as,taup,taus,kf,T_zewn,T_zco,T_ref, u) [(T_zewn(t)-x(1))/taup + ap*(x(2)-x(1))/taup + kf*(T_zco(t) - x(1))*u(t)/taup;...
%                                                          (T_zewn(t)-x(2))/taus + as*(x(1)-x(2))/taus; ...
%                                                          ((1+ap+kf*u(t))/taup)*x(3) - (as/taus)*x(4) + 2*x(1)-2*T_ref(t); ...
%                                                          ((1+as)/taus)*x(4) - (ap/taup)*x(3)];
                                                     
%fun = @(t, x) rownania(t,x,ap,as,taup,taus,kf,T_zewn,T_zco,T_ref, u);
                                                     
% Rozwi¹zywanie wstecz
T_min = t_max;
T_max = t_min;
dt = 1/msh_dens;
T_vec = T_min:-dt:T_max;
x_0 = [x0; x(end,1); x(end,2)];
y = x_0;
[t, y] = rk_4th_order_backwards(msh_dens, x_0, T_ref, T_zewn, T_zco, T_min, T_max, tauopt2);
%%
figure;hold on;grid on;
xlabel('T'); ylabel('Ts');
plot(y(:,1), y(:,2));
plot(y(end,1), y(end,2), 'ro');
%%
figure;hold on; grid on;
title('PSI')
xlabel('psi1'); ylabel('psi2');
plot(y(:,3), y(:,4));
plot(y(end,3), y(end,4), 'ro');
%% Hamiltonian
% Co to macierz B?
%ham = @(t,x,u, T_ref) x(3:4)'*[(-1-ap-kf*u)/taup ap/taup; as/taus (-1-as)/taus]*x(1:2)+x(3:4)'*[2*x(1)-2*T_ref; 0]*u;

ham = @(t,x,u, T_ref) (x(t,3:4))*[(-1-ap-kf*u(t))/taup ap/taup; as/taus (-1-as)/taus]*(x(t,1:2)')+x(t,3:4)*[2*x(t,1)-2*T_ref(t); 0]*u(t);
ham_val = zeros(1, length(y));
for i=1:length(y)
    ham_val(i) = ham(i, y, u, T_ref);
end

figure(); hold on; grid on;
title('Hamiltonian')
plot(ham_val);