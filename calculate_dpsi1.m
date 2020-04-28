function [dpsi1] = calculate_dpsi1(t, T_vec, u, T_ref)
% T_vec = [T Ts psi1 psi2]
ap = 25.96;
taup = 53.7;
as = 6.47;
taus = 153.54;
kf = 4.62;
dpsi1 = ((1+ap+kf*u)/(taup*t))*T_vec(3) - (as/(taus*t))*T_vec(4) + 2*T_vec(1)-2*T_ref;

end
