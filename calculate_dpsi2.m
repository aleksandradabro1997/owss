function [dpsi2] = calculate_dpsi2(t, T_vec)
% T_vec = [T Ts psi1 psi2]
ap = 25.96;
taup = 53.7;
as = 6.47;
taus = 153.54;
kf = 4.62;
dpsi2 = ((1+as)/(taus*t))*T_vec(4) - (ap/(taup*t))*T_vec(3);

end