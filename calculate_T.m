function [dT] = calculate_T(t, T_vec, u, T_zewn, Tzco)
% T_vec = [T Ts]
ap = 25.96;
taup = 53.7;
kf = 4.62;
dT = ((T_zewn-T_vec(1))+ ap*(T_vec(2)-T_vec(1)) + kf*u*(Tzco-T_vec(1)))/(taup*t);

end

