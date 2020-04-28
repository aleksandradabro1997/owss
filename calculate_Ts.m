function [dTs] = calculate_Ts(t, T_vec, T_zewn)
% T_vec = [T Ts]
    as = 6.47;
    taus = 153.54;
    dTs = ((T_zewn-T_vec(2)) + as*(T_vec(1)-T_vec(2)))/(taus*t);
end

