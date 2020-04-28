function [A] = generate_A(tau)
% Generuje macierz ograniczeñ na podstawie wektora ograniczen
% 0 <= tau(1) <= tau(2) <= tau(3) <= tau(4) ... <= tau(end) = t_max
vec_1 = -ones(1, length(tau));
vec_2 = ones(1, length(tau)-1);
A_part1 = diag(vec_1);
A_part2 = diag(vec_2, -1);
A = A_part1 + A_part2;

end

