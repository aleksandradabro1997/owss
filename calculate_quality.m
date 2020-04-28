function [quality] = calculate_quality(T, T_ref)
    quality = (1/2)*integral((T-T_ref)^2)
end

