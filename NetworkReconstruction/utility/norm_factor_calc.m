function norm_factor = norm_factor_calc(degree_array, alpha)

if alpha < 0
   temp = degree_array(degree_array > 0);
   norm_factor = sum(temp.^(alpha));
else
    norm_factor = sum(degree_array.^(alpha));
end