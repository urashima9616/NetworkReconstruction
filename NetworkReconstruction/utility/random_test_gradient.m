clear
rand_trial = 20;
params = zeros(2,rand_trial);
params_const = zeros(2,rand_trial);


for cindx = 1:rand_trial
   cindx
   rng(cindx)
   sanity_check_itemized
   params(1,cindx) = gradient(1,1);
   params(2,cindx) = gradient_num(1,1);
   params_const(1,cindx) = gradient_const(1,1);
   params_const(2,cindx) = gradient_const_num(1,1);
end