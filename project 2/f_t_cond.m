function val = f_t_cond(t, lambda, tau)
tempsum = 0;
d = length(lambda);
for i = 1:d
    tempsum = tempsum + lambda(i)*(t(i+1)-t(i));
end
tempprod = 1;
for i = 1:d
    tempprod = tempprod*lambda(i)^n_func(i,t,tau);
end
[f_t_val, zero_check] = f_t(t);
if ~zero_check
    val = exp(-tempsum)*tempprod*f_t_val;
%     if f_t_val == 0 || val == 0
%         disp('__')
%         disp(num2str(f_t_val))
%         disp([val exp(-tempsum) tempprod])
%         disp(lambda)
%     end
else
    val = 0;
    %disp('A')
end
end

