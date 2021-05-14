function val= n_func(i,t,tau)
val = 0;
for j = 1:length(tau)
    if (tau(j) >= t(i) && tau(j) < t(i+1))
        val = val+1;
    end
end
end

