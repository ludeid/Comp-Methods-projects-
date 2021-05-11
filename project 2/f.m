function val= f(t)
    for i = 1:length(t)-1
        if t(i)>=t(i+1)
            val = 0;
            return
        end
    end
    
    val = 1;
    for i = 1:length(t)-1
        val = val*(t(i+1) -t(i));
    end
end
