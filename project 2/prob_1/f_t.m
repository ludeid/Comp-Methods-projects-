function [val, check] = f_t(t)
if issorted(t)
    check = false;
    val = 1;
    for i = 1:length(t)-1
        val = val*(t(i+1) -t(i));
    end
else
    val = 0;
    check = true;
end
%     for i = 1:length(t)-1
%         if ~(t(i)<t(i+1))
%             val = 0;
%             check = true;
%         end
%     end
end
