function dp = d_varphi_1_2D(xi, eta, k, ol)
% ol 表示对哪个变量求偏导，xi 为1， eta 为2
if ol == 1
    if k == 1
        dp = - 1;
    elseif k == 2
        dp = 1;
    elseif k == 3
        dp = 0;
    end
elseif ol == 2
    if k == 1
        dp = -1;
    elseif k == 2
        dp = 0;
    elseif k == 3
        dp = 1;
    end
end

end

