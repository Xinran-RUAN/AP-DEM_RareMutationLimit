function w = trapz_weights(nx, dx)
%TRAPZ_WEIGHTS Trapezoidal weights for nx endpoint nodes.
w = dx * ones(nx,1);
if nx >= 2
    w(1) = dx/2;
    w(end) = dx/2;
end
end
