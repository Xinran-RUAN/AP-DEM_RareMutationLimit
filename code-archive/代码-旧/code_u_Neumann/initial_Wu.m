function [W, u] = initial_Wu(theta, x)
%% 
W = exp(1 .* (sin(pi .* theta(1:end-1)) + sin(pi .* x)).^2);
u = 0 .* theta;     
u = u - max(u);   
end

