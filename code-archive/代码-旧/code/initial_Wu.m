<<<<<<< HEAD
function [W, u] = initial_Wu(theta, x)

%W = ones(N_x+1, N_theta); % theta \in [0, 1-dtheta]
% u = (sin(pi * theta) - 1);
%u = exp(1 .* (sin(pi * theta) + 0. * sin(2*pi*theta)).^2);
W = exp(1 .* (sin(pi .* theta) + sin(pi .* x)).^2);
u = 0 .* theta;     
u = u - max(u);   
end

=======
function [W, u] = initial_Wu(theta, x)

%W = ones(N_x+1, N_theta); % theta \in [0, 1-dtheta]
% u = (sin(pi * theta) - 1);
%u = exp(1 .* (sin(pi * theta) + 0. * sin(2*pi*theta)).^2);
W = exp(1 .* (sin(pi .* theta) + sin(pi .* x)).^2);
u = 0 .* theta;     
u = u - max(u);   
end

>>>>>>> d725b20c4e0dc455060136e6e0d3a79bba8f525a
