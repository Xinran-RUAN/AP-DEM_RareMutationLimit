<<<<<<< HEAD
function D = initial_D(theta_m, theta)
if theta_m == 0.5
    D = exp(-1 .* sin(pi .* theta).^2); % theta_m = 0.5;
else
    D = exp(-1 .* (sin(pi * theta) + 0.4 * sin(2 * pi * theta)).^2); % theta_m \neq 0.5
end
end

=======
function D = initial_D(theta_m, theta)
if theta_m == 0.5
    D = exp(-1 .* sin(pi .* theta).^2); % theta_m = 0.5;
else
    D = exp(-1 .* (sin(pi * theta) + 0.4 * sin(2 * pi * theta)).^2); % theta_m \neq 0.5
end
end

>>>>>>> d725b20c4e0dc455060136e6e0d3a79bba8f525a
