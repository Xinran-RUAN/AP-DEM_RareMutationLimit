function B = prepare_part(eps, dt, d_theta, N_theta)

%% solve u µÄÏµÊı¾ØÕó
beta = eps * dt / d_theta^2;  
B = (1 + 2 * beta) * diag(ones(N_theta, 1), 0) +...
            - beta * diag(ones(N_theta-1, 1), 1) +...
            - beta * diag(ones(N_theta-1, 1), -1);
B(1, N_theta) = - beta;
B(N_theta, 1) = - beta;   

end

