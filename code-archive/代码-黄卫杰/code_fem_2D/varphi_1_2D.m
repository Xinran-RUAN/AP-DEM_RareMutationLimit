function phi = varphi_1_2D(xi, eta, k)

if k == 1
    phi = 1 - xi - eta;
elseif k == 2
    phi = xi;
elseif k == 3
    phi = eta;
end

end

