function [xi, eta, weight] = numerical_intergal(k)

if k == 22
    xi = [0, 1, 0];
    eta = [0, 0, 1];
    weight = [1/6, 1/6, 1/6];
end

if k == 23
    xi = [0, 0, 1, 1/2, 1/2, 0];
    eta = [0, 1, 0, 0, 1/2, 1/2];
    weight = [0, 0, 0, 1/6, 1/6, 1/6];
end

end

