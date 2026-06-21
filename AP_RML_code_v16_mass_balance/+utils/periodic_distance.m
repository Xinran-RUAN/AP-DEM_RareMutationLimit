function d = periodic_distance(a, b)
%PERIODIC_DISTANCE Distance on T=[0,1).
d = abs(a - b);
d = min(d, 1 - d);
end
