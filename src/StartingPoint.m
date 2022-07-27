function x_start = StartingPoint(P)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[K, n] = size(P);
x_start = zeros(n, 1);
for k = 1 : K
    i = find(P(k,:), 1);
    x_start(i) = 1;
end

end