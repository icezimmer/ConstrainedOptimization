function [momentum, momentum_coeff] = Momentum(d_old, alpha, beta)
%{
Compute the momentum. The momentum is compute so that the new point must be
    inside the triangular with vertices x (old point), d_old (old direction)
    and d_new (new direction)
INPUT:
    d_old : (vector) old direction
    alpha : (float) step size for the new direction
    beta  : (float) max value for the momentum coefficient
OUTPUT:
    momentum       : (vector) momentum
    momentum_coeff : (float) momentum coefficient
%}

momentum_coeff = min(beta, 1 - alpha);
momentum_coeff = max(0, momentum_coeff);
momentum = momentum_coeff * d_old;

end