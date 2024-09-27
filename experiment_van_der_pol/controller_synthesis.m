function K = controller_synthesis(A,B)

% weighting matrices for the LQR
Q = eye(2).*100;
R = 1;

% solver call for the ricatti equation
K = lqr(A,B,Q,R);

end