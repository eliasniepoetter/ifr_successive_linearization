function K = controller_synthesis(A,B)

Q = eye(2).*10;
R = 1;
K = lqr(A,B,Q,R);

end