function K = controller_synthesis(A,B)

Q = eye(2).*1000;
R = 1;
K = lqr(A,B,Q,R);

end