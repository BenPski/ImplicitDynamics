function [g,xi,eta,mu,lambda] = initDynamics(N)
    N = double(N);
    L = 10e-2;
    xi = repmat([0;pi/4/L;0;0;0;1],1,N);
    eta = zeros(6,N);
    g = zeros(16,N);
    g = zeros(12,N);
    G = eye(4);
    %g(:,1) = reshape(G,1,16);
    g(:,1) = flatten(G);
%     L = 10e-2;
    ds = L/(N-1);
    for i=2:N
        G = G*expm(ds*se(xi(:,i-1)));
        %g(:,i) = reshape(G,1,16);
        g(:,i) = flatten(G);
    end
    E = 1e6;
    G = E/3;
    L = 10e-2;
    D = 1e-2;
    rho = 1000;
    A = pi/4*D^2;
    I = pi/64*D^4;
    J = 2*I;
    K = diag([E*I,E*I,G*J,G*A,G*A,E*A]);
    M = rho*diag([I,I,J,A,A,A]);
    xi_ref = [0;0;0;0;0;1];
    mu = M*eta;
    lambda = K*(xi-xi_ref);
end

function out = flatten(g)
    out = [reshape(g(1:3,1:3),9,1);g(1:3,4)];
end

function out = unflatten(g)
    out = [reshape(g(1:9),3,3),g(10:12);0,0,0,1];
end
