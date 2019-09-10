function T = timingDynamics(N)
    %given the values of N to test try getting the average time for a
    %single step of the dynamics
    
    steps = 1000;
    T = zeros(1,length(N));
    
    for i=1:length(N)
        n = N(i)
        [g,xi,eta] = initDynamics(n);
        tic;
        for j=1:steps
            %[g,xi,eta] = implicit_dynamics(g,xi,eta,[-0;0;0],0.1);
            [g_s,xi_s,eta_s,b_s] = implicit_dynamics_step(g,xi,eta,[0;0;0],xi(:,1),0.1);
        end
        T(i) = toc/steps;
    end
end