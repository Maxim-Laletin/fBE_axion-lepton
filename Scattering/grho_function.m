function g_rho = grho_function(T)
% from 1803.01038
    if T >= 120E-3 % GeV

        % Coefficients from Table 1
        a = [1, 1.11724, 0.312672, -0.0468049, -0.0265004, -0.0011976, ...
             0.000182812, 0.000136436, 8.55051e-05, 1.2284e-05, 3.82259e-07, -6.87035e-09];
        b = [0.0143382, 0.0137559, 0.00292108, -0.000538533, -0.000162496, -2.87906e-05, ...
             -3.84278e-06, 2.78776e-06, 7.40342e-07, 1.1721e-07, 3.72499e-09, -6.74107e-11];
        
        % Compute log temperature
        t = log(T);
        
        % Compute polynomial sums
        num_g_rho = sum(a .* t.^(0:11));
        den_g_rho = sum(b .* t.^(0:11));
        
        
        % Compute g_rho and g_s
        g_rho = num_g_rho ./ den_g_rho;

    else
    
        % Parameter values in GeV
        me = 511e-6;
        mmu = 0.1056;
        mpi0 = 0.135;
        mpipm = 0.140;
        m1 = 0.5;
        m2 = 0.77;
        m3 = 1.2;
        m4 = 2.0;
        
        % Define helper functions
        S_fit = @(x) 1 + (7/4) * exp(-1.0419*x) .* (1 + 1.034*x + 0.456426*x.^2 + 0.0595249*x.^3);
        f_rho = @(x) exp(-1.04855*x) .* (1 + 1.03757*x + 0.508630*x.^2 + 0.0893988*x.^3);
        b_rho = @(x) exp(-1.03149*x) .* (1 + 1.03317*x + 0.398264*x.^2 + 0.0648056*x.^3);
        
        % Compute g_rho(T) using Eq. (C.3)
        g_rho = 2.030 + 1.353 * S_fit(me./T).^(4/3) + 3.495 * f_rho(me./T) + ...
                 3.446 * f_rho(mmu./T) + 1.05 * b_rho(mpi0./T) + 2.08 * b_rho(mpipm./T) + ...
                 4.165 * b_rho(m1./T) + 30.55 * b_rho(m2./T) + 89.4 * b_rho(m3./T) + 8209 * b_rho(m4./T);
    end 
end