function convolution = convolution_check(C0, CMats, rhos, S0, SMats, lambdas, t)
    % Checks the convolution intergral of the C(t) and S(t) matrices.
    % Should be the identiy matrix of t*I
    % Also returns the errors from scipy.quad intergration
    %
    % Parameters
    % ----------
    % C0 : 2D matrix
    %     Equilibrium relaxation
    % CMats : 3D matrix
    %     Relaxtion modulus coefficient
    % rhos : 1D array
    %     Relaxation time constants
    % S0 : 2D matrix
    %     Instantaneous creep modulus
    % SMats : 3D matrix
    %     Creep modulus coefficient
    % lambdas : 1D array
    %     Creep time constants
    % t : float
    %     time
    % 
    % Returns
    % -------
    % convolution : 2D matrix
    %     Matrix of the convolution of C(t) and S(t)
    
    C = @(time) modulus_at_time(C0, CMats, rhos, time, 'relax');
    S = @(time) modulus_at_time(S0, SMats, lambdas, time, 'creep');
    
    f  = @(tau, t) dot(C(t - tau), S(tau));
    
    convolution = integral(@(tau) f(tau, t), 0, t, 'ArrayValued', true);
    % Statements here; these must include putting a value in the output argument
    end % of the function