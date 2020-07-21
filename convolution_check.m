function convolution = convolution_check(C0, CMats, rhos, S0, SMats, lambdas, t)
    % Checks the convolution integral of the C(t) and S(t) matrices.
    % 
    % Performs the numerical convolution integration of the creep and
    % relaxation material functions, using the integral function. The returned
    % matrix of the time convolution is an array of the elements of the columns
    % are added together. While, if computed by hand, would be an idenity
    % matrix multiplied by the given time.
    %
    % Parameters
    % ----------
    % C0 : 2D matrix
    %     The equilibrium relaxation in a 2D array.
    % CMats : 3D matrix
    %     The relaxation modulus coefficient matrices in a 3D array. The third
    %     dimension is to access the matrix, while the first and second are the
    %     rows and columns.
    % rhos : 1D array
    %     The inverted relaxation time constants in a 1D array, in descending
    %     order.
    % S0 : 2D matrix
    %     The instantaneous creep modulus in a 2D array.
    % SMats : 3D matrix
    %     The creep modulus coefficient matrices in a 3D array. The third
    %     dimension is to access the matrix, while the first and second are the
    %     rows and columns.
    % lambdas : 1D array
    %     The inverted creep time constants in a 1D numpy array, in descending
    %     order.
    % t : float
    %     The time at which to calculate the convolution.
    % 
    % Returns
    % -------
    % convolution : 1D array
    %     An array of the convolution matrix for C(t) and S(t).
    
    C = @(time) relax_modulus(C0, CMats, rhos, time);
    S = @(time) creep_modulus(S0, SMats, lambdas, time);
    
    f  = @(tau, t) dot(C(t - tau), S(tau));
    
    convolution = integral(@(tau) f(tau, t), 0, t, 'ArrayValued', true);
    
    end % of the function