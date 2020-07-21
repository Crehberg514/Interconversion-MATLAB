function modTime = creep_modulus(M0, MMats, timeConst, time)
    % Shortcut to compute the anisotropic creep modulus, at a given time.
    %
    % A shortcut function that presets the creep flag in the modulus_at_time
    % function. Calculates the creep modulus of an anisotropic material
    % function at a given time . Requires the instantaneous matrix, coefficient
    % matrices, the inverted time constants, and the time at which the modulus
    % will be calculated.
    %
    % Parameters
    % ----------
    % M0 : 2D matrix
    %     An instantaneous modulus in a 2D array.
    % MMats : 3D matrix
    %     The coefficient moduli in a 3D array, with the third dimension 
    %     consisting of the matrices and the first and second being the row and
    %     columns of the matrix.
    % timeConst : vector
    %     The corresponding inverted time constants in a 1D array, in
    %     descending order.
    % time : float
    %     The time at which to calculate the modulus.
    %
    % Returns
    % -------
    % modTime : 2D matrix
    %     The modulus, at the given time, for the prescribed material function.
    %     Returned in a 2D array.
    %     Uses the modulus_at_time function with the proper flag preset.
    
    modTime = modulus_at_time(M0, MMats, timeConst, time, 'creep');
    
    end % of the function