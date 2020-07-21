function modTime = modulus_at_time(M0, MMats, timeConst, time, property)
    % Computes the anisotropic modulus, at a given time.
    %
    % A single function that calculates the modulus of an anisotropic material
    % function at a given time. Requires the instantaneous/equilibrium matrix, 
    % coefficient matrices, the inverted time constants, and the time at which the
    % modulus will be calculated. A property flag is also passed to the function
    % to switch between creep or relaxation calculations.
    %
    % Parameters
    % ----------
    % M0 : 2D matrix
    %     An instantaneous or equilibrium modulus in a 2D array.
    % MMats : 3D matrix
    %     The coefficient moduli in a 3D array, with the third dimension 
    %     consisting of the matrices and the first and second being the row and
    %     columns of the matrix.
    % timeConst : vector
    %     The corresponding inverted time constants in a 1D array, in
    %     descending order.
    % time : float
    %     The time at which to calculate the modulus.
    % property : string
    %     Flag to set the function to calculate the modulus of a creep material
    %     function or of a relaxation material function.
    %     Valid inputs are: 'creep' or 'relax'
    %
    % Returns
    % -------
    % modTime : 2D matrix
    %     The modulus, at the given time, for the prescribed material function.
    %     Returned in a 2D array.
    
    % Find the size of the square matrix and the number of coefficent matrices
    [~, coeffSize, numCoeff] = size(MMats);
    
    % Preallocate the modulus matrix
    modTime = zeros(coeffSize);
    
    % Calculate the relaxation modulus
    if property == 'relax'
        
        exponential = exp(-time * timeConst);
        
        for i = 1:numCoeff
            modTime = modTime + (MMats(:,:,i) * exponential(i));
        end
        
        modTime = modTime + M0;
        
    end
    
    % Calculate the creep modulus
    if property == 'creep'
        
        exponential = exp(-time * timeConst);
        
        for i = 1:numCoeff
            modTime = modTime + (MMats(:,:,i) * (1 - exponential(i)));
        end
        
        modTime = modTime + M0;
        
    end

    end % of the function