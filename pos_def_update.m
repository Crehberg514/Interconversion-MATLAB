function [M0Update, MMatsUpdate] = pos_def_update(M0, MMats)
    % Checks and corrects a matrix for positive-definiteness.
    %
    % First, checks to see if a matrix has the property of
    % positive-definiteness. If the matrix has the property, then it returns
    % the matrix. If it does not, then the nearest semi-positive-definite form
    % of the matrix is found and returned.
    %
    % Parameters
    % ----------
    % M0 : 2D matrix
    %     The instantaneous or equilibrium matrix given in a 2D array
    % MMats : 3D matrix
    %     The coefficient moduli given in a 3D array, with the third dimension
    %     addressing the matrix, and the first and second dimensions addressing
    %     the rows and columns.
    %
    % Returns
    % -------
    % M0Update : 2D Matrix
    %     Returns the original matrix if it was already positive-definite.
    %     Or returns the nearest semi-positive-definite form of the
    %     instantaneous or equilibrium matrix in a 2D numpy array.
    % MMatsUpdate : 3D Matrix
    %     Returns a 3D array, with the same constraints as given to the
    %     function.
    %     Each array returned is the original matrix if it was already
    %     positive-definite.
    %     All other matrices in the array are the nearest
    %     semi-positive-definite form of the coefficient moduli
    
    [~, ~, amountCoeff] = size(MMats);

    if is_pos_def(M0)
        M0Update = M0;
    else
        M0Update = nearest_SPD(M0);
    end
    
    MMatsUpdate = zeros(size(MMats));
    
    for i = 1:amountCoeff
        if is_pos_def(MMats(:,:,i))
            MMatsUpdate(:,:,i) = MMats(:,:,i);
        else
            MMatsUpdate(:,:,i) = nearest_SPD(MMats(:,:,i));
        end
    end
    
end % of the function