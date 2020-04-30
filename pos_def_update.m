function [M0Update, MMatsUpdate] = pos_def_update(M0, MMats)
    % Checks a matrices for positive-definitness. If matrix is not
    % positive-definite, finds nearest positive-definite matrix and
    % replaces
    %
    % Parameters
    % ----------
    % M0 : 2D matrix
    % MMats : 3D matrix
    %
    % Returns
    % -------
    % M0Update : 2D Pos Def Matrix
    % MMatsUpdate : 3D Pos Def Matrix
    
    [r,c,amountCoeff] = size(MMats);
    
    if is_pos_def(M0)
        M0Update = M0;
    else
        M0Update = nearest_SPD(M0);
    end
    
    MMatsUpdate = zeros(r,c,amountCoeff);
    
    for i = 1:amountCoeff
        if is_pos_def(MMats(:,:,i))
            MMatsUpdate(:,:,i) = MMats(:,:,i);
        else
            MMatsUpdate(:,:,i) = nearest_SPD(MMats(:,:,i));
        end
    
end % of the function