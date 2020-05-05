function modTime = modulus_at_time(M0, MMats, timeConst, time, property)
    % Gives the matrix of creep or relaxation at a given time
    %
    % Parameters
    % ----------
    % M0 : 2D matrix
    %     Instantaneous/Equilibrium modulus
    % MMats : 3D matrix
    %     Coefficient moduli
    % timeConst : vector or inverted rhos or lambdas
    %     Time constants
    % time : float
    %     Time at which to calculate property
    % property : string
    %     Switch for creep or relaxation calculations (relax or creep)
    %
    % Returns
    % -------
    % modTime : 2D matrix
    %     Modulus at a given time
    
    [~, coeffSize, numCoeff] = size(MMats);
    
    modTime = zeros(coeffSize);
    
    if property == 'relax'
        
        disp(time)
        %disp(timeConst)
        exponential = exp(-time * timeConst);
        
        for i = 1:numCoeff
            modTime = modTime + (MMats(:,:,i) * exponential(i));
        end
        
        modTime = modTime + M0;
        
    end
    
    if property == 'creep'
        
        exponential = exp(-time * timeConst);
        
        for i = 1:numCoeff
            modTime = modTime + (MMats(:,:,i) * (1 - exponential(i)));
        end
        
        modTime = modTime + M0;
        
    end

    end % of the function