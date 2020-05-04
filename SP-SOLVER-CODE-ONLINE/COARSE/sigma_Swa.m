function [sigma] = sigma_Swa(sigma1, S, PORO3D)

    % Archie function for electrical conductivity
   
    
    sigma=sigma1.*(PORO3D.^1.8).*(S.^2);

    
end