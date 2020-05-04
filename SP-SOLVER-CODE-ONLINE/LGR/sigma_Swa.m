function [sigma] = sigma_Swa(sigma1, S, PORO3D)
    % Function to calculate electrical conductivity for each gridblock as a function of Sw
    % where:    sigma   = electrical conductivity [S/m]
    %           Sw      = water saturation [-]
    %           x       = exponent and/or coefficients of Sw in sigma
    
    
    % Archie 3 function for electrical conductivity
    %sigma=x(1).*Sw.^x(2).*x(3);
    
    sigma=sigma1.*(PORO3D.^1.8).*(S.^2);
     %sigma=sigma1.*(PORO3D.^0).*(S.^2);
    
%     % Archie 4 function for electrical conductivity
%     sigma=x(1).*Sw.^x(2)+Sw.^2*x(3)+Sw.*x(4)+x(5);
    
%     % Polynomial function for electrical conductivity Archie 2
%     sigma=(Sw.^3).*x(1)+(Sw.^2).*x(2)+Sw.*x(3)+x(4);

%     % Archie 1 function for electrical conductivity Archie 2
%     sigma=x(1).*Sw.^x(2);
    
end