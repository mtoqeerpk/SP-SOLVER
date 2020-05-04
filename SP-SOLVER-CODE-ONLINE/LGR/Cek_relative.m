function [Cr_ek] = Cek_relative(Swn, m_ek)
    % Function to calculate relative EK coupling as a function of Swn
    % where:    Cr      = relative coupling coefficient [-]
    %           Swn     = normalised water saturation [-]
    %           m       = exponent and/or coefficients of Swn in Cr
    
    Cr_ek=(Swn).^m_ek;
    
end