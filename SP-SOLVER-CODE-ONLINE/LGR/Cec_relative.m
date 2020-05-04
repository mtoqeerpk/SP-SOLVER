function [Cr_ec] = Cec_relative(Swn, m_ec)
    % Function to calculate relative EC coupling as a function of Swn
    % where:    Cr      = relative coupling coefficient [-]
    %           Swn     = normalised water saturation [-]
    %           m       = exponent and/or coefficients of Swn in Cr
    

    Cr_ec=(1-Swn).^m_ec;
  
end

