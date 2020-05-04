function [Cr_te] = Cte_relative(Swn, m_te)
    % Function to calculate relative TE coupling as a function of Swn
    % where:    Cr      = relative coupling coefficient [-]
    %           Swn     = normalised water saturation [-]
    %           m       = exponent and/or coefficients of Swn in Cr
    

    Cr_te=(1-Swn).^m_te;

end

