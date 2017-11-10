function volumeIn = flowRateVolumetric( maxVelocity )
%flowRateVolumetric - Description
%
% Syntax: volumeIn = flowRateVolumetric( maxVelocity )
% maxVelocity = the max volocity of a Poiseuill flow in a rectangular cross section tube
% volumeIn    = the volumetric flow rate of this Poiseuill flow
%
% Long description
% the Poiseuill flow can be described in three equivalent coefficient:
% 1. Max velocity               : the velocity in the center of flow's profile/intersection/corss section
% 3. volumetric flow rate       : the fluid's volume enter in a tube per unit time
% 2. averange velocity          : the volumetric flow rate divide the area of intersection of inlet/ intersection
% This function aimes at calculating the flow rate basing on the known max velocity
% *****************************************************************************************************************
        % Setup Build-in variabels  
        Width     = 50e-6;    % [m]
        Height    = 50e-6;    % [m]
        b         = Width / 2;
        c         = Height / 2;
        viscosity = 1.002e-3;
        n_order   = 100;

        % first step: calculate the pression derivative on main flow direction (x direction)
        y_0  = 0;
        z_0  = 0;
        first_term = 1 - (z_0/c)^2 ;
        series_term = 0;
        for k = 1:1:n_order
                series_term = series_term + ...,
                               4*( (-1)^k / aK(k)^3 ) * ...,
                               ( cosh(aK(k)*y_0/c) / cosh(aK(k)*b/c) ) * ...,
                               cos(aK(k)*z_0/c);
        end
        % Px = d(Pression)/dx
        Px = maxVelocity / ( -1/(2*viscosity)*(c^2)*(first_term + series_term) );

        %calculate the flow rate
        series_term = 0;
        for k = 1:1:n_order
                series_term = series_term + ( tanh(aK(k)*b/c) / (aK(k)^5) );
        end
        Q = -4/(3*viscosity)*Px*b*(c^3)*( 1 - 6*c/b*series_term );
        volumeIn = Q;
end