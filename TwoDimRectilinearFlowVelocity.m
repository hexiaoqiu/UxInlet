function u_x = TwoDimRectilinearFlowVelocity(y, z, U_max, Width, Height, n_order)
%myFun - Description
%
% Syntax: dimensional velocity on inlet = VelocityInletSquare(input)
% Inout arguments: z        ->  z coordinate dimensional height
%                  y        ->  y coordinate dimensional width
%                  U_max    ->  the center line velicity
%                  Width    ->  the side length of the square channel
%                  Height   ->  the radius of the particle
%                  n_order  ->  the order of the precison of the formula (by default is 10)
%                  
% the function is meant to calculate one point's velocity of Steady Two-Dimensional Rectilinear flow
% Steady Two-Dimensional Rectilinear flow = 2-Dimensional Poiseuill Flow 
%                                         = fully developed laminar flow in a tube

% first step: check the input arguments if they are licit  
        if (U_max < 0)||(Width < 0)||(Height < 0)
                u_x = -1;
                fprintf('The input argument is invalid! U_max or Width or Height is smaller than 0 \n')
                return
        end

        if abs(y) > (Width/2) 
                u_x = -1;
                fprintf('The input argument is invalid! y coordinate is out of channel! \n')
                return   
        end        
        
        if abs(z) > (Height/2) 
                u_x = -1;
                fprintf('The input argument is invalid! y coordinate is out of channel! \n')                
                return   
        end

        if (abs(z) == Height/2)||(abs(y) == Width/2)
                u_x = 0;
                fprintf('y or z coordinates are on the wall! \n')
                return 
        end

        if n_order < 10
                n_order = 10;
        end
        
% second step: calculating the undimensional cooridnate (the origin is at the center of tube's intersection)
        b = Width / 2;
        c = Height / 2;
        viscosity = 1.002e-3;
        
% third step: calculating the Px (pression x direction derivative)
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

        Px = U_max / ( -1/(2*viscosity)*(c^2)*(first_term + series_term) );

% fourth step: calculating the series term in parentheses
        first_term = 1 - (z/c)^2 ;
        series_term = 0;
        for k = 1:1:n_order
                series_term = series_term + ...,
                               4*( (-1)^k / aK(k)^3 ) * ...,
                               ( cosh(aK(k)*y/c) / cosh(aK(k)*b/c) ) * ...,
                                cos(aK(k)*z/c);
        end
        u_x = -1/(2*viscosity) * Px * (c^2) * (first_term + series_term);
        
end