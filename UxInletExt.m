function u_X = UxInletExt(y,z,maxVelocity,yLength,zLength)
% UxInlet - Description
% This func aimed at offer a precise boundary condition for COMSOL's CFD module
% In order to comply the convention of COMSOL's matlab extern function, the input and output
% are vectors with the same size
% the function is meant to calculate one point's velocity of Steady Two-Dimensional Rectilinear flow
% Steady Two-Dimensional Rectilinear flow = 2-Dimensional Poiseuill Flow 
%                                         = fully developed laminar flow in a tube
%
% Syntax: dimensional velocity on inlet = VelocityInletSquare(input)
% Input arguments: 
%                  z        ->  z coordinate dimensional height
%                  y        ->  y coordinate dimensional width
% Output arguments:
%                  u_X      -> the main flow direction's velocity in the position y z
% Build in Variables:
%                  U_max    ->  the center line velicity
%                  Width    ->  the side length of the square channel
%                  Height   ->  the radius of the particle
%                  n_order  ->  the order of the precison of the formula (by default is 10)
%                  


% the external function offered by Matlab and invoked by COMSOL must take vectors as input arguments 
        if size(y) ~= size(z)
                fprintf('the size of input is not consistent!')
                return;
        end
        % get the size of input vectors
        outputSize = size(y);
        numLine    = outputSize(1);
        numColomn  = outputSize(2); % the number of colomn is always 1

        nLoop = numLine;
        u_X = zeros(outputSize);

% Setup Build-in variabels  
        U_max   = maxVelocity(1);      % [m/s]
        Width   = yLength(1);    % [m]
        Height  = zLength(1);    % [m]
        n_order = 100;
% call func to calculate the velocity for each position
        for index = 1:1:nLoop
                y_i = y(index);                
                z_i = z(index);
                % y_i = y(index) - Width/2;                
                % z_i = z(index) - Height/2;
                u_X(index) = TwoDimRectilinearFlowVelocity(y_i, z_i, U_max, Width, Height, n_order);
        end

end