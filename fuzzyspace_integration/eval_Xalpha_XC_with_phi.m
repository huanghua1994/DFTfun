function XC = eval_Xalpha_XC_with_phi(natom, nbf, phi, ipw, D)
% Evaluate XC matrix using Xalpha functional and phi
% TODO: find references of formulas used in this function
% Input parameters:
%   nbf   : Total number of basis functions
%   phi   : Values of basis function at integral points
%   ipw   : Integral point weights
%   D     : Density matrix
% Output parameter:
%   XC : Exchange-correlation matrix

    nintp = size(phi, 1);
    XC = zeros(nbf, nbf);
    
    % rho_g = \sum_{u, v} phi_{g, u} * D_{u, v} * D_{g, v} is the electron density 
    rho = zeros(nintp, 1);
    for i = 1 : nbf
        rho_j = zeros(nintp, 1);
        for j = 1 : nbf
            rho_j = rho_j + D(i, j) * phi(:, j);
        end
        rho = rho + phi(:, i) .* rho_j;
    end
    
    % Slater exchange kernel: V_x[rho] = -3/4 (3/pi)^(1/3) rho^(1/3)
    % Xalpha kernel = 3/2 * alpha * V_x[rho], alpha = 0.7 includes correlation term
    % But don't know why we need to use (6/pi)^(1/3) instead of (3/pi)^(1/3)
    ker = -0.7 * (9/8) * ((6/pi)^(1/3)) * rho.^(1/3); 
    
    % ???
    ker = ker .* ipw;
    for i = 1 : nbf
        for j = i : nbf
            XC(i, j) = sum(phi(:, i) .* phi(:, j) .* ker);
            XC(j, i) = XC(i, j);
        end
    end
end