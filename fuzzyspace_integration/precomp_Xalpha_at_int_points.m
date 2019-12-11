function rhos = precomp_Xalpha_at_int_points(int_points, atom_xyz, nbf, bf_coef, bf_alpha, bf_exp, bf_center, bf_nprim)
    [ncenters, ~] = size(atom_xyz);

    rhos = cell(ncenters, nbf, nbf);
    for iatom = 1 : ncenters
        % Shift the integral point center to atom center
        x = int_points(:, 1) + atom_xyz(iatom, 1);
        y = int_points(:, 2) + atom_xyz(iatom, 2);
        z = int_points(:, 3) + atom_xyz(iatom, 3);
        
        rho0 = 0;
        for mu = 1 : nbf
        for nu = 1 : nbf
            % Exponents of each basis function pair
            xlA = bf_exp(mu, 1);
            ymA = bf_exp(mu, 2);
            znA = bf_exp(mu, 3);
            xlB = bf_exp(nu, 1);
            ymB = bf_exp(nu, 2);
            znB = bf_exp(nu, 3);
            
            % Center of each basis function pair
            xA = bf_center(mu, 1);
            yA = bf_center(mu, 2);
            zA = bf_center(mu, 3);
            xB = bf_center(nu, 1);
            yB = bf_center(nu, 2);
            zB = bf_center(nu, 3);
            
            poly_coef1 = (x-xA).^xlA .* (y-yA).^ymA .* (z-zA).^znA;
            poly_coef2 = (x-xB).^xlB .* (y-yB).^ymB .* (z-zB).^znB;
            
            % Psi1: mu-th basis function's value at location (x, y, z)
            % Psi2: nu-th basis function's value at location (x, y, z)
            Psi1 = zeros(length(x), bf_nprim(mu));
            Psi2 = zeros(length(x), bf_nprim(nu));
            for p = 1 : bf_nprim(mu)
                alpha1 = bf_alpha(mu, p);
                Psi1(:, p) = (bf_coef(mu,p) .* poly_coef1) .* exp(-alpha1*(x-xA).^2) .* exp(-alpha1*(y-yA).^2) .* exp(-alpha1*(z-zA).^2);
            end
            for q = 1 : bf_nprim(nu)
                alpha2 = bf_alpha(nu, q);
                Psi2(:, q) = (bf_coef(nu,q) .* poly_coef2) .* exp(-alpha2*(x-xB).^2) .* exp(-alpha2*(y-yB).^2) .* exp(-alpha2*(z-zB).^2);
            end
            
            psi = 0;
            for p = 1 : bf_nprim(mu)
                for q = 1 : bf_nprim(nu)
                    Psi12 = Psi1(:, p) .* Psi2(:, q);
                    psi   = psi  + Psi12;
                end
            end
            
            rhos{iatom,mu,nu} = psi;
        end
        end
    end
end