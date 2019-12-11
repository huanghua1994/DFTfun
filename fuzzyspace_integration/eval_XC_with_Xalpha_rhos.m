function XC = eval_XC_with_Xalpha_rhos(natom, nbf, rhos, int_weights, D)
    XC = zeros(nbf, nbf);
    for k = 1 : natom
    
        rho0 = 0;
        for i = 1 : nbf
        for j = 1 : nbf
            rho0 = rho0 + D(i, j) * rhos{k, i, j};
        end
        end
        rho0 = rho0.^(1/3);
        
        for i = 1 : nbf
        for j = i : nbf
            rhos_kij  = rhos{k, i, j} .* rho0;
            XC(i, j) = XC(i, j) + sum(int_weights{k} .* rhos_kij);
            XC(j, i) = XC(i, j);
        end
        end
    end
    XC = (-9/8) * ((3/pi)^(1/3)) * 0.7 * XC; 
end