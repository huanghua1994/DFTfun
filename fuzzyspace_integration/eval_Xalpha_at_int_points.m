function rhos = eval_Xalpha_at_int_points(int_points, atom_xyz, coef, alpha, shapematrix, centers, K, L, P)
	[ncenters, ~] = size(atom_xyz);

	for iatom = 1 : ncenters
		% 将中心从原点移动到原子中心
		points(:, 1) = int_points(:, 1) + atom_xyz(iatom, 1);
		points(:, 2) = int_points(:, 2) + atom_xyz(iatom, 2);
		points(:, 3) = int_points(:, 3) + atom_xyz(iatom, 3);
		
		rho0 = get_densitynck(points(:,1),points(:,2),points(:,3),coef,alpha,shapematrix,centers,K,L,P);
		rho0 = rho0.^(1/3);
	   
		for mu = 1 : K
			for nu = 1 : K
				rhos{iatom,mu,nu} = rho0 .* get_density(mu,nu,points(:,1),points(:,2),points(:,3),coef,alpha,shapematrix,centers,K,L,P);
			end
		end 
	end
end