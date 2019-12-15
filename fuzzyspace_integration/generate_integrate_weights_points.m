function [ip, ipw] = generate_integrate_weights_points(atom_xyz)
% Get numerical integral points and their fuzzy weights for XC integral
% Ref: [JCP 88, 2547], doi: 10.1063/1.454033, and http://sobereva.com/69
% Input parameter:
%   atom_xyz : Atom coordinates
% Output parameters:
%   ip  : Numerical integral points of all atoms
%   ipw : Fuzzy weights of numerical integral points

    % The values of nCheby and nLebedev are chosen by the original author
    rm       = 1;    % A parameter, see my_cheb2_becke
    nCheby   = 75;   % Number of radial direction integral points + 1
    nLebedev = 302;  % Number of points on a sphere

    % TODO: prune grid points for different atoms, so the number of integral points
    % per atom could be different and the total number of grid points can be reduced

    % Generate Lebedev points on a sphere
    Lebedev_points = getLebedevSphere(nLebedev);
    
    % Generate radial direction points & weights
    [rad_r, rad_w] = my_cheb2_becke(nCheby, rm);
    
    % Combine radial direction points with sphere points
    n_radcut = sum(rad_r < 10);
    ip_atom  = zeros(n_radcut * nLebedev, 3);
    ipw_atom = zeros(n_radcut * nLebedev, 1);
    orig_ip  = [Lebedev_points.x, Lebedev_points.y, Lebedev_points.z];
    for ir = 1 : n_radcut
        spos = (ir-1) * nLebedev + 1;
        epos = ir * nLebedev;
        ip_atom(spos : epos, 1 : 3) = orig_ip * rad_r(ir);
        ipw_atom(spos : epos) = Lebedev_points.w * rad_w(ir);
    end
    nintp_atom = n_radcut * nLebedev;  % Total number of integral points
    natom = size(atom_xyz, 1);
    ip    = zeros(nintp_atom * natom, 3);
    ipw   = zeros(nintp_atom * natom, 1);
    
    % W_mat(i, j, k): fuzzy weight of integral point i to atom pair (j, k)
    % W_mat(i, j, k): 第 i 个积分节点关于原子对 (j, k) 的模糊划分权重
    W_mat = zeros(nintp_atom, natom, natom);
    dist  = squareform(pdist(atom_xyz));
    for iatom = 1 : natom
        % Shift the integral point to atom
        rnowx = ip_atom(:, 1) + atom_xyz(iatom, 1); 
        rnowy = ip_atom(:, 2) + atom_xyz(iatom, 2); 
        rnowz = ip_atom(:, 3) + atom_xyz(iatom, 3); 
        
        % Enumerate each atom pair (j, k)
        for jatom = 1 : natom 
            dxj = rnowx - atom_xyz(jatom, 1);
            dyj = rnowy - atom_xyz(jatom, 2);
            dzj = rnowz - atom_xyz(jatom, 3);
            dij = sqrt(dxj.^2 + dyj.^2 + dzj.^2);
            for katom = 1 : natom
                if katom ~= jatom
                    dxk = rnowx - atom_xyz(katom, 1);
                    dyk = rnowy - atom_xyz(katom, 2);
                    dzk = rnowz - atom_xyz(katom, 3);
                    dik = sqrt(dxk.^2 + dyk.^2 + dzk.^2);
                    smu = (dij - dik) / dist(jatom, katom);
                    
                    % s(d(i,j)) = 0.5 * (1 - p(p(p(d(i,j)))))
                    for k = 1 : 3
                        smu = 1.5 * smu - 0.5 * smu.^3;
                    end
                    W_mat(:, jatom, katom) = (0.5 * (1 - smu));
                else
                    W_mat(:, jatom, katom) = 1;
                end
           end
        end
        
        % \prod_{k} W_mat(:, j, k) is the actual weight of integral points
        % belonging to atom k. Normalizing it gives us the fuzzy weight. 
        % 积分节点上第 j 个原子的实际权重等于在此处原子 j 与所有其他原子
        % k 之前权重的乘积。归一化后得到以当前原子中心的积分格点属于当前
        % 原子的模糊划分权重。
        pvec = ones(nintp_atom, natom);
        for i = 1 : natom  
            for j = 1 : natom     
                pvec(:, i) = pvec(:, i) .* W_mat(:, i, j);
            end  
        end
        sum_pvec = sum(pvec, 2);
        
        sidx = (iatom-1) * nintp_atom + 1;
        eidx = iatom * nintp_atom;
        ip(sidx : eidx, 1) = rnowx;
        ip(sidx : eidx, 2) = rnowy;
        ip(sidx : eidx, 3) = rnowz;
        ipw(sidx : eidx) = ipw_atom .* pvec(:, iatom) ./ sum_pvec;
        
    end  % end of iatom loop
end  % end function    