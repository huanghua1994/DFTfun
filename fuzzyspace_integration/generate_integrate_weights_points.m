function [int_points, int_weights] = generate_integrate_weights_points(atom_xyz)
% atom_xyz    : 各原子中心坐标
% int_points  : 各积分格点坐标
% int_weights : 各积分格点权重

    [ncenters, ~] = size(atom_xyz);
    
    % 原作者设的参数值
    P        = 5;    % 原作者设的 P 值，不知为何是 5
    nCheby   = 75;   % 径向使用 75 个节点
    nLebedev = 302;  % 球面使用 302 个 Lebedev 节点    

    % 生成 Lebedev 节点
    Lebedev_points = getLebedevSphere(nLebedev);
    
    % 生成径向与球面积分节点
    [rad_r, rad_w] = my_cheb2_becke(nCheby, P);
    
    % 将径向节点与球面节点进行组合
    rad_cut_num  = sum(rad_r < 10);
    int_weights0 = zeros(rad_cut_num * nLebedev, 1);
    int_points   = zeros(rad_cut_num * nLebedev, 3);
    orig_points  = [Lebedev_points.x, Lebedev_points.y, Lebedev_points.z];
    for ir = 1 : rad_cut_num
        spos = ir * nLebedev - nLebedev + 1;
        epos = ir * nLebedev;
        int_weights0(spos : epos) = Lebedev_points.w * rad_w(ir);
        int_points(spos : epos, 1 : 3) = orig_points * rad_r(ir);
    end
    
    % W_mat(i, j, k) 表示当前的第 i 个积分节点关于原子对 j, k 的模糊划分归属
    W_mat = zeros(rad_cut_num * nLebedev, ncenters, ncenters);
    
    dist = squareform(pdist(atom_xyz));
    for iatom = 1 : ncenters
        % 将积分格点中心从原点移动到原子中心
        rnowx = int_points(:, 1) + atom_xyz(iatom, 1); 
        rnowy = int_points(:, 2) + atom_xyz(iatom, 2); 
        rnowz = int_points(:, 3) + atom_xyz(iatom, 3); 
        
        for jatom = 1 : ncenters %遍历每个原子对 jk，计算空间划分权重 W_i
            dxj = rnowx - atom_xyz(jatom, 1);
            dyj = rnowy - atom_xyz(jatom, 2);
            dzj = rnowz - atom_xyz(jatom, 3);
            dij = sqrt(dxj.^2 + dyj.^2 + dzj.^2);
            for katom = 1 : ncenters
                if katom ~= jatom
                    dxk = rnowx - atom_xyz(katom, 1);
                    dyk = rnowy - atom_xyz(katom, 2);
                    dzk = rnowz - atom_xyz(katom, 3);
                    dik = sqrt(dxk.^2 + dyk.^2 + dzk.^2);
                    
                    mu  = (dij - dik) / dist(jatom, katom);
                    
                    % 迭代三次，得到 s(v(i,j)) = 0.5 * (1 - p(p(p(v(i,j)))))
                    smu = 1.5 * mu  - 0.5 *  mu.^3;
                    smu = 1.5 * smu - 0.5 * smu.^3;
                    smu = 1.5 * smu - 0.5 * smu.^3;

                    W_mat(:, jatom, katom) = (0.5 * (1 - smu)); %储存遮蔽矩阵jk
                else
                    W_mat(:, jatom, katom) = 1;
                end
           end
        end

        % 对每个积分格点 i，遮蔽矩阵 W_mat(i, j, :) 的乘积即为格点 i 属于原子 j 的权重
        pvec = ones(radius_cut_num * nLebedev, ncenters);
        for iw = 1 : ncenters  
            for iz = 1 : ncenters     
                pvec(:, iw) = pvec(:, iw) .* W_mat(:, iw, iz);
            end  
        end

        % 对 pvec 求和以归一化
        sum_pvec = zeros(radius_cut_num * nLebedev, 1);
        for ix = 1 : ncenters
            sum_pvec = pvec(:, ix) + sum_pvec;
        end
        
        % 得到以当前原子（iatom）为中心的积分格点属于当前原子的权重
        int_weights{iatom} = int_weights0 .* pvec(:, iatom) ./ sum_pvec;
    end

end % end function    
 