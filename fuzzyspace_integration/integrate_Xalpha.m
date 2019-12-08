function [weightss,rhos] = integrate_Xalpha(xyz,d,spreads,shapematrix,centers,K,L,P)

[ncenters,~]=size(xyz);
    
for iatom=1:ncenters
	for jatom=1:ncenters
		if iatom~=jatom
			dist(iatom,jatom)=2.07;
		end
	end
end

% 径向全部使用 75 个节点
[r0,w1] = cheby2(75,[-1,1]); % r0 == x_i, Cheb 节点, w1 是权重

% 作者在这里将 P 取为了 5
r1= 5* (1+r0)./(1-r0);   % R_i = P \frac{1 + x_i}{1 - x_i}
w1= w1'.*5*2./(r0-1).^2; % dR = \frac{2P}{(1-x)^2} dx

r1=r1(1:end-1);
w1=w1(1:end-1);

% 球面全部使用 302 个节点
surface_points=getLebedevSphere(302);

weights=[0];
points0=[0,0,0];
originalpoints=[surface_points.x,surface_points.y,surface_points.z];
% 将 Lebedev 节点与 R_i 节点组合，如果 R_i 大于 10 bohr 则忽略之
for ix=1:length(r1)
	rnow=r1(ix);  
	if rnow<10 
		points0=[points0;originalpoints*rnow];
		weights=[weights;w1(ix)*(surface_points.w)*rnow^2];
	end
end
points0=points0(2:end,:);
weights=weights(2:end);
    
spaceweight=cell(ncenters,ncenters);
totalresult=0;

for iatom=1:ncenters %遍历各个原子
    
	% 将中心从原点移动到原子中心
	for abc=1:3%xyz
		points(:,abc)=points0(:,abc)+xyz(iatom,abc);
	end
    rnowx=points(:,1);
    rnowy=points(:,2);
    rnowz=points(:,3);
    
      
    for jatom=1:ncenters %遍历每个原子对 jk，计算空间划分权重 W_i
        ri=sqrt((rnowx-xyz(jatom,1)).^2+(rnowy-xyz(jatom,2)).^2+(rnowz-xyz(jatom,3)).^2);
        iz=1;
        for katom=1:ncenters
            if katom~=jatom
                rj=sqrt((rnowx-xyz(katom,1)).^2+(rnowy-xyz(katom,2)).^2+(rnowz-xyz(katom,3)).^2);      
				rmu=(ri-rj)/dist(jatom,katom); %判断积分格点是否在 j k 中间线的位置 
				
				% 迭代三次，得到 s(v(i,j)) = 0.5 * (1 - p(p(p(v(i,j)))))
				partfun = @(r) 1.5*r-0.5*r.^3;
				tmp1=partfun(partfun(partfun(rmu)));
				spaceweight{jatom,katom}=(0.5*(1-tmp1)); %储存遮蔽矩阵jk
				iz=iz+1;
            else
                spaceweight{jatom,katom}=1;
            end
       end
    end

    for ix=1:ncenters
        pvec{ix}=1;
	end

    for iw=1:ncenters  
        for iz=1:ncenters     
            tmp=spaceweight{iw,iz};
            pvec{iw}=pvec{iw}.*tmp; %对每个积分格点，遮蔽矩阵的乘积=权重
        end  
    end
   
    oneatomresult=0;
	
	% 对 pvec 求和以归一化
    sum_pvec=0;
    for ix=1:ncenters
        sum_pvec=(pvec{ix}+sum_pvec);
    end
	weightss{iatom}=weights.*pvec{iatom}./sum_pvec;
    
	oneatomresult= oneatomresult+sum(weights.*pvec{iatom}./sum_pvec.*get_densitynck(points(:,1),points(:,2),points(:,3),d,spreads,shapematrix,centers,K,L,P));

    rho0= get_densitynck(points(:,1),points(:,2),points(:,3),d,spreads,shapematrix,centers,K,L,P);
   
	for mu=1:K
		for nu=1:K
			rhos{iatom,mu,nu}= (rho0.^(1/3)).*get_density(mu,nu,points(:,1),points(:,2),points(:,3),d,spreads,shapematrix,centers,K,L,P);
		end
	end 
	totalresult=totalresult+ oneatomresult;
end

end % end function    
 