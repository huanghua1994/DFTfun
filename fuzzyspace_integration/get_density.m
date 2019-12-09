function psi = get_density(mu, nu, x, y, z, coef, alpha, shapematrix, centers, K, L, P)
%  specific orbital's density
    psi=0;         

    shapeA=shapematrix(mu,:);
    shapeB=shapematrix(nu,:);
    xlA=shapeA(1);
    ymA=shapeA(2);
    znA=shapeA(3);
    
    xlB=shapeB(1);
    ymB=shapeB(2);
    znB=shapeB(3);

    xA = centers(mu,1);
    yA = centers(mu,2);
    zA = centers(mu,3);
    xB = centers(nu,1);
    yB = centers(nu,2);
    zB = centers(nu,3);
    for p=1:L(mu) % 遍历两个一个STO对应的高斯函数
        for q=1:L(nu) % 遍历两个一个STO对应的高斯函数             
            alpha1 = alpha(mu,p);
            alpha2 = alpha(nu,q);
           
            Psi1=(coef(mu,p).*(x-xA).^xlA.*(y-yA).^ymA.*(z-zA).^znA).*(exp(-alpha1*(x-xA).^2).*exp(-alpha1*(y-yA).^2).*exp(-alpha1*(z-zA).^2));
            Psi2=(coef(nu,q).*(x-xB).^xlB.*(y-yB).^ymB.*(z-zB).^znB).*(exp(-alpha2*(x-xB).^2).*exp(-alpha2*(y-yB).^2).*exp(-alpha2*(z-zB).^2));     
  
            psi=psi+Psi1.*Psi2;
        end
    end
end