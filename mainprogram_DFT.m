%******************************************************************
%  DFT demo code using Gaussian basis sets
%  
%  Copyright (C) 2014-2015 xiangrufan@GitHub <1034534198@qq.com>
%  Released in MIT License 
%  
%  Revised: Hua Huang <huangh223@gatech.edu>, 2019
%  
%******************************************************************

clear all
Farraytmp=load('Farray.mat');
Farray=Farraytmp.Fvunum;


deletelimit=0;
[species,xyz] = findgeomgjf('gaussian_testjob\CH3MO.gjf'); % using reorient in gaussian before using this code, this can accelerate the calculation a lot. 
xyz=xyz*1.889725989;
[spreads,d,shapematrix,centers,Nelec,Nnuc,nucchg,K,L]=initialization_HF(species,xyz,'STO2G'); % The code only support STO2G and 321G basis set

internucEnergy=0;
for ix =1:Nnuc
    for iy=(ix+1):(Nnuc)
        distance = xyz(ix,:)- xyz(iy,:);
        internucEnergy=internucEnergy+nucchg(ix)*nucchg(iy)/(distance*distance')^0.5;
    end
end   
fprintf('internucEnergy = %d\n', internucEnergy);

%%
% Compute all single electron operations
S = zeros(K,K);
T = zeros(K,K);
V = zeros(K,K);
ERI_diag = zeros(K,K);
flag2 = zeros(K,K);

%Calculate the parts of the Fock matrix hamiltonian:
for mu = 1 : K 
for nu = 1 : K 
    if (flag2(mu, nu) ~= 0), continue; end
    
    for p=1:L(mu)
    for q=1:L(nu)
        RA = [centers(mu,1) centers(mu,2) centers(mu,3)];
        RB = [centers(nu,1) centers(nu,2) centers(nu,3)];
        alpha = spreads(mu,p);
        beta = spreads(nu,q);
        
        S(mu,nu) = S(mu,nu) + d(mu,p)*d(nu,q)*...
        myoverlap3(alpha,beta,shapematrix(mu,:),shapematrix(nu,:),RA,RB );
        
        T(mu,nu) = T(mu,nu) + d( mu,p)*d( nu,q)*...
        mykinetic(alpha,beta,shapematrix(mu,:),shapematrix(nu,:),RA,RB );
        
        for i = 1 : Nnuc
            RC = xyz(i,:);
            V(mu,nu) = V(mu,nu) + d(mu,p)*d(nu,q)*...
            nucchg(i)*mynuc_elec3(alpha,beta,shapematrix(mu,:),shapematrix(nu,:),RA,RB,RC,Farray);
        end
    end
    end
    
    % Screening
    if abs(T(mu,nu))<deletelimit
        T(mu,nu)=0;
    end
    if abs(V(mu,nu))<deletelimit
        V(mu,nu)=0;
    end
    if abs(S(mu,nu))<deletelimit
        S(mu,nu)=0;
    end
    
    % Force symmetry
    S(nu,mu)=S(mu,nu);
    T(nu,mu)=T(mu,nu);
    V(nu,mu)=V(mu,nu);
    flag2(mu,nu)=1;
    flag2(nu,mu)=1;
end
end
[orbs, Dorbs] = eig(S);
[orbs, Dorbs] = sorteig(orbs, Dorbs);    
X = zeros(size(orbs));
for i = 1 : K
    X(:, i) = orbs(:, i) / (Dorbs(i, i)^0.5);
end

%%
% compute angular momentum part of the 2-electron integral
% size of this temporary Efunc result matrix depend on the maximum
% angular momentum and contraction
efuncresultx=zeros(6,6,5,K,K);
efuncresulty=zeros(6,6,5,K,K);
efuncresultz=zeros(6,6,5,K,K);
for phi1=1:K
for phi2=1:K
    shape1=squeeze(shapematrix(phi1,:));
    shape2=squeeze(shapematrix(phi2,:));
    L1=shape1(1);
    L2=shape2(1);
    m1=shape1(2);
    m2=shape2(2);
    n1=shape1(3);
    n2=shape2(3);
    A=[centers(phi1,1) centers(phi1,2) centers(phi1,3)];
    B=[centers(phi2,1) centers(phi2,2) centers(phi2,3)];
    Ax=A(1);
    Bx=B(1);
    Ay=A(2);
    By=B(2);
    Az=A(3);
    Bz=B(3);

    alphax1 = repmat(spreads( phi1,1:L(phi1))',[1,L(phi2)]) ;
    alphax2 = repmat((spreads( phi2,1:L(phi2))),[L(phi1),1]) ;
    for t1=0:(L1+L2)
        efuncresultx(1:L(phi1),1:L(phi2),t1+1,phi1,phi2)=Efunc(t1,L1,L2,alphax1,alphax2,Ax,Bx);
    end
    for u1=0:(m1+m2)
        efuncresulty(1:L(phi1),1:L(phi2),u1+1,phi1,phi2)=Efunc(u1,m1,m2,alphax1,alphax2,Ay,By);
    end    
    for v1= 0:(n1+n2)
        efuncresultz(1:L(phi1),1:L(phi2),v1+1,phi1,phi2)=Efunc(v1,n1,n2,alphax1,alphax2,Az,Bz);
    end
end
end

for ix =0:4
    for iy=0:4
        for iz=0:4
            effuncall(:,:,ix+5*(iy)+25*(iz)+1,:,:) = efuncresultx(:,:,ix+1,:,:).*efuncresulty(:,:,iy+1,:,:).*efuncresultz(:,:,iz+1,:,:);
        end
    end
end

%% 
% Rename some variables 
natom     = Nnuc;         % Total number of atoms
nocc      = Nelec / 2;    % Total number of occupied orbitals
atom_xyz  = xyz;          % Atom coordinates
nbf       = K;            % Total number of basis functions
bf_coef   = d;            % coef  terms  of basis functions
bf_alpha  = spreads;      % alpha terms  of basis functions
bf_exp    = shapematrix;  % Polynomial exponents terms of basis functions
bf_center = centers;      % Center of basis functions
bf_nprim  = L;            % Number of primitive functions in each basis function
Hcore     = T + V;        % Core Hamiltonian

%%
% Two-electron repulsive integral
% Note: in real-world calculation, we cannot store the 4D ERI tensor,
% shell quartets need to be computed in each SCF iteration repeatly
ERI  = zeros(nbf, nbf, nbf, nbf);
flag = zeros(nbf, nbf, nbf, nbf);
tic;
for i = 1 : nbf
for j = 1 : nbf
for k = 1 : nbf
for l = 1 : nbf
    % 8-way symmetry property
    if (flag(i, j, k, l) == 1), continue; end

    % primitive screening
    if (abs(ERI(i, j, k, l)) < 1e-14)
        ERI(i, j, k, l) = 0;
    end

    ERI(i, j, k, l) = get2elec_tablevE(i,j,k,l,centers,spreads,shapematrix,d,L,Farray,effuncall);
    
    ERI(j, i, k, l) = ERI(i, j, k, l);
    ERI(j, i, l, k) = ERI(i, j, k, l);
    ERI(i, j, l, k) = ERI(i, j, k, l);
    ERI(l, k, i, j) = ERI(i, j, k, l);
    ERI(k, l, i, j) = ERI(i, j, k, l);
    ERI(l, k, j, i) = ERI(i, j, k, l);
    ERI(k, l, j, i) = ERI(i, j, k, l);
 
    flag(j, i, k, l) = 1;
    flag(j, i, l, k) = 1;
    flag(i, j, l, k) = 1;
    flag(l, k, i, j) = 1;
    flag(k, l, i, j) = 1;
    flag(l, k, j, i) = 1;
    flag(k, l, j, i) = 1;
end
end
end   
end
ut = toc;
fprintf('ERI tensor calculation = %.3f (s)\n', ut);

%% 
% Precompute the values of basis functions at XC integral points
% Note: in real-world calculation, we cannot store the 3D rho tensor
% for large molecules, rho need to be computed in each SCF iteration repeatly
tic;
[int_points, int_weights] = generate_integrate_weights_points(atom_xyz);
rho = calc_bf_value_at_int_points(int_points, atom_xyz, nbf, bf_coef, bf_alpha, bf_exp, bf_center, bf_nprim);
ut = toc;
fprintf('Precompute bf values at integral points = %.3f (s)\n', ut);

%% 
% SCF iteration

% Initial guess for density matrix: use core Hamiltonian
Hprime = X' * Hcore * X;
[Cprime, diag] = eig(Hprime);
[Cprime, diag] = sorteig(Cprime, diag);
C = X * Cprime;
C = C(:, 1 : nocc);
D = C * C';

% SCF iteration
iter   = 0;
deltaE = 1; 
while (iter < 100)
    tic;
    iter = iter + 1;
    
    % Constrct the Coulomb matrix
    J = zeros(nbf, nbf);
    for i = 1 : nbf
    for j = 1 : nbf
        for k = 1 : nbf
        for l = 1 : nbf
            J(i, j) = J(i, j) + D(k, l) * ERI(i, j, l, k);
        end
        end
    end
    end
    
    % Constrct the exchang-correlation matrix
    XC = eval_Xalpha_XC_with_rho(natom, nbf, rho, int_weights, D);
    
    % Constrct the complete Fock matrix
    H = 2 * J + XC;
    F = Hcore + H;
    Fprime = X' * F * X;
    
    % Construct density matrix using eigendecomposition
    [Cprime, diag] = eig(Fprime);
    [Cprime, diag] = sorteig(Cprime, diag);
    C = X * Cprime;
    C = C(:, 1 : nocc);
    D = C * C';
    
    % Density matrix mixing
    Dnew = D;
    if (iter > 1)
        mix_ratio = 0.7;
        D = Dold * mix_ratio + Dnew * (1 - mix_ratio);
    end
    if (iter > 10)
        mix_ratio = 0.5;
        D = Dold * mix_ratio + Dnew * (1 - mix_ratio);
    end
    Dold = D;
    
    % Calculate new energy
    energy = sum(sum(D .* (F + XC + Hcore)));
    energy = energy + internucEnergy;
    if (iter > 2)
        deltaE = energyold - energy;
    end
    if (abs(deltaE) < 1e-11)
        break
    end
    energyold = energy;
    ut = toc;
    fprintf('Iter %2d, energy = %d, deltaE = %e, %.3f (s)\n', iter, energy, deltaE, ut);
end

 