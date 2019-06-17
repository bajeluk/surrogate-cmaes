% Based on the article by Wu (1995)
% 'The interaction of water waves with a group of submerged spheres'
% All variables have the same notations as in the paper
% Infinite depth
% Calculation of hydrodynamic coefficients and excitation force
% 04/06/2015 Algorithms is changed to solve everything in one big matrix (now coupled coefficients are calculated in a right way)
% Based on \Debugging\ArraySubmerged_v9.m
% 12/06/2015 Changed boundaries for the inetgrals evaluation in quadgk
% 19/06/2015 Modified for parallel computation, now frequency and
% wavenumber are not in the wave structure, but separate variables
% 31/08/2015 Added varargin
%           - varargin{1} - flag for the form of the output variables 0 - (l, i, j, k), 1 - as (numSphere*3 x numSphere*3) matrix
% 11/09/2015 Corrected the calculation of the excitation force

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT, Depends on the varargin{1}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Flag = 0
% Coefficients of the sphere l in mode i due to the motion of sphere k in mode j
% addedMass = myu(l,i,k,j)
% dampingCoef = lambda(l,i,k,j)
% excForce = F_exc(l,j)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Flag = 1, N - numSphere*3
% addedMass = A(N, N)
% dampingCoef = B(N, N)
% excForce = F(N, 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [addedMass, dampingCoef, excForce] = arraySubmergedSphereParfor(array, wave, frequency, waveNumber, numApprox, varargin)

if nargin == 5
    flag2D = 0;
else
    flag2D = varargin{1};
end

numSphere = array.number;
ro = wave.waterDensity;
% K = wave.waveNumber;
% w = wave.frequency;
K = waveNumber;
w = frequency;
alpha = wave.angle;

radiusSphere = array.radius;
xSphere = array.sphereCoordinate(1,:);
ySphere = array.sphereCoordinate(2,:);
zSphere = array.sphereCoordinate(3,:);

num_nm = (1+numApprox)/2*numApprox;
num_lnm = numSphere*num_nm;

% Matrices identification
coef_AB = zeros(2*num_lnm);
fh = zeros(2*num_lnm, numSphere, 3);
fhX = zeros(2*num_lnm, 1);
myu = zeros(numSphere, 3, numSphere, 3);
lambda = zeros(numSphere, 3, numSphere, 3);
F_exc = zeros(numSphere, 3);

M = containers.Map('KeyType','double','ValueType','any');
fctrM = containers.Map('KeyType','double','ValueType','double');
lgM = containers.Map('KeyType','double','ValueType','any');

for l = 1:numSphere
    
    sphere_l = [xSphere(l) ySphere(l) zSphere(l)]';
    a_l = radiusSphere(l);
    z_l = zSphere(l);
    y_l = ySphere(l);
    x_l = xSphere(l);

    exp_t = exp(K*z_l+1i*K*(x_l*cos(alpha)+y_l*sin(alpha)));
    
    for n = 0:numApprox-1

        for m = 0:n
            % After Eq.(12)
            if m == 0
                eps_m = 1;
            else
                eps_m = 2;
            end

            % Index (row) of the element for coefficients in front of A and B
            ind_lnm1 = num_nm*(l-1)+(1+n)/2*n+m+1;  % Eq.(19a)
            ind_lnm2 = num_lnm+num_nm*(l-1)+(1+n)/2*n+m+1;  % Eq.(19b)

            coef_AB(ind_lnm1,ind_lnm1) = coef_AB(ind_lnm1,ind_lnm1) -(n+1)/a_l;
            coef_AB(ind_lnm2,ind_lnm2) = coef_AB(ind_lnm2,ind_lnm2) -(n+1)/a_l;

            for k = 1:numSphere
                % Eq.(20)
                fh(ind_lnm1,k,1) = -(n==1)*(m==1)*(k==l);
                fh(ind_lnm1,k,2) = 0;
                fh(ind_lnm1,k,3) =  (n==1)*(m==0)*(k==l);
                fh(ind_lnm2,k,1) = 0;
                fh(ind_lnm2,k,2) = -(n==1)*(m==1)*(k==l);
                fh(ind_lnm2,k,3) = 0;
            end

            % Eq.(29a)-(29b)
            if isKey(fctrM,n+m)
                fnpm = fctrM(n+m);
            else
                fnpm = factorial(n+m);
                fctrM(n+m) = fnpm;
            end
            fhX(ind_lnm1,1) = n*eps_m*w*(-1i)^(m+1)*exp_t*cos(m*alpha)*(K*a_l)^(n-1)/fnpm;
            fhX(ind_lnm2,1) = n*eps_m*w*(-1i)^(m+1)*exp_t*sin(m*alpha)*(K*a_l)^(n-1)/fnpm;

            for lam = 1:numSphere
                a_lam = radiusSphere(lam);
                sphere_lam = [xSphere(lam) ySphere(lam) zSphere(lam)]';
                vect = sphere_l - sphere_lam;

                % Eq.(17)
                z_lam = zSphere(lam);
                y_lam = ySphere(lam);
                x_lam = xSphere(lam);
                r_laml = norm(vect);
                Rlaml = sqrt((x_l-x_lam)^2+(y_l-y_lam)^2);
                
%                 if Rlaml > 500
%                     continue;
%                 end
                
                theta_laml = acos(vect(3)/r_laml);
                beta_laml = atan2(vect(2),vect(1));

                for nd = 0:numApprox-1 % n dash
                    % cached calculation
                    ind_lamnd = num_nm*(lam-1)+(1+nd)/2*nd;
                    coef_d = eps_m*n*a_lam^(nd+1)*a_l^(n-1);
                    
                    for md = 0:nd % m dash

                        % Index (column) of the element for coefficients in front of A and B
                        ind_lamndmdA = ind_lamnd+md+1;
                        ind_lamndmdB = num_lnm+ind_lamnd+md+1;

                        mdmm=md-m;
                        mdpm=md+m;
                        
                        % hashcode for the key
                        mparas = 173*(173*(173*(23*(n+nd)) + (z_l+z_lam)) + mdmm) + Rlaml;
                        pparas = 173*(173*(173*(23*(n+nd)) + (z_l+z_lam)) + mdpm) + Rlaml;
                        

%                             funM = @(x)x.^(n+nd).*(x+K)./(x-K).*exp(x*(z_l+z_lam)).*besselj(m-md, x*Rlaml);
%                             funP = @(x)x.^(n+nd).*(x+K)./(x-K).*exp(x*(z_l+z_lam)).*besselj(mdpm, x*Rlaml);
                        
                        
                        if isKey(M, mparas)
                            int_mMmd = M(mparas);
                        else
                            funM = @(x)x.^(n+nd).*(x+K)./(x-K).*exp(x*(z_l+z_lam)).*besselj(mdmm, x*Rlaml);
                            int_mMmd1 = quadgk(funM, 0, K*1.5, 'Waypoints', [K*0.7, K+0.1*K*1i, K*1.25]);
                            int_mMmd2 = integral(funM, K*1.5, inf);                                         % changed 12/06/2015
                            int_mMmd = int_mMmd1 + int_mMmd2; 
                            
%                             int_mMmd = apxIntegral(funM, 0, inf, 'Waypoints', [K*0.7, K+0.1*K*1i, K*1.25]);	% changed 12/06/2015
                            M(mparas) = int_mMmd;
                        end
                        
                        if (mdmm) == (mdpm)
                            int_mPmd = int_mMmd;
                        else
                            if isKey(M, pparas)
                                int_mPmd = M(pparas);
                            else
                                funP = @(x)x.^(n+nd).*(x+K)./(x-K).*exp(x*(z_l+z_lam)).*besselj(mdpm, x*Rlaml);
                                int_mPmd1 = quadgk(funP, 0, K*1.5, 'Waypoints', [K*0.7, K+0.1*K*1i, K*1.25]);
                                int_mPmd2 = integral(funP, K*1.5, inf);   
                                int_mPmd = int_mPmd1 + int_mPmd2;
%                                 int_mPmd = apxIntegral(funP, 0, inf, 'Waypoints', [K*0.7, K+0.1*K*1i, K*1.25]);	% changed 12/06/2015
                                M(pparas)=int_mPmd;
                            end

                        end
                        


%                             I11 = (-1)^m*pi*( cos((m-md)*beta_laml)*int_mMmd+(-1)^md*cos((mdpm)*beta_laml)*int_mPmd);
%                             I12 = (-1)^m*pi*(-sin((m-md)*beta_laml)*int_mMmd+(-1)^md*sin((mdpm)*beta_laml)*int_mPmd);
%                             I21 = (-1)^m*pi*( sin((m-md)*beta_laml)*int_mMmd+(-1)^md*sin((mdpm)*beta_laml)*int_mPmd);
%                             I22 = (-1)^m*pi*( cos((m-md)*beta_laml)*int_mMmd-(-1)^md*cos((mdpm)*beta_laml)*int_mPmd);
                        I11 = (-1)^md*pi*( cos((mdmm)*beta_laml)*int_mMmd+(-1)^m*cos((mdpm)*beta_laml)*int_mPmd);
                        I12 = (-1)^md*pi*(-sin((mdmm)*beta_laml)*int_mMmd+(-1)^m*sin((mdpm)*beta_laml)*int_mPmd);
                        I21 = (-1)^md*pi*( sin((mdmm)*beta_laml)*int_mMmd+(-1)^m*sin((mdpm)*beta_laml)*int_mPmd);
                        I22 = (-1)^md*pi*( cos((mdmm)*beta_laml)*int_mMmd-(-1)^m*cos((mdpm)*beta_laml)*int_mPmd);

                        if lam ~= l
                            %hashing
                            lgp = 173*(23*(n+nd)) + theta_laml;
                            
                            if isKey(lgM, lgp)
                                Pnnd = lgM(lgp);
                            else
                                Pnnd = legendre(n+nd,cos(theta_laml));
                                lgM(lgp)=Pnnd;
                            end

                            if mdpm > n+nd
                                PmPmd = 0;
                            else
                                PmPmd = Pnnd(mdpm+1);   % P_{n+nd}^{mdpm}
                            end

                            if (mdmm+1) > 0
                                PmdMm = Pnnd(mdmm+1); % P_{n+nd}^{mdmm}
                            else
                                if abs(mdmm)>nd+n
                                    PmdMm = 0;
                                else
                                    if isKey(fctrM,n+nd-(m-md))
                                        nndmmd = fctrM(n+nd-(m-md));
                                    else
                                        nndmmd = factorial(n+nd-(m-md));
                                        fctrM(n+nd-(m-md)) = nndmmd;
                                    end
                                    if isKey(fctrM,n+nd+(m-md))
                                        fnpm = fctrM(n+nd+(m-md));
                                    else
                                        fnpm = factorial(n+nd+(m-md));
                                        fctrM(n+nd+(m-md)) = fnpm;
                                    end
                                    PmdMm = (-1)^(m-md)*nndmmd/fnpm*Pnnd(abs(m-md)+1);
                                end
                            end

                            if isKey(fctrM,nd-md)
                                ndmd = fctrM(nd-md);
                            else
                                ndmd = factorial(nd-md);
                                fctrM(nd-md) = ndmd;
                            end
                            if isKey(fctrM,n+m)
                                nm = fctrM(n+m);
                            else
                                nm = factorial(n+m);
                                fctrM(n+m) = nm;
                            end
                            
                            coef0 = (-1)^(n+m)*coef_d/(2*r_laml^(n+nd+1)*ndmd*nm);
                            

                            if n+nd-m-md < 0
                                fact1 = 0;
                            else
                                if isKey(fctrM,n+nd-m-md)
                                    fact1 = fctrM(n+nd-m-md);
                                else
                                    fact1 = factorial(n+nd-m-md);
                                    fctrM(n+nd-m-md) = fact1;
                                end
%                                 fact1 = factorial(n+nd-m-md);
                            end
                            
                            if isKey(fctrM,n+nd-md+m)
                                fact2 = fctrM(n+nd-md+m);
                            else
                                fact2 = factorial(n+nd-md+m);
                                fctrM(n+nd-md+m) = fact2;
                            end
                            
%                             fact2 = cachedfactorial(fctrM,(n+nd-md+m));

                            coef_AB(ind_lnm1,ind_lamndmdA) = coef_AB(ind_lnm1,ind_lamndmdA)...
                                + coef0*((-1)^m*fact1*PmPmd*cos((mdpm)*beta_laml)+fact2*PmdMm*cos((m-md)*beta_laml));

                            coef_AB(ind_lnm1,ind_lamndmdB) = coef_AB(ind_lnm1,ind_lamndmdB)...
                                + coef0*((-1)^m*fact1*PmPmd*sin((mdpm)*beta_laml)-fact2*PmdMm*sin((m-md)*beta_laml));

                            coef_AB(ind_lnm2,ind_lamndmdA) = coef_AB(ind_lnm2,ind_lamndmdA)...
                                + coef0*((-1)^m*fact1*PmPmd*sin((mdpm)*beta_laml)+fact2*PmdMm*sin((m-md)*beta_laml));

                            coef_AB(ind_lnm2,ind_lamndmdB) = coef_AB(ind_lnm2,ind_lamndmdB)...
                                - coef0*((-1)^m*fact1*PmPmd*cos((mdpm)*beta_laml)-fact2*PmdMm*cos((m-md)*beta_laml));
                        end
                        
                        if isKey(fctrM,nd-md)
                            ndmd = fctrM(nd-md);
                        else
                            ndmd = factorial(nd-md);
                            fctrM(nd-md) = ndmd;
                        end
                        if isKey(fctrM,n+m)
                            nm = fctrM(n+m);
                        else
                            nm = factorial(n+m);
                            fctrM(n+m) = nm;
                        end
                        
                        coef1 = coef_d*(-1)^m/(2*pi*ndmd*nm);

                        coef_AB(ind_lnm1,ind_lamndmdA) = coef_AB(ind_lnm1,ind_lamndmdA) + coef1*I11;
                        coef_AB(ind_lnm1,ind_lamndmdB) = coef_AB(ind_lnm1,ind_lamndmdB) + coef1*I21;
                        coef_AB(ind_lnm2,ind_lamndmdA) = coef_AB(ind_lnm2,ind_lamndmdA) + coef1*I12;
                        coef_AB(ind_lnm2,ind_lamndmdB) = coef_AB(ind_lnm2,ind_lamndmdB) + coef1*I22;
                    end % md
                end % nd
            end % lam
        end % m
    end % n
end % l

% disp(M.keys);
% disp(M.values);


for k = 1:numSphere

    fh_k(:,:) = fh(:,k,:);
    AB = coef_AB\fh_k;        

    for l = 1:numSphere
        a_l = radiusSphere(l);
        A11 = AB(num_nm*(l-1)+3,:);
        A01 = AB(num_nm*(l-1)+2,:);
        B11 = AB(num_lnm+num_nm*(l-1)+3,:);

        % Eq.(23a)-(23c)
        myu(l,1,k,:) =  real( 4/3*ro*pi*a_l^2*(3*A11-a_l*(k==l)*[1 0 0]));
        myu(l,2,k,:) =  real( 4/3*ro*pi*a_l^2*(3*B11-a_l*(k==l)*[0 1 0]));
        myu(l,3,k,:) =  real(-4/3*ro*pi*a_l^2*(3*A01+a_l*(k==l)*[0 0 1]));

        lambda(l,1,k,:) = -imag( 4/3*ro*pi*a_l^2*(3*A11-a_l*(k==l)*[1 0 0]))*w;
        lambda(l,2,k,:) = -imag( 4/3*ro*pi*a_l^2*(3*B11-a_l*(k==l)*[0 1 0]))*w;
        lambda(l,3,k,:) = -imag(-4/3*ro*pi*a_l^2*(3*A01+a_l*(k==l)*[0 0 1]))*w;
    end % l
end % k

ABX = coef_AB\fhX;

ind1 = 3:num_nm:num_lnm;	% A11(7)
ind2 = ind1 + num_lnm;      % B11(7)
ind3 = ind1-1;              % A01(7)

%% modified 11/09/2015
for l = 1:numSphere
    a_l = radiusSphere(l);
    % Eq.(32a)-(32c)
    F_exc(l,1) = -4*ro*1i*w*pi*a_l^2*ABX(ind1(l));
    F_exc(l,2) = -4*ro*1i*w*pi*a_l^2*ABX(ind2(l));
    F_exc(l,3) =  4*ro*1i*w*pi*a_l^2*ABX(ind3(l));
end

%%

if flag2D == 0
    addedMass = myu;
    dampingCoef = lambda;
    excForce = F_exc;
else
    A = zeros(numSphere*3);
    B = zeros(numSphere*3);
    X = zeros(numSphere*3, 1);
    % Reshaping hydrodynamic coefficients in 2D matrices
    for ii = 1:numSphere
        for jj = 1:numSphere
            A(3*(ii-1)+1:3*ii,3*(jj-1)+1:3*jj) = squeeze(myu(ii,:,jj,:));
            B(3*(ii-1)+1:3*ii,3*(jj-1)+1:3*jj) = squeeze(lambda(ii,:,jj,:));
        end
        X(3*(ii-1)+1:3*ii,1) = F_exc(ii,:).';
    end
    addedMass = A;
    dampingCoef = B;
    excForce = X;
end
end

