function uMain = PEsolvePML(PEpars)

% equation iU_x + b*U_{zz} + Vu + zeta = 0

% U = U(x,z) -- unknown function
% |- b(x,z)     --\
% |- V(x,z)     ---} equation coefficients
% |- zeta(x,z)  --/

x = PEpars.x;
z = PEpars.z;
dx = x(2)-x(1);
dz = z(2)-z(1);
nx = length(x);
nz = length(z);

%% PML
% PML is basically imaginary part of b!
% we replace b with b1 = b + Q, where Q is 
% Q = 0 for z_l < z < z_r
% Q = Q(z) = -i\sigma (z - z_r)^3, z > z_r
% z_r is the right boundary of the domain
% and similarly for z < z_l


sigma = 400;
if isfield(PEpars,'sigma');
    sigma = PEpars.sigma;
end;

lThickness = 100;
if isfield(PEpars,'lThickness');
    lThickness = PEpars.lThickness;
end;

nPML = max( 1,fix(lThickness/dz) );

if nPML > 1
    lShape(1:nPML,1) = - sigma*1i*(0:nPML-1).^3/((nPML - 1)^3);
    PML(1:nz,1) = 0;
    PML(nz-nPML+1:nz,1) = lShape(1:nPML,1);
    PML(1:nPML,1) = flipdim(lShape(1:nPML,1), 1);
end;

% figure;
% plot(z,imag(PML),'linewidth',2);
% title('Q(z) -- imaginary part of b_1(z)')


%%

V = PEpars.V;
b = PEpars.b;
zeta = PEpars.zeta;




uMain(1:nz,1:nx) = 0;
uMain(1:nz,1) = PEpars.ic;

% equation solution

vRHS(1:nz,1) = 0;       % right hand side vector


for ii = 2:nx
    
    % assembling system natrix and right hand side vector
    
    %% PML construction
    
    
    
    bPML(1:nz,1) = b(ii) + PML(1:nz,1);
%     if ii == 2
%         figure;
%         plot(z,b(ii)*( z*0 + 1  ),'linewidth',2,'color','red');
%         hold on;
%         plot(z,real(bPML(1:nz,1)),'linewidth',2,'color','blue','linestyle',':');
%         
%         plot(z,imag(bPML(1:nz,1)),'linewidth',2,'color','green');
%         
%         
%         legend('b(z) -- original b','\Re(b_{PML}(z)) -- real part of b_{PML} (same as original b(z))','\Im(b_{PML}(z)) -- imaginary part to adsorb waves');
%         
%         
%     end;
    
    %% END PML construction

    dLower(1:nz,1) = -1i*bPML(1:nz,1)/(dz^2);
    dUpper(1:nz,1) = -1i*bPML(1:nz,1)/(dz^2);
    dMain(1:nz,1) = 2*1i*bPML(1:nz,1)/(dz^2) + 2/dx - 1i*V(1:nz,ii);
    
    Mc = spdiags([dLower dMain dUpper], -1:1, nz, nz);

    %vRHS(1:nz,1) = [0; uMain(1:nz-1,ii-1)]*1i*b(ii)/(dz^2) + (-2*1i*b(ii)/(dz^2) + 2/dx + 1i*V(1:nz,ii) ).*uMain(1:nz,ii-1) + [uMain(2:nz,ii-1); 0]*1i*b(ii)/(dz^2);
    
    vRHS(1:nz,1) = (-2*1i*bPML(1:nz,1)/(dz^2) + 2/dx + 1i*V(1:nz,ii) ).*uMain(1:nz,ii-1);
    vRHS(2:nz,1) = vRHS(2:nz,1) + 1i*uMain(1:nz-1,ii-1).*bPML(2:nz,1)/(dz^2);
    vRHS(1:nz-1,1) = vRHS(1:nz-1,1) + 1i*uMain(2:nz,ii-1).*bPML(1:nz-1,1)/(dz^2);
    
    if ~(isempty(zeta))
        vRHS(1:nz,1) = ( vRHS(1:nz,1) + 1i*zeta(1:nz,ii) + 1i*zeta(1:nz,ii-1) );
    end;
    
    
    
    % solving linear system
    %-----------------------
    
    uMain(1:nz,ii) = Mc\vRHS;
    
    %-----------------------
        
   
end;



