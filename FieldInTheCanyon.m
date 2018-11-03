clear variables;
%close all;
clc;

%% Defining params

%  grid:
zm = 1200;       % domain size, m       
xm = 30000;
ym = 2000;

nz = 1201;       % number of grid points
ny = 1000;
nx = 5001;


% Bottom parameters

cw = 1500;      % sound speed in water, m/s
cb = 1700;      % sound speed in bottom, m/s
rhow = 1;       % density in water, g/cm^3
rhob = 1.7;       % density in bottom, g/cm^3
attb = 0.5;     % bottom attenuation, dB/wavelength
H0 = 55;
% h_0, m
H1 = 2000;

% Canyon parameters

DH = 5;
sigma = 0.005;

% Source

ys = 0;
zs = 10;        % source depth, m
freq = 40;      % frequency, Hz

% Receiver

zr = 10;        % receiver depth, m
yr = 0;


% Wavenumbers (Pekeris disp. relation is satisfied up to 10^(-10) )
% if cw, cb, rhow, rhob, h0 or freq are modified, kj must be recalculated!

K = [     0.162157190084758 % 0.161297680458947
   0.147830944703158
   0.147804814006921
   0.147761256373868
   0.147700262144574
   0.147621817898414
   0.147525906564526
   0.147412507576635
   0.147281597081880
   0.147133148216671
   0.146967131466566
   0.146783515131399
   0.146582265922675
   0.146363349726304
   0.146126732571386
   0.145872381853698
   0.145600267869848
   0.145310365723182
   0.145002657659782
   0.144677135874035 ];
    





%% Computing the field

% some auxillary variables...

gammaw = 1/rhow;
gammab = 1/rhob;

eta = 1/(40*pi*log10(exp(1)));
nub = 2*1i*eta*attb*(2*pi*freq/cb)^2;
omeg = 2*pi*freq;


x = linspace(0,xm,nx);
y = linspace(-ym,ym,ny);
z = linspace(0,zm,nz);

nmod = length(K);

Pxz(1:nx,1:nz)=0;
Pxy(1:nx,1:ny)=0;
Px(1:nx)=0;

phi(1:nz) = 0;

aPtrack(1:nx) = 0;      % acoustic field along the track -- for y=0, z=zr
aPxy(1:nx,1:ny) = 0;    % acoustic field P(x,y,z_r) at z = zr -- 2D array

Xx = repmat(x',1,ny);
Yy = repmat(y,nx,1);


izBot = find(z<=H0,1,'last');

%% 1) Modal functions and their derivatives

B_jl(1:nmod,1:nmod) = 0;


phih0(1:nmod) = 0;
phih0z(1:nmod) = 0;
phizs(1:nmod) = 0;
phizr(1:nmod) = 0;

phi(1:nz,1:nmod) = 0;

for jj = 1:nmod
    k_j = K(jj);
    
    kw = omeg/cw;
    kb = omeg/cb;
    
    if k_j > kb
        
        % TRAPPED MODE: sin in water, exp in bottom
        
        kvw = sqrt( kw^2 - k_j^2 );
        kvb = sqrt( k_j^2 - kb^2 );
        
        % computing the norm to normalize mode functions
        
        % I1 = \int_0^{H_0}(\sin(k_w^v z))^2/\rho_w dz
        
        I1 = H0/(2*rhow) - sin(2*kvw*H0)/(4*kvw*rhow);
        
        % I3 = \int_{H0}^{H1}(  a_3\exp(k^v_b (z-H_0) ) + b_3\exp( - k^v_b (z-H_0) ) )/\rho_b dz
        
        %a3 = ( sin(kvw*H0) + rhob*kvw*cos(kvw*H0)/(rhow*kvb) )/(2*rhob);
        
        a3 = 0;
        
        b3 = ( sin(kvw*H0) - rhob*kvw*cos(kvw*H0)/(rhow*kvb) )/(2*rhob);
        
        
        I3 = 2*a3*b3*(H1-H0) + 0.5*(b3^2 - a3^2)/kvb + ( (a3^2)*exp(2*kvb*(H1-H0) ) - (b3^2)*exp(2*kvb*(H0-H1) ) )/(2*kvb);
        
        Nrm = I1 + I3;
        
        
        phi(1:izBot,jj) = sin(kvw*z(1:izBot));
        phi(izBot+1:nz,jj) = a3*rhob*exp( kvb*( z(izBot+1:nz) - H0) ) + b3*rhob*exp( -kvb*( z(izBot+1:nz) - H0) );
        
    else
        
        % BOTTOM MODE: sin in water, sin in bottom
        
        kvw = sqrt( kw^2 - k_j^2 );
        kvb = sqrt( kb^2 - k_j^2 );
        
        % computing the norm to normalize mode functions
        
        % I1 = \int_0^{H_0}(\sin(k_w^v z))^2/\rho_w dz
        
        I1 = H0/(2*rhow) - sin(2*kvw*H0)/(4*kvw*rhow);
        
        % I2 = \int_{H0}^{H1}(  a_3\exp(k^v_b (z-H_0) ) + b_3\exp( - k^v_b (z-H_0) ) )/\rho_b dz
        
        I2 = ( (H1 - H0)/(2*rhob) - sin(2*kvb*(H1 - H0) )/(4*kvb*rhob) )*( ( sin(kvw*H0) )^2 )/( ( sin(kvb*(H0-H1)) )^2 )   ;
        
        Nrm = I1 + I2;
        
        phi(1:izBot,jj) = sin(kvw*z(1:izBot));
        phi(izBot+1:nz,jj) =  sin( kvb*(z(izBot+1:nz) - H1) )*( sin(kvw*H0) )/( sin(kvb*(H0-H1)) ) ;
        
    end;
    
    Nrm = sqrt(Nrm);
    
    phih0(jj) = sin( kvw*H0 )/Nrm;        % \psi_j(h_0)
    phih0z(jj) = kvw*cos( kvw*H0 )/Nrm;   % derivative \psi_{jz}(h_0)
    phizs(jj) = sin( kvw*zs )/Nrm;        % \psi_j(z_s)
    phizr(jj) = sin( kvw*zr )/Nrm;        % \psi_j(z_r)
    
    
    phi(1:nz,jj) = phi(1:nz,jj)/Nrm;
    
end;


figure;
hold all;
for jj=1:nmod
    plot(z, phi(1:nz,jj) );
end;


%% 2) Matrix Bjl -- mode interaction


for jj = 1:nmod
    for ll = 1:nmod
                    
        B_jl(jj,ll) =  ( (phih0(jj)*phih0(ll))*( gammaw*(omeg/cw)^2 - gammab*(omeg/cb)^2 ) +...
            (phih0(jj)*phih0(ll))*(gammab-gammaw)*K(ll)^2   -  (rhob - rhow)*phih0z(jj)*phih0z(ll)*(gammaw)^2 );    
    end;
end;

%disp(B_jl);

b11 = B_jl(1,1);

%% 3) Now we solve modal equation for A_1: compute \psi_{1l}(y_k) for all y

s1 = 0.5*(sqrt( 1 + 4*DH*b11/(sigma^2) ) - 1);
M1 = fix(s1);

Pm = 1;

psi(1:M1+1,1:ny) = 0;
psiys(1:M1+1) = 0;
psiyr(1:M1+1) = 0;
kH(1:M1+1) = 0;

for mm = 0:M1
    

    kappaHg = s1 - mm;
    aHg = kappaHg + s1 + 1;
    bHg = kappaHg + 1;
    
    
    coeffHg = Pm.*fliplr(cumprod([1 (0:mm-1)+aHg])./cumprod([1 (0:mm-1)+bHg]));
    
    coeffHg2 = conv(coeffHg,coeffHg);
    
    coeffHg2Xi = zeros(1,2*mm+1);
    zXi = 1;
    
    for tt = 1:2*mm + 1
        coeffHg2Xi(2*mm+1-tt+1:2*mm+1) = coeffHg2Xi(2*mm+1-tt+1:2*mm+1) + coeffHg2(2*mm+1-tt+1)*zXi;
        zXi = conv(zXi,[-0.5 0.5]);
    end;
    
    PsiNrm = 0;
    coeffHg2Xi = fliplr(coeffHg2Xi);
    
    for tt = 0:mm
        PsiNrm = PsiNrm + coeffHg2Xi(2*tt+1)*gamma(tt+0.5)*gamma(kappaHg)/gamma(kappaHg+tt+0.5);
    end;
    
    PsiNrm = sqrt(PsiNrm/sigma);
    
    psi(mm+1,1:ny) = (cosh(sigma*y).^(-kappaHg)).*polyval(coeffHg,0.5*(1-tanh(sigma*y)))/PsiNrm;
    
    psiys(mm+1) = (cosh(sigma*ys)^(-kappaHg))*polyval(coeffHg,0.5*(1-tanh(sigma*ys)))/PsiNrm;
    
    psiyr(mm+1) = (cosh(sigma*yr)^(-kappaHg))*polyval(coeffHg,0.5*(1-tanh(sigma*yr)))/PsiNrm;
    
    kH(mm+1) = sqrt(K(1)^2 + (sigma*kappaHg)^2 - b11*DH);
    
%     coeffHg2Xi
%     syms u v
%     aa = (subs(hypergeom([-mm aHg], bHg, u),u,0.5*(1-v)))^2;
%     if mm>1
%     expand(aa)
%     else
%         aa
%     end
    
    Pm = conv(Pm,[-1 1]);
    
    
end;



aFieldXY(1:nx,1:ny) = 0;

aFieldX(1:nx) = 0;


for kk = 1:M1+1
    
    q(1:nx,1) = 0.5*1i*psiys(kk)*phizs(1)*exp( 1i*kH(kk)*abs(x) )/kH(kk);
    
    aFieldXY(1:nx,1:ny) = aFieldXY(1:nx,1:ny) + phizr(1)*repmat(psi(kk,:),nx,1).*repmat(q,1,ny);
    
end;

set(0, 'DefaultAxesFontSize', 16, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontSize', 16, 'DefaultTextFontName', 'Arial');   

figure;
hold all;
for kk=1:M1+1
    plot(y,psi(kk,:));
    disp( sum(psi(kk,:).*psi(kk,:)*(y(2)-y(1))) );
end;


figure;
imagesc(x/1000,y/1000,20*log10(abs(4*pi*aFieldXY)).');
ylim([-1 1]);
set(gca,'XTick',0:5:20);
set(gca,'YTick',-1:0.5:1);
xlabel('x, km');
ylabel('y, km');
%title('P(x,y,z_s), y_s = 0');
grid on;
caxis([-90 -55]);
colorbar('YTick',[-90 -80 -70 -60 -50]);



dlmwrite('TLx.txt',[x.'/1000 20*log10(abs(4*pi*aFieldXY(1:nx,ceil(ny/2))))],'\t');

