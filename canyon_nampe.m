close all;
clear variables;
clc;

% Pavel Petrov, 2018.11.03
% A solution of the canyon problem (Canyon Test case 1) using Collins-style
% mode parabolic equations we use formula (5.320) from 
% Jensen, Porter, Kuperman, Schmidt "Computational ocean acoustics" (2011)
% hereafter referred to as JPKS2011
% we make a (paraxial/narrow-angle) parabolic equation out of it 
% and solve it numerically


set(0, 'DefaultAxesFontSize', 16, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontSize', 16, 'DefaultTextFontName', 'Arial'); 

% source

freq = 40;
zs = 10;

omeg = 2*pi*freq;

%domain & bathymetry

x = 0:2:30000;
y = -2000:1:2000;

nx = length(x);
ny = length(y);
iy0 = ceil(ny/2);


h0 = 50;
dh = 5;
sigm = 0.005;

hy = h0 + dh*( (sech(sigm*y)).^2 );

figure;
plot(y,hy,'color','black','linewidth',2);
set(gca,'YDir','reverse');
grid on;

hmax = ceil(max(hy));
hmin = floor(min(hy));

hs = hmin:1:hmax;
nh = length(hs);


% media parameters

cw = 1500; % sound speed in the water
cb = 1700; % -- in the bottom

rhow = 1;   % density in the water
rhob = 1.7; % -- in the bottom

gammaw = 1/rhow;
gammab = 1/rhob;



betab = 0; % attenuation in dB/lambda

%% Solve spectral problem in order to obtain the values of k_j(y) for all y

% we find minimal and maximal values of h (hmin/hmax), then we solve the
% spectral problem for all values between hmin and hmax with a step of 1 m.
% then we do linear interpolation of onto all elements of vector y (h(y)).

% spectral problem solution opts

dz0 = 0.25;
opts.Tgr = 3;
opts.Ngr = 3;
opts.nmod = -1;
opts.Hb = 200;
opts.BotBC = 'D';

kh = [];

for ii = 1:nh
    
    hi = hs(ii);
    
    MP = [[0    cw  cw  rhow    rhow    0   0];
          [hi   cw  cb  rhow    rhob    0   betab]
        ];
    
    % solve spectral problem
    
    if ii == 1
        [krs, wmode, dwmode] = ac_modesr(dz0,MP,freq,opts);
        
        nmod = length(krs);
        kh = zeros(nh,nmod); 
        
        
        phihs = zeros(size(kh));
        phihzs = zeros(size(kh));
        
        z = dz0*(0:size(wmode,1)-1);
        
        izs = find(z>=zs,1,'first');
    else
        
        [krs, wmode, dwmode] = ac_modesr(dz0,MP,freq,opts);
    end;
    
    % improve accuracy by using some perturbation theory 
    
    % this is simplest way: use Eq.(5.176) from JPKS2011
    wnum_im_part = ModesAttCoeffs(dz0,freq,krs,wmode,MP);
    
    wnum_im_part1 = PekerisCorrection(krs,MP,freq);
    
    
    disp('h=');disp(hi);
    
    % how accurate is just a real k_j obtained from spactral problem
    % solution?
    
    disp('k_0='); disp(krs(1:nmod));
    err2 = ModesAccuracyCheckPekerisComplex(krs(1:nmod),MP,freq);
    disp(err2);
    
    % how accurate is correction according to (5.176)?
    % let us substitute complex wavenumber into Pekeris dispersion relation
    disp('k_int='); disp(krs(1:nmod) + 1i*wnum_im_part(1:nmod));
    err = ModesAccuracyCheckPekerisComplex(krs(1:nmod) + 1i*wnum_im_part(1:nmod),MP,freq);
    disp(err);
    
    
    % if necessary, we can use a correction by perturbation theory for
    % the Pekeris dispersion relation
    disp('k_ptb='); disp(sqrt(krs(1:nmod).^2 + wnum_im_part1(1:nmod)));
    err1 = ModesAccuracyCheckPekerisComplex(sqrt(krs(1:nmod).^2 + wnum_im_part1(1:nmod)),MP,freq);
    disp(err1);
    
    % one more step, should be very accurate, check "err" output
    
    wnum_im_part2 = PekerisCorrection(sqrt(krs(1:nmod).^2 + wnum_im_part1(1:nmod)),MP,freq);
    disp('k_ptb2='); disp(sqrt(krs(1:nmod).^2 + wnum_im_part1(1:nmod) + wnum_im_part2(1:nmod) ));
    err3 = ModesAccuracyCheckPekerisComplex(sqrt(krs(1:nmod).^2 + wnum_im_part1(1:nmod) + wnum_im_part2(1:nmod) ),MP,freq);
    disp(err3);
    
    % save the wavenumbers into kh (k_j=k_j(h))
    
    % simplest way, just a correction according to (5.176)
    kh(ii,1:nmod) = krs(1:nmod) + 1i*wnum_im_part(1:nmod);
    
    % a more accurate way: few steps of the perturbation theory for the
    % Pekeris dispersion relation
    kh(ii,1:nmod) = sqrt(krs(1:nmod).^2 + wnum_im_part1(1:nmod) + wnum_im_part2(1:nmod) );
    
    phihs(ii,1:nmod) = wmode(izs,1:nmod);
    
    
    
end;

%% Computing the field by solving mode parabolic equations (MPEs) one by one
P(1:ny,1:nx) = 0;
Py0(1,1:nx) = 0;

for jj = 1:nmod
    
    % interpolate to obtain k_j=k_j(y) from k_j(h) and h(y)
    
    if length(hs) > 1
        ky = interp1(hs,kh(:,jj),hy);
    else
        ky(1:ny) = kh(1,jj);
    end;

    
    % reference wavenumber
    
    k0 = real(ky(iy0));
    
    
    % preparing input for PESolve routine
    
    PEpars.x = x;
    PEpars.z = y;
    
    PEpars.b(1:nx) = 1/(2*k0);
    PEpars.V(1:ny,1:nx) = repmat(((ky.^2 - k0^2)/(2*k0)).',1,nx);
    PEpars.zeta(1:ny,1:nx) = 0;
    
    PEpars.rightBC = 'transparent';
    PEpars.leftBC = 'transparent';
    PEpars.lThickness = 1000;
    PEpars.sigma = 300;
    
    w = 1/k0;
    phisrc = phihs(end,jj);
    
    % initial condition (or starter) approximating point source, see
    % Petrov, Sturm, JASA 2016, V.139, p.1343
    
    PEpars.ic = (phisrc/(2*sqrt(pi)))*exp( - ((y-y(iy0)).^2)/(w^2) );  
    
    A = PEsolvePML(PEpars);
    
    % the modes should be also y-dependent, but here it doesn't matter
    
    P = P + A.*repmat(exp(1i*k0*x),ny,1)*phisrc;
    
    
    % summing up all modes to obtain the total field
    Py0 = Py0 + A(iy0,1:nx).*exp(1i*k0*x)*phisrc;
    
    
end;


figure;
imagesc(x/1000,y/1000,trlo(4*pi*P));
caxis([-90 -55]);
ylim([-1 1]);
grid on;
xlabel('x, km');
ylabel('y, km');


figure;
plot(x/1000,trlo(4*pi*Py0));
ylim([-80 -40]);
xlabel('x, km');




dlmwrite('TL_can_AMPE.txt',[x.'/1000 trlo(4*pi*Py0).'],'\t');
