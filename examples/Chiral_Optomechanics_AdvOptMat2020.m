close all;
clearvars; clc

% Reproduces figures in "Optomechanics of chiral dielectric metasurfaces",
%               Adv. Optical Mater. 2020, 8, 1901507 (2020). DOI: 10.1002/adom.201901507 
%
%  Simone Zanotto, Pisa, oct. 2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is free software distributed under the BSD licence (see the 
%  containing folder).
% However, shall the results obtained through this code be included 
%  in an academic publication, we kindly ask you to cite the source 
%  website and, if applicable, the following paper:
%
% Simone Zanotto, Alessandro Tredicucci, Daniel Navarro-Urrios, 
% Marco Cecchini, Giorgio Biasiol, Davide Mencarelli, Luca Pierantoni,
% Alessandro Pitanti, "Optomechanics of chiral dielectric metasurfaces",
%               Adv. Optical Mater. 2020, 8, 1901507 (2020). DOI: 10.1002/adom.201901507 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('PPML_root')    % to be replaced by the proper path 

epsOx   = (1.5)^2;%
dOx     = 5;
dGaAs   = 210.5; %

lambdav = [1520:2:1572];

a   = 1128*1544.5/1548;          
ff1 = 0.661;
ff2 = 0.817;
ff3 = 0.372;
ff4 = 0.296;
%%%%%%%%%%%

theta = 0.01;     % angle in degrees (never set 0)
phi   = 90;       % angle in degrees 

parfor i = 1:length(lambdav)
% Computation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i
lambda = lambdav(i);

halfnpw = 5;       
    
epsGaAs = 3.38^2;
    
L = 3;     % number of internal layers 

%            air  |  oxide  |  patt. GaAs |  oxide |   air          
f1      = [          ff1       ff1           ff1            ]; % fraction of B in A
f2      = [          ff2       ff2           ff2            ]; % fraction of B in A
f3      = [          ff3       ff3           ff3            ]; % fraction of B in A
f4      = [          ff4       ff4           ff4            ]; % fraction of B in A
d       = [    0     dOx       dGaAs         dOx     0      ]; % nm
                                       
epsA    = [           epsOx    epsGaAs       epsOx          ]; % material A 
epsB    = [           1        1             1              ]; % material B
sigma2d = [       0         0            0           0      ];

epssup = 1;                                
epssub = 1;


k0   = 2*pi/lambda;         % wavevector in nm ^-1
kparx = k0*sin(theta*pi/180)*cos(phi*pi/180);
kpary = k0*sin(theta*pi/180)*sin(phi*pi/180);

% calling the scattering matrix subroutine for reflection
S = ZSM_2d_Lshape(a,L,...
   epssup,epssub,epsA,epsB,sigma2d,f1,f2,f3,f4,d,...
   halfnpw,k0,kparx,kpary);
Sv{i} = S;

end


clearvars a d dGaAs dOx epsA epsB epsGaAs epsOx epssup epssub f1 f2 f3 ff1 ff2 ff3
clearvars halfnpw i j k0 k_GaAs kparx kpary L lambda_GaAs n_GaAs lambda S
%%
dspac = 1455; 
dsub = 336180; %dsub = 360060; 
res = 0.01;
[R11,  T16, lambdavref] = trilay(Sv,lambdav,dspac,dsub,res);
dspac = 1456; % 1 nm diffferential step
[R11d, T16d, waste]  = trilay(Sv,lambdav,dspac,dsub,res);

U = [1 1; 1 -1]*(1/sqrt(2));

pols = {'R', 'L', 'H', 'V', 'HV'};
alphas = {[1; 1i]/sqrt(2), [1; -1i]/sqrt(2), [0; 1],      [1; 0],      [0; 1]    };
Fs     = {diag([1,1]),     diag([1,1]),      diag([1,1]), diag([1,1]), diag([1,0]) };

for i = 1:length(lambdavref)
R  = (U')*R11{i}*U; 
Rp = (U')*(R11d{i} - R11{i})*U; 

Lambda = (1/sqrt(2))*[1,  1; 1i, -1i];
Rcirc = (Lambda')*R*Lambda;

for j = 1:5
alpha = alphas{j}; F = Fs{j};
rho(i,j)   = (alpha')*(R')*F*R*alpha;
sigma(i,j) = (alpha')*( (R')*F*Rp + (Rp')*F*R )*alpha; 
end
end

sig = 40; % numero di step di lambdavref che stanno in una sigma
gausvec = [1:1:(6*sig+1)];
gaus = (exp(-0.5*((gausvec-3*sig-1)/sig).^2))/sig/sqrt(2*pi);

figure(1)
RHV = real(rho(:,5));
RHVsmo = conv(RHV,gaus,'same');
plot(lambdavref,RHV,'g'); hold on;
plot(lambdavref,RHVsmo,'b'); hold on;
title(['Cross polarized reflection'])
xlabel('Wavelength(nm)')
axis([1520 1560 0 0.22]); grid on


figure(2)
CD = real(rho(:,2) - rho(:,1));
CDsmo = conv(CD,gaus,'same');
plot(lambdavref,lambdavref*0,'k'); hold on;
plot(lambdavref,CD,'g'); hold on;
plot(lambdavref,CDsmo,'b'); hold on;
title(['Circular dichroism (CD)'])
xlabel('Wavelength(nm)')
axis([1520 1560 0 0.22]); grid on


figure(3)
plot([1520 1560],[0 0],'k'); hold on;
plot(lambdavref,real(sigma(:,3)),'b'); hold on;
plot(lambdavref,real(sigma(:,5)),'c'); hold on;
axis([1520 1560 -1e-2 2e-2]); grid on
title(['Intensity modulation at detector'])
xlabel('Wavelength(nm)')

figure(4)
CD = real(sigma(:,2) - sigma(:,1));
CDsmo = conv(CD,gaus,'same');
plot(lambdavref,lambdavref*0,'k'); hold on;
plot(lambdavref,CD,'g'); hold on;
plot(lambdavref,CDsmo,'b'); hold on;
title(['Modulation of CD'])
xlabel('Wavelength(nm)')
axis([1520 1560 -3e-3 3e-3]); grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [RR11, TT16, lambdavref] = trilay(Sv,lambdav,dspac,dsub,res)

Transm12v = zeros(length(lambdav),2,2);
Refl22v   = zeros(length(lambdav),2,2);
Transm21v = zeros(length(lambdav),2,2);
Refl11v   = zeros(length(lambdav),2,2);
for i = 1:length(lambdav)

Svv = diag([-1 1 -1  -1 ]) * Sv{i} * diag([-1 -1 -1 1]); % 

Transm12v(i,1,1) =  Svv(3,1);
Transm12v(i,1,2) =  Svv(3,2);
Transm12v(i,2,1) =  Svv(4,1);
Transm12v(i,2,2) =  Svv(4,2);
                           
Refl22v(i,1,1) =  Svv(3,3);
Refl22v(i,1,2) =  Svv(3,4);
Refl22v(i,2,1) =  Svv(4,3);
Refl22v(i,2,2) =  Svv(4,4);

Transm21v(i,1,1) =  Svv(1,3);
Transm21v(i,1,2) =  Svv(1,4);
Transm21v(i,2,1) =  Svv(2,3);
Transm21v(i,2,2) =  Svv(2,4);
                           
Refl11v(i,1,1) =  Svv(1,1);
Refl11v(i,1,2) =  Svv(1,2);
Refl11v(i,2,1) =  Svv(2,1);
Refl11v(i,2,2) =  Svv(2,2);


end

lambdavref = [min(lambdav):res:max(lambdav)];
for i = 1:length(lambdavref)

dn = -0.0014;    
epsGaAs = 3.38^2;
    
id = eye(2); 
n = sqrt(epsGaAs);
T43 = [id*(1+n), -id*(1-n);  -id*(1-n), id*(1+n) ];
T65 = [id*(1+n),  id*(1-n);   id*(1-n), id*(1+n) ];

phisub = 2*pi*n*dsub/lambdavref(i);
phispac = 2*pi*dspac/lambdavref(i);

T54   = [id*exp(1i*phisub), id*0; id*0, id*exp(-1i*phisub)];
T32   = [id*exp(1i*phispac), id*0; id*0, id*exp(-1i*phispac)];
T = (1/4/n)*T65*T54*T43*T32;
alpha = T(1:2,1:2);
beta  = T(1:2,3:4);
gamma = T(3:4,1:2);
delta = T(3:4,3:4);

Transm12 = reshape(interp1(lambdav,Transm12v,lambdavref(i)),[2,2]);
Refl22   = reshape(interp1(lambdav,Refl22v  ,lambdavref(i)),[2,2]);
Transm21 = reshape(interp1(lambdav,Transm21v,lambdavref(i)),[2,2]);
Refl11   = reshape(interp1(lambdav,Refl11v  ,lambdavref(i)),[2,2]);
TT16{i} = (alpha-beta*(delta\gamma))*( (id + Refl22*(delta\gamma)) \Transm12 );
RR11{i}   = Refl11 - Transm21*(delta\gamma)*( (id + Refl22*(delta\gamma)) \Transm12 );


end
end