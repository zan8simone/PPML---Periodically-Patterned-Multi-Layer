%%%
function epar = epar_1d(a,L,...
   epssup,epssub,epsA,epsxxB,epsxyB,epsyyB,epszB,sigma,...
   f,d,halfnpw,k0,kparx,kpary,pol)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is free software distributed under the BSD licence (see the 
%  containing folder).
% However, shall the results obtained through this code be included 
%  in an academic publication, we kindly ask you to cite the source 
%  website.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the electric field Fourier components
%
% 1-d structure with anisotropic material
%
% Simone Zanotto, Orsay Oct. - Dec 2012; Firenze Feb.-Dec 2016
% Gaia Da Prato, Pisa Feb-May 2021
%
% -------------------------------------------------------------------------
%                              INPUTS
%
% a              -> periodicity 
% L              -> number of layers excluded superstrate and substrate
% epssup, epssub -> super-and sub-strate permittivities 
% epsA           -> permittivity of material A in the internal layers
%                                         [vector with L components]
% epsxxB         -> dielectric tensor xx-component for material B in the internal layers    
%                                         [vector with L components]
% epsxyB         -> dielectric tensor xy-component for material B in the internal layers    
%                                         [vector with L components]
% epsyyB         -> dielectric tensor yy-component for material B in the internal layers    
%                                         [vector with L components]
% epszB          -> dielectric tensor z-component for material B in the internal layers  
%                                         [vector with L components]
% sigma          -> conductivity of interfaces
%                                         [vector with (L+1) components]
% f              -> duty cycle in the internal layers
%                                         [vector with L components]
% d              -> thicknesses of the layers, 
%                   including super- and sub-strate
%                                         [vector with (L+2) components]
% halfnpw         -> half number of harmonic waves used in calculation
% k0              -> wavevector in vacuum  (= 2*pi/lambda0)
% kparx           -> x-projection of the incident wavevector 
% kpary           -> y-projection of the incident wavevector 
% pol             -> polarization state of the incident light	[‘s’ or ‘p’]
%
% -------------------------------------------------------------------------
%                             OUTPUTS
%
% epar           -> Fourier components of the in-plane electric field. See manual for details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% NOTES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% *has been tested extensively only when superstrate and substrate have equal 
%     epsilon and thickness
% 
% *the lengths (a,d) and inverse lengths (k0,kpar) must be set in the same
%  units (e.g. microns and inverse microns, respectively)
%
% *halfnpw = 0 sets the number of harmonic waves to 1, e.g. the scattering
%  matrix reduces to the ordinary 4x4 formalism for unpatterned multilayers
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Building the lattice
npw = (2*halfnpw+1);
nx = -halfnpw:halfnpw; 



for i = 1:npw 
for j = 1:npw
nx_ij(i,j) = (nx(i)-nx(j));
end
end

kx = (2*pi/a)*nx + kparx;
ky =        0*nx + kpary;


% Initializing the work variables
q   = zeros(2*npw,2*npw,L+2);
phi = zeros(2*npw,2*npw,L+2);
A   = zeros(2*npw,2*npw,L+2);

% Eigensolutions for super- and substrate
q(1:npw,1:npw,1)     = diag(sqrt_whittaker(epssup*k0^2 - kx.*kx - ky.*ky));
q(npw+1:2*npw,npw+1:2*npw,1) = q(1:npw,1:npw,1);
phi(:,:,1) = eye(2*npw);
eps = epssup;

kstorto = [(diag(ky)/eps)*diag(ky), -(diag(ky)/eps)*diag(kx);...
          -(diag(kx)/eps)*diag(ky), (diag(kx)/eps)*diag(kx)];
A(:,:,1)   = ((eye(2*npw)*k0^2 - kstorto)*phi(:,:,1))/q(:,:,1);

q(1:npw,1:npw,L+2)     = diag(sqrt_whittaker(epssub*k0^2 - kx.*kx - ky.*ky));
q(npw+1:2*npw,npw+1:2*npw,L+2) = q(1:npw,1:npw,L+2);
phi(:,:,L+2) = eye(2*npw);
eps = epssub;

kstorto = [(diag(ky)/eps)*diag(ky), -(diag(ky)/eps)*diag(kx);...
          -(diag(kx)/eps)*diag(ky), (diag(kx)/eps)*diag(kx)];
A(:,:,L+2)   = ((eye(2*npw)*k0^2 - kstorto)*phi(:,:,L+2))/q(:,:,L+2);


% Eigensolutions for internal layers
if L > 0
for l = 1:L

	
	F = f(l)*sinc(f(l)*nx_ij);
   
	epsxx = (epsxxB(l) - epsA(l))*F + epsA(l)*eye(npw);
	epsxy =  epsxyB(l)*F; 
	epsyy = (epsyyB(l) - epsA(l))*F + epsA(l)*eye(npw);
	epsz  = (epszB(l)  - epsA(l))*F + epsA(l)*eye(npw);
   
   % Eigensolution calculation
   kstorto = [(diag(ky)/epsz)*diag(ky), -(diag(ky)/epsz)*diag(kx);...
          -(diag(kx)/epsz)*diag(ky), (diag(kx)/epsz)*diag(kx)];
   kdritto = [diag(kx)*diag(kx), diag(kx)*diag(ky);...
              diag(ky)*diag(kx), diag(ky)*diag(ky)];
   epsilone = [epsyy, -epsxy; -epsxy, epsxx];
   [phhi,qq] = eig(epsilone*(eye(2*npw)*k0^2 - kstorto) ...
                    - kdritto);
   q(:,:,l+1) = diag(sqrt_whittaker(diag(qq)));
   phi(:,:,l+1) = phhi;
   A(:,:,l+1) = (eye(2*npw)*k0^2 - kstorto)*phhi/q(:,:,l+1);
   
end
end

% Forward propagation 
S1 = zeros(2*npw,2*npw,L+2); S2 = zeros(2*npw,2*npw,L+2);
S1(:,:,1) = eye(2*npw);
for l = 1:L+1
[S1(:,:,l+1),S2(:,:,l+1)] = smpropag_fw_cond(S1(:,:,l),S2(:,:,l),...
                                phi(:,:,l),phi(:,:,l+1),...
                                A(:,:,l),A(:,:,l+1),...
                                exp(1i*diag(q(:,:,l)*d(l))),exp(1i*diag(q(:,:,l+1)*d(l+1))),...
								sigma(l)*376.730/k0);
end

% Backward propagation
S3 = zeros(2*npw,2*npw,L+2); S4 = zeros(2*npw,2*npw,L+2);
S4(:,:,L+2) = eye(2*npw);
for ll = 1:L+1
    l = L+3-ll;
[S3(:,:,l-1),S4(:,:,l-1)] = smpropag_bw_cond(S3(:,:,l),S4(:,:,l),...
                                phi(:,:,l),phi(:,:,l-1),...
                                A(:,:,l),A(:,:,l-1),...
                                exp(1i*diag(q(:,:,l)*d(l))),exp(1i*diag(q(:,:,l-1)*d(l-1))),...
								sigma(l-1)*376.730/k0);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






% incidence on the 0 order from superstrate % 
as = zeros(2*npw,L+2); bs = zeros(2*npw,L+2);

M = [kparx*kpary,  epssup*k0^2 - kparx^2; ...
    -epssup*k0^2 + kpary^2, -kparx*kpary]/(epssup*q(halfnpw+1,halfnpw+1,1));
if (pol == 's')
abar = linsolve(M,[-kpary; kparx]); %          s-pol
elseif (pol == 'p')
abar = linsolve(M,[kparx; kpary]); %          p-pol
else
error('RTA:PolUnknown','Invalid value for pol')
end
Exyz = vertcat(M*abar,[kpary, -kparx]*abar/epssup); 
norm_E = (Exyz')*Exyz;
abar = abar/sqrt(norm_E);

as(halfnpw+1,1)     = abar(1)*exp(-1i*q(halfnpw+1,halfnpw+1,1)*d(1)); % inizialization of a vectors
as(npw+halfnpw+1,1) = abar(2)*exp(-1i*q(halfnpw+1,halfnpw+1,1)*d(1)); %

clear abar Exyz M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% propagation with waves from superstrate 
for l = 1:L+2
    as(:,l) = (eye(2*npw) - S2(:,:,l)*S3(:,:,l))\(S1(:,:,l)*as(:,1));
    bs(:,l) = (eye(2*npw) - S3(:,:,l)*S2(:,:,l))\(S3(:,:,l)*S1(:,:,l)*as(:,1));
end

% epar nel substrato
l = L+2;
epar =   A(:,:,l)*(as(:,l) - diag(exp(1i*diag(q(:,:,l)*d(l))))*bs(:,l));
