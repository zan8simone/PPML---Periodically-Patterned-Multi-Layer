function [x,z,Ex,Ez,Sz] = field_1d_tm(a,L,...
   epssup,epssub,epsxA,epszA,epsxB,epszB,sigma,f,d,...
   halfnpw,k0,kpar,nx,nz)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is free software distributed under the BSD license (see the 
%  containing folder).
% However, shall the results obtained through this code be included 
%  in an academic publication, we kindly ask you to cite the source 
%  website.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes fields excited by incident plane wave.
%
% 1-d structure, TM polarization.
%
% Simone Zanotto, Orsay Oct. - Dec 2012; Firenze Feb. 2016
% -------------------------------------------------------------------------
%                              INPUTS
%
% a              -> periodicity 
% L              -> number of layers excluded superstrate and substrate
% epssup, epssub -> super-and sub-strate permittivities 
% epsxA          -> in-plane component of dielectric tensor 
%                   for material A in the internal layers
%                                         [vector with L components]
% epsxB          -> idem, for material B  
% epszA, epszB   -> idem, out-of-plane components
% sigma          -> conductivity of interfaces
%                                         [vector with (L+1) components]
% f              -> fraction of material B in the internal layers
%                         ("dut ycycle")  [vector with L components]
% d              -> thicknesses of the layers, 
%                   including super- and sub-strate
%                                         [vector with (L+2) components]
% halfnpw        -> half number of harmonic waves used in calculation
% k0             -> wavevector in vacuum  (= 2*pi/lambda0)
% kpar           -> x-projection of the incident wavevector 
% nx             -> number of points along x for field display
% nz             -> number of points along z for field display, specified
%                   for each internal layer [vector with (L+2) components]
%
% -------------------------------------------------------------------------
%                             OUTPUTS
%
% x, z           -> vectors of x and z points where field is computed
% Ex, Ez         -> E field x- and z-components  [functions of x and z]
% Sz             -> Poynting vector z-component
%                    averaged over the unit cell [function of z] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% NOTES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% *the lengths (a,d) and inverse lengths (k0,kpar) must be set in the same
%  units (e.g. microns and inverse microns, respectively)
%
% *epssup must be real (otherwise the incident waves are ill-defined)
%
% *E fields are complex. 
%  The real, time-dependent fields are Re(E*exp(-i omega t)) [omega = c*k0]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n = -halfnpw:halfnpw;
npw = size(n,2);
kx = (2*pi/a)*n + kpar;

% Inizializziamo le variabili di lavoro
q   = zeros(npw,npw,L+2);
phi = zeros(npw,npw,L+2);
A   = zeros(npw,npw,L+2);
epsx = zeros(npw,npw,L+2);
etaz = zeros(npw,npw,L+2);

% Autosoluzioni per super- e sub-strato
q(:,:,1)   = diag(sqrt_whittaker(epssup*k0^2 - kx.*kx));
phi(:,:,1) = eye(npw);
A(:,:,1)   = diag(k0^2 - kx.*kx/epssup)/q(:,:,1);
epsx(:,:,1) = eye(npw)/epssup;
etaz(:,:,1) = eye(npw)*epssup;
q(:,:,L+2)   = diag(sqrt_whittaker(epssub*k0^2 - kx.*kx));
phi(:,:,L+2) = eye(npw);
A(:,:,L+2)   = diag(k0^2 - kx.*kx/epssub)/q(:,:,L+2);
epsx(:,:,L+2) = eye(npw)/epssub;
etaz(:,:,L+2) = eye(npw)*epssub;

% Autosoluzioni degli strati interni 
if L > 0
for l = 1:L

   F = zeros(npw);
   
            for i = 1:npw
            for j = 1:npw
                if (i == j)
                F(i,j) = f(l);
                else
                F(i,j) = sin(pi*f(l)*(n(i)-n(j)))/(pi*(n(i)-n(j)));
                end
            end
            end
 
   epsx(:,:,l+1) = (1/epsxB(l) - 1/epsxA(l))*F + 1/epsxA(l)*eye(npw);
   etaz(:,:,l+1) = (  epszB(l) -   epszA(l))*F +   epszA(l)*eye(npw);
   
   % Calcoliamo le autosoluzioni
   [phhi,qq] = eig(epsx(:,:,l+1)\(k0^2*eye(npw) - diag(kx)*(etaz(:,:,l+1)\diag(kx))));
   q(:,:,l+1) = diag(sqrt_whittaker(diag(qq)));
   phi(:,:,l+1) = phhi;
   A(:,:,l+1) = (k0^2*eye(npw) - diag(kx)*(etaz(:,:,l+1)\diag(kx)))*phhi/q(:,:,l+1);
   
end
end

% Forward propagation 
S1 = zeros(npw,npw,L+2); S2 = zeros(npw,npw,L+2);
S1(:,:,1) = eye(npw);
for l = 1:L+1
[S1(:,:,l+1),S2(:,:,l+1)] = smpropag_fw_cond(S1(:,:,l),S2(:,:,l),...
                                phi(:,:,l),phi(:,:,l+1),...
                                A(:,:,l),A(:,:,l+1),...
                                exp(1i*diag(q(:,:,l)*d(l))),exp(1i*diag(q(:,:,l+1)*d(l+1))),...
								sigma(l)*376.730/k0);
end

% Backward propagation
S3 = zeros(npw,npw,L+2); S4 = zeros(npw,npw,L+2);
S4(:,:,L+2) = eye(npw);
for ll = 1:L+1
    l = L+3-ll;
[S3(:,:,l-1),S4(:,:,l-1)] = smpropag_bw_cond(S3(:,:,l),S4(:,:,l),...
                                phi(:,:,l),phi(:,:,l-1),...
                                A(:,:,l),A(:,:,l-1),...
                                exp(1i*diag(q(:,:,l)*d(l))),exp(1i*diag(q(:,:,l-1)*d(l-1))),...
								sigma(l-1)*376.730/k0);
end

% incidence on the 0 order from superstrate
as = zeros(npw,L+2); bs = zeros(npw,L+2);

% unit electric field, phase zero at x = 0, z = d(1) %%%%%%%%%%%%%%%%%%%%%%
as(halfnpw+1,1) = exp(-1i*q(halfnpw+1,halfnpw+1,1)*d(1))*sqrt(epssup)/k0; 
% incident flux
incflux = q(halfnpw+1,halfnpw+1,1)/k0^2; 


for l = 1:L+2
    as(:,l) = (eye(npw) - S2(:,:,l)*S3(:,:,l))\(S1(:,:,l)*as(:,1));
    bs(:,l) = (eye(npw) - S3(:,:,l)*S2(:,:,l))\(S3(:,:,l)*S1(:,:,l)*as(:,1));
end



%%%%%%%%%%%

Ez = zeros(nx,sum(nz,2));
Ex = zeros(nx,sum(nz,2));
Sz = zeros(1,sum(nz,2));

x = linspace(-a/2,a/2,nx);

z0 = 0;
nz0 = 0;
z = zeros(0);
for l = 1:L+2 
zz = linspace(0,d(l),nz(l)); % running z in each layer
z = horzcat(z,zz+z0); 

for iz = 1:nz(l)
    hy = phi(:,:,l)*(diag(exp(1i*diag(q(:,:,l))* zz(iz)      ))*as(:,l) + ...
                     diag(exp(1i*diag(q(:,:,l))*(d(l)-zz(iz))))*bs(:,l)); 
    ex =   A(:,:,l)*(diag(exp(1i*diag(q(:,:,l))* zz(iz)      ))*as(:,l) - ...
                     diag(exp(1i*diag(q(:,:,l))*(d(l)-zz(iz))))*bs(:,l));     
    ez = -etaz(:,:,l)\diag(kx)*hy;
    
    for ix = 1:nx
    Ez(ix,iz+nz0) = (exp(1i*kx*x(ix)))*ez;
    Ex(ix,iz+nz0) = (exp(1i*kx*x(ix)))*ex;
    end
    
    Sz(iz+nz0) = real((hy')*ex)/incflux;
end

z0 = z0+d(l); nz0 = nz0 + nz(l);
end


