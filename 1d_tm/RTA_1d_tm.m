%%%
function [RR,TT,AA] = RTA_1d_tm(a,L,...
    epssup,epssub,epsxA,epszA,epsxB,epszB,sigma,f,d,...
    halfnpw,k0,kpar)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is free software distributed under the BSD licence (see the
%  containing folder).
% However, shall the results obtained through this code be included
%  in an academic publication, we kindly ask you to cite the source
%  website.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates reflection, transmission and absorption.
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
% f              -> duty cycle in the internal layers
%                                         [vector with L components]
% d              -> thicknesses of the layers,
%                   including super- and sub-strate
%                                         [vector with (L+2) components]
% halfnpw        -> half number of harmonic waves used in calculation
% k0             -> wavevector in vacuum  (= 2*pi/lambda0)
% kpar           -> x-projection of the incident wavevector
%
% -------------------------------------------------------------------------
%                             OUTPUTS
%
% RR             -> Reflectance      (reflected flux/incident flux)
% TT             -> Transmittance    (transmitted flux/incident flux)
% AA             -> Absorbance for each internal layer & adjacent conducting interface
%                                    (absorbed power/incident flux)
%                                              [vector with L components]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% NOTES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% *the lengths (a,d) and inverse lengths (k0,kpar) must be set in the same
%  units (e.g. microns and inverse microns, respectively)
%
% *epssup must be real (otherwise the incident waves are ill-defined)
%
% *epssub can be either real or complex.
%  In the first case, TT is the transmittance towards the far field.
%  In the second case, TT is the absorbance in the substrate...
%
% *Above diffraction thresholds, RR and TT contain also the contribution of
% diffracted beams
%
% *halfnpw = 0 sets the number of harmonic waves to 1, e.g. the scattering
%  matrix reduces to the ordinary 2x2 formalism for unpatterned multilayers
%
% *thickness of superstrate and substrate do not influence the result.
%  in effect, super- and substrate are semi-infinite.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n = -halfnpw:halfnpw;
npw = size(n,2);
kx = (2*pi/a)*n + kpar;

% Initializing the work variables
q   = zeros(npw,npw,L+2);
phi = zeros(npw,npw,L+2);
A   = zeros(npw,npw,L+2);
epsx = zeros(npw,npw,L+2);
etaz = zeros(npw,npw,L+2);

% Eigensolutions for super- and substrate
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

% Eigensolutions for internal layers
if L > 0
    for l = 1:L

        F = zeros(npw);

        for i = 1:npw % Fourier transform of square centered at the origin
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

        % Eigensolution calculation
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

% Flux at the beginning of internal layers + beginning of the substrate
flux = zeros(L+1,1);
for l = 2:L+2
    hy = phi(:,:,l)*(as(:,l) + diag(exp(1i*diag(q(:,:,l)*d(l))))*bs(:,l));
    ex =   A(:,:,l)*(as(:,l) - diag(exp(1i*diag(q(:,:,l)*d(l))))*bs(:,l));

    hy = hy.'; ex = ex.';
    flux(l-1) = real(hy*ex')/incflux;
end
RR = 1-flux(1);
TT = flux(L+1);


AA = zeros(L,1);
if L>0
    for l = 1:L
        AA(l) = flux(l)-flux(l+1);
    end
end

if abs(sum(AA)+RR+TT-1) > 1e-5 % Check energy conservation
    error('RTA:EnNotCons','Energy not conserved')
end
