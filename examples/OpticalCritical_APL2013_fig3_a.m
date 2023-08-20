close all;
clear; clc

% Calculates and plots the angularly-resolved reflectance
% theoretical counterpart of Fig. 3a of "Optical critical coupling into
%                highly confining metal-insulator-metal
%                optical resonators" Appl. Phys. Lett. 103, 091110, (2013)
%
% Simone Zanotto, Firenze, feb. 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is free software distributed under the BSD licence (see the
%  containing folder).
% However, shall the results obtained through this code be included
%  in an academic publication, we kindly ask you to cite the source
%  website and, if applicable, the following paper:
%
% J.-M- Manceau, S. Zanotto, I. Sagnes, G. Beaudoin, and R. Colombelli,
%  "Optical critical coupling into highly confining metal-insulator-metal
%   optical resonators" Appl. Phys. Lett. 103, 091110, (2015)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('PPML_root')    % to be replaced by the proper path

% Computation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numnu = 41;   % Number of frequency points
nu = linspace(2,4.5,numnu);      % frequencies in THz
numtheta = 55;
theta = linspace(13,67,numtheta);     % angle in degrees (never set 0)

halfnpw = 20;         % good for convergence in this problem
% NOTE: this parameter is computation-dependent.
% Be careful in drawing conclusions!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = zeros(numnu, numtheta);
for i = 1:numnu
    for j = 1:numtheta
        % main cycles over energy and angles
        fprintf('Frequency: %i, Theta: %i.\n', i, j) 

        % dispersive dielectric constants
        % mind the signs of imag parts (misprint in the APL)
        wn = nu(i)*100/3;        % wavenumber in cm^-1 !!!!!!
        epsAu = 1 - (7.25e4)^2/( wn^2 + 1i*216*wn ); % Ordal Au
        epsGaAs = 11*(1 + (292^2 - 268^2)/(268^2 - wn^2 - 1i*2.4*wn)); % GaAs with phonon
        om = 2*pi*nu(i)*1e12;       % ang freq in Hz
        epsGaAs_dop = 11 - (2.6e24*1.6e-19^2/(8.8e-12*0.067*9.1e-31*(om^2 + 1i*om/100e-15))); % doped GaAs

        a = 30;  % microns
        L = 4;   % number of internal layers
        h = 4.67;   % etching depth

        %        superstr.  |  patt. Au  |  patt. GaAs  |  unp. GaAs  | dop. GaAs | subst.
        f     = [              .77          .77            1            1                       ]; % fraction of B in A
        d     = [    30        .2           h              (9-h)        1                 30    ]; % microns

        epsxA = [              1            1              1            1                       ]; % material A
        epszA = [              1            1              1            1                       ]; %
        epsxB = [              epsAu        epsGaAs        epsGaAs      epsGaAs_dop             ]; % material B
        epszB = [              epsAu        epsGaAs        epsGaAs      epsGaAs_dop             ]; %
        sigma = [          0            0              0            0                0          ];

        epssup = 1;
        epssub = epsAu;

        k0   = 2*pi*nu(i)*0.01/3;         % wavevector in micron ^-1
        kpar = k0*sin(theta(j)*pi/180);

        % calling the scattering matrix subroutine for reflection
        [RR,TT,AA] = RTA_1d_tm(a,L,...
            epssup,epssub,epsxA,epszA,epsxB,epszB,sigma,f,d,...
            halfnpw,k0,kpar);
        R(i,j) = RR;

    end
end
%%

h = pcolor(theta,nu,R);
set(h,'edgecolor','none')
xlabel('incidence angle (deg)')
ylabel('frequency (THz)')