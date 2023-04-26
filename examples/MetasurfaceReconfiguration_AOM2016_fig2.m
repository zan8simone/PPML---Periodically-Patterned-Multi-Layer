close all;
clear all; clc

% Calculates and plots the polarization-resolved reflectance spectra 
%   at x = 0 of Fig. 2d in "Metasurface reconfiguration through lithium ion 
%   intercalation in a transition metal oxide", 
%   Advanced Optical Materials, 2016
%
% Simone Zanotto, Firenze, nov. 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is free software distributed under the BSD licence (see the 
%  containing folder).
% However, shall the results obtained through this code be included 
%  in an academic publication, we kindly ask you to cite the source 
%  website and, if applicable, the following paper:
%
% Simone Zanotto et al., "Metasurface reconfiguration through lithium ion 
%   intercalation in a transition metal oxide", 
%   Advanced Optical Materials", (2016)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('PPML_root')    % to be replaced by the proper path 

numlambda = 51;                              % Number of lambda points
lambda = linspace(1200,1700,numlambda);      % lambda in nm 

theta = 15;     % polar angle in degrees (never set 0)
phi   = 0;      % azimuthal angle in degrees 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Some structure parameters
a = 850;      % lattice spacing

halfnpw = 10; % plane wave truncation. 
              % WARNING! this is a computation-dependent parameter. 
              % inaccurate simulation may result from wrong settings

% material properties              
epsOxide = 5;         % approximate non-dispersive values
epsPt = -25  + 1i*70;  %
epsAl = -180 + 1i*35;  %

epssup = 1; 
epssub = epsPt;   % super- and sub-strate permeabilities 
L = 2;            % number of internal layers (patterned Al, switching oxide)
    
% main cycle over wavelength
parfor i = 1:numlambda
i

%           air  |  antennas | sw oxide  |   Pt   
fx    = [           510/850,       1                 ]; % values from SEM measurement
fy    = [           270/850,       1                 ]; % values from SEM measurement
d     = [1500       150,         500           200   ]; % 
epsA  = [           1            1                   ]; % material A 
epsB  = [           epsAl       epsOxide             ]; % material B
sigma = [        0           0           0           ];

k0   = 2*pi/lambda(i);         % wavevector in nm ^-1
kparx = k0*sin(theta*pi/180)*cos(phi*pi/180);
kpary = k0*sin(theta*pi/180)*sin(phi*pi/180);

S = ZSM_2d_rect(a,a,L,...
   epssup,epssub,epsA,epsB,sigma,fx,fy,d,...
   halfnpw,k0,kparx,kpary);
rss(i)   = S(1,1);
rsp(i)   = S(2,1);
rps(i)   = S(1,2);
rpp(i)   = S(2,2);
end

%%

plot(lambda,abs(rss).^2,'k','LineWidth',0.5)
hold on
plot(lambda,abs(rpp).^2,'k','LineWidth',2)