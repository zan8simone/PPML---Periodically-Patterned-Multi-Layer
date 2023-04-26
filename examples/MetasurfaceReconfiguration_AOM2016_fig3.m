close all;
clear all; clc
 
% Calculates and plots the handedness-preserving reflectance spectra 
%   at x = 0 of Fig. 3f in "Metasurface reconfiguration through lithium ion 
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

numlambda = 21;                              % Number of lambda points
lambda = linspace(1200,1700,numlambda);      % lambda in nm 

theta = 15;     % polar angle in degrees (never set 0)
phi   = 0;      % azimuthal angle in degrees 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Some structure parameters
a = 820;      % lattice spacing

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
for i = 1:numlambda
i

%           air  |  antennas | sw oxide  |   Pt           
f1    = [           0.76          1                 ]; % fraction of B in A
f2    = [           0.60          1                 ]; % fraction of B in A
f3    = [           0.30          1                 ]; % fraction of B in A
d     = [1500       150           480           200 ]; % nm
epsA  = [           1             1                 ]; % material A 
epsB  = [           epsAl         epsOxide          ]; % material B
sigma = [      0             0           0          ];

k0   = 2*pi/lambda(i);         % wavevector in nm ^-1
kparx = k0*sin(theta*pi/180)*cos(phi*pi/180);
kpary = k0*sin(theta*pi/180)*sin(phi*pi/180);

% calling PPML function "rxx_2d_Lshape"
S = ZSM_2d_Lshape(a,L,...
   epssup,epssub,epsA,epsB,sigma,f1,f2,f3,f3,d,...
   halfnpw,k0,kparx,kpary);
R{i} = S(1:2,1:2);
end

clearvars -except R lambda

%%


Lambda = (1/sqrt(2))*[1,  1; -1i, 1i];

for i = 1:length(lambda)
Rcirc  = (Lambda')*R{i}*Lambda;
RRR(i) = abs( Rcirc(1,1) ).^2 ;
RRL(i) = abs( Rcirc(2,1) ).^2 ;
RLR(i) = abs( Rcirc(1,2) ).^2 ;
RLL(i) = abs( Rcirc(2,2) ).^2 ;
end
Rtot = RRR+RRL+RLR+RLL;

plot(lambda,(RRR+RRL-RLR-RLL)./Rtot,'b'); hold on;
plot(lambda,(RRR+RLL)./Rtot,'r'); hold on;
