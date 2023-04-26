close all;
clear all
          
% Calculates and plots the joint absorption spectrum  
% in Fig. 4e of "Universal lineshapes at the crossover between weak
%                and strong critical coupling in Fano-resonant
%                coupled oscillators" Sci. Rep. 2016; 6: 24592
%
% Simone Zanotto, Firenze, feb. 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is free software distributed under the BSD licence (see the 
%  containing folder).
% However, shall the results obtained through this code be included 
%  in an academic publication, we kindly ask you to cite the source 
%  website and, if applicable, the following paper:
%
% Simone Zanotto, Alessandro Tredicucci, "Universal lineshapes 
%  at the crossover between weak and strong critical coupling 
%  in Fano-resonant coupled oscillators", Sci. Rep. 2016; 6: 24592
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('PPML_root')    % to be replaced by the proper path 

numom = 201;                  % number of points in energy
om = linspace(100,200,numom); % photon energy (hbar omega) in meV
k0 = om*(2*pi/1240);          % vacuum wavenumber in inverse micron
theta = 0.1;                  % inc. angle. NEVER SET ZERO.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Some structure parameters
a = 3.2;  % lattice spacing

halfnpw = 10; % plane wave truncation. 
              % WARNING! this is a computation-dependent parameter. 
              % inaccurate simulation may result from wrong settings

% material properties              
epsAu = -4000 + 300*1i; % Gold at 10 um

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MQW parameters
%%%% (see Appl. Phys. Lett. 105, 081105 (2014) )
epsbg = 10.05;    % GaAs at 10 um
n2deg = 3e11;     % doping in cm-2
omp   = sqrt(4*pi*4.8e-10^2*n2deg/(epsbg*0.067*9.1e-28*21.8e-7))*6.58e-13; %
                    % eff. MQW bulk plasma freq, meV. :: 21.8 = Lb + Lw (nm) 
om12    = 150;             % meV
gamma12 = 3;               % meV HWHM 
epszn = epsbg./(1-omp^2./(om12^2-om.^2-1i*2*gamma12*om)); % see 


% main cycle over photon energy
for i = 1:size(om,2)
    i

epssup = 1; 
epssub = 1;   % super- and sub-strate permeabilities (air)
L = 3;        % number of internal layers
    
    %           GaAs   | metal pett. | GaAs
f     = [       0        0.73          0            ]; % internal layer fill fraction  (fraction of material B in A)
d     = [0      0.5      0.05          1.3        0 ]; % layer thicknesses
epsxA = [       epsbg    epsbg         epsbg        ]; 
epszA = [       epsbg    epsbg         epszn(i)     ];   
epsxB = [       epsbg    epsAu         epsbg        ];     
epszB = [       epsbg    epsAu         epszn(i)     ];   
sigma = [   0          0         0             0    ];

% calling the PPML SM_1d_tm function
[rl,rr,tlr,trl] = SM_1d_tm(a,L,...
   epssup,epssub,epsxA,epszA,epsxB,epszB,sigma,f,d,...
   halfnpw,k0(i),k0(i)*sin(theta*pi/180)); 
t_vec(i) = tlr;
rl_vec(i) = rl;
rr_vec(i) = rr;
end

%%% Calculating the coherent absorption parameters
%%%% see "Interferometric control of absorption in
%             thin plasmonic metamaterials: general
%             two port theory and broadband operation"
%       Optics Express 23, 9202 (2015); DOI:10.1364/OE.23.009202
T =  abs(t_vec).^2;
Rl = abs(rl_vec).^2;
Rr = abs(rr_vec).^2;

Al = 1-T-Rl;
Ar = 1-T-Rr;

detSmq = abs(rl_vec.*rr_vec - t_vec.^2).^2;
Am = sqrt((1-Al).*(1-Ar)-detSmq);
xi = 0.7;
Ajoint_max_1 = 0.5*(1+xi)*Al + 0.5*(1-xi)*Ar + sqrt(1-xi^2)*Am;
Ajoint_min_1 = 0.5*(1+xi)*Al + 0.5*(1-xi)*Ar - sqrt(1-xi^2)*Am;


plot(om,Ajoint_max_1,'k.')
hold on;
plot(om,Ajoint_min_1,'k.')
