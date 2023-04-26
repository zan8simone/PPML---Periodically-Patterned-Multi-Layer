function [S3r,S4r] = smpropag_bw_cond(S3,S4,p1,p2,A1,A2,f1,f2,sigma)

% Backward propagation of scattering matrix. 
% Includes the effect of a conducting interface
% ----- railway Novara-Pisa, 20 feb 17 -------

J1 = ( A1\A2 - sigma*(p1\A2) + p1\p2 )/2;
J2 = (-A1\A2 + sigma*(p1\A2) + p1\p2 )/2;
J3 = (-A1\A2 - sigma*(p1\A2) + p1\p2 )/2;
J4 = ( A1\A2 + sigma*(p1\A2) + p1\p2 )/2;

S3r = (J4 - diag(f1)*S3*J2)\(diag(f1)*S3*J1 - J3)*diag(f2);
S4r = (J4 - diag(f1)*S3*J2)\(diag(f1)*S4);

end