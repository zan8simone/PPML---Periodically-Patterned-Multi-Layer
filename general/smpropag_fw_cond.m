function [S1r,S2r] = smpropag_fw_cond(S1,S2,p1,p2,A1,A2,f1,f2,sigma)

% Forward propagation of scattering matrix. 
% Includes the effect of a conducting interface
% ----- railway Novara-Pisa, 20 feb 17 -------

I1 = ( A1\A2 + sigma*(p1\A2) + p1\p2 )/2;
I2 = (-A1\A2 - sigma*(p1\A2) + p1\p2 )/2;
I3 = (-A1\A2 + sigma*(p1\A2) + p1\p2 )/2;
I4 = ( A1\A2 - sigma*(p1\A2) + p1\p2 )/2;

S1r = (I1 - diag(f1)*S2*I3)\(diag(f1)*S1);
S2r = (I1 - diag(f1)*S2*I3)\(diag(f1)*S2*I4 - I2)*diag(f2);

end