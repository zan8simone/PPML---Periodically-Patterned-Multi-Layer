function [F, F1, F2, F3, F4, G1, G2, G13, G24] = fill_F(halfnpw_r,npw_r,nx_r,ny_r,f1_r,f2_r,f3_r,f4_r)



% *_r : parameters for which a request is done
% F, ... G... : matrices corresponding to the required parameters

% Fdata: .mat containing   (input)                 halfnpw npw nx ny f1 f2 f3 f4    
%                    and   (correspondnig output ) F1 F2 ... ... G1 ... ...
%                                                     


load Fdata



if isequal(halfnpw_r,halfnpw) && isequal(npw_r,npw) && isequal(nx_r,nx) && isequal(ny_r,ny) && isequal(f1_r,f1) && isequal(f2_r,f2) && isequal(f3_r,f3) && isequal(f4_r,f4)
	
    return
	
else

    halfnpw = halfnpw_r; npw = npw_r; nx = nx_r; ny = ny_r; f1 = f1_r; f2 = f2_r; f3 = f3_r; f4 = f4_r; 
    [F, F1, F2, F3, F4, G1, G2, G13, G24] = fill_new_F(halfnpw,npw,nx,ny,f1,f2,f3,f4);
    clearvars -except F F1 F2 F3 F4 G1 G2 G13 G24 halfnpw npw nx ny f1 f2 f3 f4
    save Fdata	
	
end
