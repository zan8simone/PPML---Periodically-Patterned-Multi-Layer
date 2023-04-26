function [F, F1, F2, F3, F4, G1, G2, G13, G24] = fill_new_F(halfnpw,npw,nx,ny,f1,f2,f3,f4)



   F = zeros(npw);
   
% Quite crazy ordinary Fourier transforms of an L shaped inclusion    
   for i = 1:npw 
   for j = 1:npw
   F(i,j) = (f1-f3)*f4*sinc((f1-f3)*(nx(i)-nx(j)))*sinc(f4*(ny(i)-ny(j))) ...
                 *exp(1i*pi*(f1+f3-1)*(nx(i)-nx(j)))*exp(1i*pi*(f4-1)*(ny(i)-ny(j))) ... 
           + (f2-f4)*f3*sinc(f3*(nx(i)-nx(j)))*sinc((f2-f4)*(ny(i)-ny(j))) ...
                 *exp(1i*pi*(f3-1)*(nx(i)-nx(j)))*exp(1i*pi*(f2+f4-1)*(ny(i)-ny(j))) ... 
           + f3*f4*sinc(f3*(nx(i)-nx(j)))*sinc(f4*(ny(i)-ny(j))) ...
                 *exp(1i*pi*(f3-1)*(nx(i)-nx(j)))*exp(1i*pi*(f4-1)*(ny(i)-ny(j))); 
   
   end
   end
   
   
   % Crazy crossed mixed Fourier transforms of an L shaped inclusion   
   F1 = zeros(2*halfnpw+1);
   F2 = zeros(2*halfnpw+1);
   F3 = zeros(2*halfnpw+1);
   F4 = zeros(2*halfnpw+1);
   G1 = zeros(2*halfnpw+1);
   G2 = zeros(2*halfnpw+1);
   G13 = zeros(2*halfnpw+1);
   G24 = zeros(2*halfnpw+1);
   for i = 1:2*halfnpw+1
   for j = 1:2*halfnpw+1
   F1(i,j) = f1*sinc(f1*(i-j))*exp(1i*pi*(f1-1)*(i-j));
   F2(i,j) = f2*sinc(f2*(i-j))*exp(1i*pi*(f2-1)*(i-j));
   F3(i,j) = f3*sinc(f3*(i-j))*exp(1i*pi*(f3-1)*(i-j));
   F4(i,j) = f4*sinc(f4*(i-j))*exp(1i*pi*(f4-1)*(i-j));
   
   G1(i,j) = (1-f1)*sinc((1-f1)*(i-j))*exp(1i*pi*f1*(i-j));
   G2(i,j) = (1-f2)*sinc((1-f2)*(i-j))*exp(1i*pi*f2*(i-j));

   G13(i,j) = (f1-f3)*sinc((f1-f3)*(i-j))*exp(1i*pi*(f1+f3-1)*(i-j));
   G24(i,j) = (f2-f4)*sinc((f2-f4)*(i-j))*exp(1i*pi*(f2+f4-1)*(i-j));
   end
   end
   
   
   
   end