function q = sqrt_whittaker(qq)
% Square root of a complex number, 
% choosing the solution with positive imaginary part.

for i = 1:size(qq,1)
for j = 1:size(qq,2)
q(i,j) = sqrt(qq(i,j));
    
if(real(q(i,j))>=0)
   if(imag(q(i,j))<0)
   q(i,j) = -q(i,j);
   end
else
   if(imag(q(i,j))<=0)
   q(i,j) = -q(i,j);
   end 
end

end
end

