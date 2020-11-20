function [ Vec ] = NestR(R,M,N )
% Sort non-repeated rows via two-level nested array covariance matrix R
% R: M+N dimension 
% M: inner element no.
% N: outer element no.
% D: sensor spacing

Y = reshape(R,(M+N)^2,1);
d = 0.5;
Sp1 = [0:M-1]*d;    % inner sensor location
Sp2 = [M*d:(M+1)*d:(N*(M+1)-1)*d];  % outer sensor location
Sp = [Sp1,Sp2];
co = round(angle(kron(exp(-j*Sp'/100),exp(j*Sp'/100)))*2*100);
Nc = 2*N*(M+1) - 1;
p = 0;
for count = -((Nc+1)/2 - 1):((Nc+1)/2 - 1)
    p=p+1;
    temp = find(co == count);
    temp1 = temp(1);
    Z(p) = Y(temp1);
end
Vec = Z;


end

