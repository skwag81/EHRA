function [Riskn] = RISK(pga,F,H)

% Number of data points
n = length(pga);
% Combine (Integrated by Trapezoidal rule)
IP = 0; h = 1e-3;
% Use forward for first point
INT1 = F(1)*abs((interp1(pga,H,pga(1)+h)-interp1(pga,H,pga(1)))/(h));
INT2 = F(2)*abs((interp1(pga,H,pga(2)+h)-interp1(pga,H,pga(2)))/(h));
dpga = pga(2)-pga(1);
IP = IP +dpga*(INT1+INT2)/2;
% Use central difference method for mid-points
for i =2:n-2
    INT1 = F(i)*abs((interp1(pga,H,pga(i)+h)-interp1(pga,H,pga(i)-h))/(2*h));
    INT2 = F(i+1)*abs((interp1(pga,H,pga(i+1)+h)-interp1(pga,H,pga(i+1)-h))/(2*h));
    dpga = pga(i+1)-pga(i);
    IP = IP +dpga*(INT1+INT2)/2;
end
% Use backward for last point
INT1 = F(n-1)*abs((interp1(pga,H,pga(n-1))-interp1(pga,H,pga(n-1)-h))/(h));
INT2 = F(n)*abs((interp1(pga,H,pga(n))-interp1(pga,H,pga(n)-h))/(h));
dpga = pga(n)-pga(n-1);
IP = IP +dpga*(INT1+INT2)/2;
   
Riskn = IP;

return