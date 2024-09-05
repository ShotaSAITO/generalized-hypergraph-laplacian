function x = XiInv(psi,p, is_inf)
%x = abs(psi).^(1-p).*sign(psi);
% if ~exist('is_inf', 'var') 
%     is_inf = true; 
% end
% if is_inf == true
%     x(abs(x)>10^10) = Inf;
% end
x = abs(psi).^(1/(p-1)).*sign(psi);
if p == 1
   p = 1.00001;
   x = abs(psi).^(1/(p-1)).*sign(psi); 
end
end