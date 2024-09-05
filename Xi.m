function x = Xi(psi,p)
x = abs(psi).^(p-1).*sign(psi);
end