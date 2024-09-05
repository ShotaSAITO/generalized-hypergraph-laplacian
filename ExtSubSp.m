function [psi_out,c] = ExtSubSp(psi_in, EigVecs,p)
   fun = @(c)SubSp(psi_in, EigVecs, p, c);
   c0 = rand(size(EigVecs,2),1);
   options = optimoptions(@fminunc,'MaxIterations',500);
   [c,~] = fminunc(fun,c0,options);
   psi_out = psi_in;
   for i = 1:size(EigVecs,2)
         psi_out = psi_out - c(i,1) * EigVecs(:,i);
   end
end


function n = SubSp(psi_in, EigVecs, p, c)
     for i = 1:size(EigVecs,2)
         psi_in = psi_in - c(i,1) * EigVecs(:,i);
     end
     n = norm(psi_in,p);
end