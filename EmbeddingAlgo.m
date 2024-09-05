function  [EigVecs, EigVals] = EmbeddingAlgo(diff_J,Psi,k,p)

%set initial diff
Psi_out = Psi*0;
epsilon = abs(sum(sum(Psi_out(:,2:k) - Psi(:,2:k)))/sum(sum(Psi(:,2:k))));

%set initial condition
EigVecs = Psi;

Psi = Psi(:,2:k);
Psi_out = Xi(Psi,p)./vecnorm(Xi(Psi,p));

%Psi_out = XiInv(Psi,p)./vecnorm(XiInv(Psi,p));


J = Psi*0;
EigVals = zeros(1,k);
EigValsSum_out = sum(EigVals);
EigValsSum = 100;%%Initialize
epsilon_obj = abs(abs(EigValsSum_out - EigValsSum)/EigValsSum_out);
iter = 0;
while epsilon > 0.00001 && epsilon_obj > 0.00001 && iter< 1000
    iter = iter  + 1;
    Psi = Psi_out;
    EigValsSum = EigValsSum_out;
    for i = 1:k-1
        psi_tmp = Psi(:,i);
        psi = psi_tmp;
        [j,lambda] = diff_J(psi);
        EigVals(1,i+1) = lambda;
        J(:,i) = j; 
    end
    G = J - Psi * J'*Psi;
    alpha = 0.002 * sum(sum(abs(Psi)))/sum(sum(abs(G)));
    Psi_out = Psi - alpha*G;
    EigValsSum_out = sum(EigVals);
    epsilon_obj = abs(abs(EigValsSum_out - EigValsSum)/EigValsSum_out);
    Psi_out = Psi_out./vecnorm(Psi_out,p);
    epsilon = abs(sum(sum(abs(Psi_out - Psi)))/sum(sum(abs(Psi))));
    fprintf('%i\n', epsilon)
    fprintf('%i, epsilon, %i\n', sum(EigVals),epsilon_obj)
    fprintf('%i, iter\n',iter)
end

Psi_out = XiInv(Psi_out,p)./vecnorm(XiInv(Psi_out,p));
%Psi_out = Xi(Psi_out,p)./vecnorm(Xi(Psi_out,p));

EigVecs(:,2:k) = Psi_out;