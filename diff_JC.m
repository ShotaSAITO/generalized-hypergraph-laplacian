function [diff_J,lambda] = diff_JC(psi,p,H,De,Dn,w,A,psi_1)

Dn_root = Dn^(-1/2);

%%Proposed
De = sparse(De - w);
%De_root = sparse(diag(  (sum(H,1) - 1).^ (-1/2)  ));
%W = sparse(De^(-1) * w);
%A = H*De^(-1)*w*H';
%A = A-diag(diag(A));


[~,c] = ExtSubSp(Xi(psi,p),Xi(psi_1,p),p);

%%computing rayleigh quotient
delta_original = H' * Dn_root * psi;
u = delta_original(:,ones(1,size(H',2)));
YY = (Dn_root * psi)';
u_ = (De + eye(size(De,1))) * YY(ones(1,size(H',1)),:);

%fDe = diag(factorial(diag((De + eye(size(De,1))))));

delta = (De + w)^(-1) * De^(-1) * (w * H'.*abs(u - u_)).^p;


LL_left_original = sum((De + w)^(-1)*De^(-1) * w * H'.* abs(u - u_).^(p-1),2);
LL_left = LL_left_original(:,ones(1,size(H',2)));
LL_right = De^(-1) * w * H'.* abs(u - u_).^(p-1);
%LL_right = (De - eye(size(De,1),size(De,2))) * De^(-1) * w * H'.* (u - u_).^(p-1) * Dn_root;

LL = Dn_root*sum((LL_left - LL_right))';

s_p = sum(sum(delta));
lambda = s_p/norm(psi,p)^p;

n = norm(psi - c*psi_1,p);

diff_J = (1/n^p - s_p/n^(2*p)).*(LL - abs(psi - c*psi_1).^(p-1).*abs(psi).^(-p));

%nabla = sparse(De^(-1) * H' * delta_v_p);
%nabla = H*(De + eye(size(De,1)) )^(-1)*(w.* nabla)*H' - diag(diag(H*(De + eye(size(De,1)) )^(-1)*(w.* nabla)*H'));

%Lp = Dn_root * ( - H*De^(-1)*(w.* nabla)*H' + diag(delta_v_p) * A + A * diag(delta_v_p) - diag(diag(A)) * diag(delta_v_p)) * Dn_root;
%Wp =  - H*De^(-1)*(w.* nabla)*H' + sparse(diag(delta_v_p)) * A + A * diag(delta_v_p) - sparse(diag(diag(A))) * diag(delta_v_p);
%Wp =  - nabla +  sparse(diag(delta_v_p)) * A + A * diag(delta_v_p) - sparse(diag(diag(A))) * diag(delta_v_p);
%Dp = diag(sum(Wp));
%Lp = Dn_root * (Dp - Wp) * Dn_root;

%lambda = psi'*Lp*psi/norm(psi,p)^p;

%diff_J = p*(Lp*psi - lambda*Xi(psi,p))/norm(psi,p)^p;


%%computing rayleigh quotient
delta_original = H' * Dn_root * Xi(psi,p);
u = delta_original(:,ones(1,size(H',2)));
YY = (Dn_root * Xi(psi,p))';
u_ = (De + eye(size(De,1))) * YY(ones(1,size(H',1)),:);

%fDe = diag(factorial(diag((De + eye(size(De,1))))));

delta = (De + w)^(-1) * De^(-1) * (w * H'.*abs(u - u_)).^p;
s_p = sum(sum(delta));
lambda = s_p/norm(Xi(psi,p),p)^p;


end