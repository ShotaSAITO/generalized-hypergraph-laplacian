[H,Y] = cancer;
k = max(Y);
K=2;
EigVecs = zeros(size(H,1),k);
[A,Dn,De,w,L] = HypLap(H,'Saito');
[EigVecsInit,Lambda] = eig(L);
EigVecsInit = EigVecsInit(:,1:K);
ntrials = 10;
ps = 1.1:0.1:3;
accs = zeros(ntrials,length(ps));


for p_i = 1:length(ps)
    p = ps(p_i);
    %EigVecsInit = randorth(size(H,1),K);
    r = @(psi)diff_JC(psi,p,H,De,Dn,w,A,XiInv(EigVecsInit(:,1),p));
    [EigVecs,EigVals] = EmbeddingAlgo(r,EigVecsInit,K,p);
    id=find(EigVals>=min(EigVals),2,'first');
    ix = id(end); % is the index of the smallest m-th element in X
    EigVec = EigVecs(:,ix);
    for trial = 1:ntrials
        a = kmeans(EigVec,k);
        acc = AccMeasure(a,Y)
        accs(trial,p_i) = acc;
    end
end
