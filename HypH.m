function [H] = HypH(X)

%%%%%
%   Usage [H] = HypH(X)
% 
%   Input: X
%     X is a features x node matrix
%     This is returning features into making hyperedge
%
%%%%%
tmp = 0;
for i = 1:size(X,1)
    tmp = tmp + max(X(i,:));
end

H = zeros(size(X,2),tmp);
tmp = 0;
for i = 1:size(X,1)
    %jj = max(X(i,:));
    jj = unique(X(i,:));
    for j = 1:size(jj,2)
        tmp = tmp + 1;
        H(X(i,:) == jj(j),tmp) = 1;
    end
end

%%to prevent isolated nodes
a = sum(H,1) <= 1;
H(:,a) = [];

end