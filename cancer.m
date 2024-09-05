function [H,Y] = cancer
addpath('../lib')
[X,Y] = cancer_dataset;
Y = Y(1,:)'+1;
H = HypH(X);
end