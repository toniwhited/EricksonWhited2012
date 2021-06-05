function [u, X, xerror, IVerror, fe] = genvars(nt, Zetadim, mism, mu, sig, Setup)
%FUNCTION:  genvars
%Version: 9/15/2009
%Observations random generation
%INPUTS:    nt - total observarions
%           xZdim - dimension of xi (mismeasured) + Zeta (perfectly measured)
%           mu, sig - parameters for the normal distribution
%           Setup - type of distribution used to generate variables
%CALL:      chosedist
%OUTPUTS:   random vectors u, X, xerror, IVerror, and fe
%IMPORTANT: Generates X and Zetas according to Whited (NORMALIZED)

xZdim = Zetadim+mism;
dist_u = Setup(1);
dist_xerror = Setup(2);
dist_X = Setup(3);
dist_IVerror = Setup(4);
dist_fix = Setup(5);
df = [Setup(6) Setup(7)];


miuX = mu * ones(xZdim, 1); miuxerror = mu * ones(mism, 1);
sigmaX = sig * eye(xZdim); sigmaxerror = sig * eye(mism);    %zero correlation

u = chosedist(dist_u, mu, sig, nt, df);
xerror = chosedist(dist_xerror, miuxerror, sigmaxerror, nt, df);
X = chosedist(dist_X, miuX, sigmaX, nt, df);
IVerror = chosedist(dist_IVerror, miuxerror, sigmaxerror, nt, df);
fe = chosedist(dist_fix, mu, sig, nt, df);
%ENDFUNCTION

function A = chosedist(tipo, media, varianca, tam, df)
% standardize to guarantee mean 0 and variance 1

switch tipo 
    case 1
        A = mvnrnd(media, varianca, tam);
    case 2
        A = exp( mvnrnd(media, varianca, tam) );
    case 3
        A = mvtrnd(varianca, df(1), tam);
    case 4
        A = chi2rnd(df(1),[tam size(varianca,1)]);
    case 5
        A = frnd(df(1), df(2), [tam size(varianca,1)]);
end
a = length(A);  
mean_A = mean(A,1); mean_A = repmat(mean_A, a, 1);
std_A = std(A); std_A = repmat(std_A, a,1); 
A = (A-mean_A)./std_A;