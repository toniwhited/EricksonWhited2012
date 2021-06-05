function [Xnames, Zetadim] = genNames(Xdim, mism, WD_order)
%FUNCTION:  genNames
%Version:   9/15/2009
%Creates names for the variables used by Whitted estimator
%INPUTS:    Xdim - parameters dimension
%           mism - number of mismeasured variables
%           WD_ORDER - 3+WD_ORDER is the max moment order for Whited
%           estimation (also WD_ORDER gives the # of estimators)
%OUTPUTS:   Xnames - names for variables (Whited procedures)
%           Zetadim - number of perfectly measured variables
%%%%%%%%%
Xnames = cell(1,Xdim);
name_mism = 1:mism;
Xnames{mism+1} = 'Intercep';
if Xdim >= 2
    for i=1:mism
        Xnames{i} = strcat( 'xaux',int2str(i) );
    end
    for i=mism+2:Xdim
        Xnames{i} = strcat('zeta(:,',int2str( i-(mism+1) ),')');
    end
    if Xdim == mism+1 & WD_order
        disp('Model is Y = Constant + xis + u (xis mismeasured are the only variable on Whitted procedure)')
    end
else if Xdim == 1
        disp('Model is Y = Constant + u')
        return;
    end
end
Zetadim = length(Xnames)-(mism+1);
%endfunction
