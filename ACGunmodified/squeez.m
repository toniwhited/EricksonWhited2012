% End of function



% -------------------------------------------------------
% Subroutine to compute squeezes.
%
%  This subroutine compares the values of the objective function at the
%  points c+s*dc and c+0.5*s*dc, with dc = the proposed change in the vector
%  of parameters, and with step length s initially at 1. s is halved until
%  minus the objective function stops declining.
%
function[s_c,s_iter] = squeez(s_c,s_dc,emom,maxsqez,obj,estim,neq,w);

    s_c1=s_c - s_dc;  s_lm=1/2;  s_itr=1;
      s_f1 = deff(s_c1,emom,estim,neq,0);

    lc1 = s_f1'*w*s_f1;
    s_f1 = 0;
    iic  = 0;
  while s_itr <= maxsqez;

    s_c2=s_c-s_lm*s_dc;
      s_f2 = deff(s_c2,emom,estim,neq,0);
    lc2 = s_f2'*w*s_f2;
    clear s_f2;

    if lc1 <= lc2 & lc1 <= obj;

       s_c = s_c1;  s_iter = s_itr-1;  return;

    else;

       s_c1=s_c2;  s_lm=s_lm/2;  lc1=lc2;
       s_itr=s_itr+1;

    end;

  end;
s_c = s_c2;  s_iter = s_itr-1;
return;
