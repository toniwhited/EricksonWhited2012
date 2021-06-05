new ,60000;
ttt = time;
ddd = date;
call sysstate(14,1e-128);
fe = 0;
asset = 1;

@ Step 1: ------Name of the Gauss data set.----------------------------------@

    dataset = "wdset";
    vnames = getname(dataset);

@ Step 2: ------Name of the output file.-------------------------------------@

    if asset == 0;

    if fe==0;
    output file = idtest.out reset;
    else;
    output file = idtestFE.out reset;
    endif;

    else;


    if fe==0;
    output file = idtestA.out reset;
    else;
    output file = idtestFEA.out reset;
    endif;

    endif;


    if asset == 0;
@ Step 3: ------Name of the left-hand-side variable.-------------------------@

    let yname = ik;

@ Step 4: ------Name of the mismeasured variable.----------------------------@

    let xname = q;

@ Step 4: ------Names of the perfectly measured variables.-------------------@

    let zname = intercep cfk;

    else;
@ Step 3: ------Name of the left-hand-side variable.-------------------------@

    let yname = ia;

@ Step 4: ------Name of the mismeasured variable.----------------------------@

    let xname = mb;

@ Step 4: ------Names of the perfectly measured variables.-------------------@

    let zname = intercep cfa;
    endif;
@ Step 5: ------Set the number of years.-------------------------------------@

    nyr = 42;

@ Step 6: ------Set the first year of data.----------------------------------@

    fyear = 1967;

@ Step 7: ------Set the number of observations in the largest cross-section.-@

    nco = 1900;

@ Step 8:-------Set the output mode.-----------------------------------------@

    latex = 1;

@ Step 9:---------Set the maximum number of iterations.----------------------@

    maxiter = 299;

@ Step 10: -----------------Set the limit on the number of squeezes.---------@

    maxsqez=240;

@ Step 11: -----------------Set the tolerance on the minimization routine.---@

    tol=1e-9;

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*                           End of user input.                             */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
    open f1=^dataset;
    alldta = readr(f1,rowsf(f1));
outwidth  256;


@------------Number of perfectly measured regressors.--------------@

    nz = rows(zname)-1;

@---------------Initialize the counters.---------------------------@

   csave = zeros(nyr,2);

@-------Set the number of equations and parameters.----------------@

nc = 2;
neq = 2;

    @ ------- Now we enter the main body of the program. ----------- @
cls;

   if nz > 0;
   allvname = yname|xname|zname[2:rows(zname)];
   else;
   allvname = yname|xname;
   endif;
   indx = indcv(allvname,vnames);
   yindx = indx[1,1];
   xindx = indx[2,1];
   if nz > 0;
   zindx = indx[3:rows(indx),1];
   endif;
   let yrname = fyear;
   yrindx = indcv(yrname,vnames);



   if fe==1;
   newdta = zeros(1,cols(alldta));
   gvks = unique(alldta[.,1],1);
   for ii(1,rows(gvks),1);
       dtai = selif(alldta,alldta[.,1].==gvks[ii,1]);
      if rows(dtai)>1;
         dtai[.,indx] = dtai[.,indx]-meanc(packr(dtai[.,indx]))';
         newdta = newdta|dtai;
      endif;
   endfor;
   alldta = newdta[2:rows(newdta),.];

   endif;


     @alldta = selif(alldta,sumc((alldta[.,indx]./=miss(1,1))').==0);
     @

  year = 1;
  do while year <= nyr;

    @ --------- Read in the data and define the variables. ----------@

      dta=selif(alldta,(alldta[.,yrindx].==(fyear+year-1)));

      n=rows(dta);                n;;

    y = dta[.,yindx];
    x = dta[.,xindx];
    if nz > 0;
    z = dta[.,zindx];
    endif;


    icept=ones(n,1);
    if nz > 0;
    za = z;
    iz = icept~za;
    else;
    iz = icept;
    endif;

    muy = inv( moment(iz,0) )* (iz'y);
    mux = inv( moment(iz,0) )* (iz'x);

    y_d = y - iz*muy;
    x_d = x - iz*mux;


    @ -------u-- Construct the moments for subsequent estimation. ----@
    y2  = y_d^2;
    x2  = x_d^2;
    yx  = x_d.*y_d;
    y2x = (y_d^2).*x_d;
    yx2 = (x_d^2).*y_d;

    Ey2x = meanc(y2x);
    Eyx2 = meanc(yx2);       year;;Ey2x;;Eyx2;


    inEzz   =   invpd(moment(iz,0)/n);
    Ey2_z   =   meanc(y2   .* iz);
    Eyx_z   =   meanc(yx   .* iz);
    Ex2_z   =   meanc(x2   .* iz);
    zy_d    =   iz.*y_d;
    zx_d    =   iz.*x_d;

    @ ------- Initialize the parameters. --------------------------- @


    c = zeros(nc,1);

    c[1,1] = Ey2x;
    c[2,1] = Eyx2;


     @ --------The program creates the weighting matrix. ------------@
     clear w;

  w = optw(c);

     @ ------------- Start of iteration loop ------------------------@
iter = 1;
    tol=1e-9;                 @ Set the convergence criterion. @
    dc=1;                     @ Initialize the step length. @
do until abs(dc) < tol;
    clear g, f, sqz;

    f=deff(c);                @ The program jumps to the subroutine that
                                defines the f vector. @

    g=grad(c);                @ The program jumps to the subroutine that
                                computes analytic derivatives. The matrix
                                of partials of f with respect to the
                                parameters is called g. @

    obj = f'w*f;              @ This computes the value of the objective
                                function.@

    gwg= g'w*g;               @ This uses the GAUSS-NEWTON method
                                to compute full step dc. @
    gwf = g'w*f;


    onemore = 1;
    solveit:
    trap 1;
    dc = solpd(gwf,gwg);
    trap 0;
    if scalerr(dc) > 0;
      print "didn't invert! ";; "iter = ";; iter;
      if onemore == 1;
        call sysstate(14,1e-128);
        onemore = 0;
        goto solveit;
      else;
        print "whoops! ";; "iter = ";; iter;
        end;
      endif;
    endif;


 if maxsqez > 0;              @ This jumps to the subroutines that
                                adjust the step length. @
   { c_new,sqz } = squeeze(c,dc);
 else;
   c_new=c - dc;
 endif;


 dc=c_new-c;                     @ Update variables for the next iteration. @
 c=c_new;

 iter=iter +1;
                                         @ Quit iterating if necessary. @
 if iter >= maxiter; goto escape; endif;
 if key == 27; goto escape; endif;       @ Hitting the escape key will stop
                                           the iteration loop.  @
endo;
@ ------ End of iteration loop ------------------------------------ @
escape:

@ Compute t-ratios, etc.                                      @

  w = optw(c);
  g=grad(c);
  gwg = g'w*g;

clear vc;
vc=invpd(gwg)/n;
stderr=sqrt(diag(vc));
t=c./stderr;
df=n-nc;
pvt=2*cdftc(abs(t),df);

nww = n*c'w*c;
test=cdfchic(nww,2);

csave[year,1] = nww;
csave[year,2] = test;

year = year + 1;
endo;

/*--------------This business makes a LaTeX table------------*/

format /rd 8,3;
if latex == 1;

for jj(1, nyr, 1);
  "& ";; csave[jj,1];; "  ";;  " \\" "\\";
  "&(";; csave[jj,2];; ") ";;  " \\" "\\";
endfor;
meanc(csave[.,2].<0.05);

else;

for jj(1, nyr, 1);
  "  ";; csave[jj,1];; "  ";
  " (";; csave[jj,2];; ") ";
endfor;

endif;

goto eop;  @ This skips all of the subroutines.        @


@ ----------------- Subroutines follow --------------- @

@ Subroutine to define the f vector. @

proc deff(a);
local f;

f = zeros(neq,1);

f[1,1] = Ey2x - a[1,1];
f[2,1] = Eyx2 - a[2,1];

retp(f);
clear f;
endp;

@  ------------------------------------------------------- @
@ Subroutine to compute gradient matrix. @

proc grad(a);

retp(-eye(nc));

endp;

@  ------------------------------------------------------- @
@ Subroutine to compute squeezes.

  This subroutine compares the values of the objective function at the
  points c+s*dc and c+0.5*s*dc, with dc = the proposed change in the vector
  of parameters, and with step length s initially at 1. s is halved until
  minus the objective function stops declining.
@
proc (2) = squeeze(s_c,s_dc);
local s_c1,s_lm,s_itr,lc1,s_c2,lc2,s_f1,s_f2;


    s_c1=s_c - s_dc; s_lm=1/2; s_itr=1;
    s_f1 = deff(s_c1);
    lc1 = s_f1'w*s_f1;
    clear s_f1;

  do until s_itr > maxsqez;

    s_c2=s_c-s_lm*s_dc;
    s_f2 = deff(s_c2);
    lc2 = s_f2'w*s_f2;
    clear s_f2;

    if lc1 <= lc2 and lc1 <= obj;

       retp(s_c1,s_itr-1); goto eoproc;

    else;

       s_c1=s_c2; s_lm=s_lm/2; s_c2=s_c - s_lm*s_dc; lc1=lc2;
       s_itr=s_itr+1;

    endif;

  endo;

retp(s_c2,s_itr-1);
eoproc:

endp;

@  ------------------------------------------------------ @
@ Subroutine to compute optimal weighting matrix. @

proc optw(a);

local f,ff;

f = zeros(n,neq);

f[.,1] = y2x - a[1,1] + (-2*Eyx_z'inEzz*zy_d' - Ey2_z'inEzz*zx_d')';
f[.,2] = yx2 - a[2,1] + (-Ex2_z'inEzz*zy_d' - 2*Eyx_z'inEzz*zx_d')';

ff = moment(f,0);

clear f;

    retp(invpd(ff./n));
    clear ff;

endp;
"Start time:  ";; timestr(ttt);
"End time:    ";; timestr(0);
"Start date:  ";; datestr(ddd);
"End date:    ";; datestr(0);

eop:
end;   @ This ends the program and closes all files.   @
@ ----------------------- END OF PROGRAM ------------------------------- @
