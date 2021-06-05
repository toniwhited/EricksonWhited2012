new ,60000;
ttt = time;
ddd = date;
call sysstate(14,1e-128);
output file=sumstata.out reset;
outwidth 256;
cls;
indata = "wdset";

numfirm = 0;
@--------Open GAUSS data set. -------------------------------------@

    closeall f1;
    open f1=^indata;


vnames = getname(indata);
let addon =  gvkey fyear ia mb cfa ik q cfk;
let addon =  gvkey fyear ik q cfk;
let addon =  gvkey fyear ia mb cfa;

x = readr(f1,rowsf(f1));

alldta = x[.,indcv(addon,vnames)];
alldta = selif(alldta,alldta[.,2].>1966 .and alldta[.,2].<2009);

format /rd 12,4;
"Summary stats";
$addon';
"Max    ";; maxc(alldta)';
"Min    ";;minc(alldta)';
"Std    ";;stdc(alldta)';
"Mean   ";;meanc(alldta)';
"Median ";;median(alldta)';
"N      ";; rows(alldta);
print;
let aidx = 3 4 5;
"Correlation Matrix";
cmx = corrx(alldta[.,aidx]);
$addon[aidx]';
cmx;

print;
"Covariance Matrix";
cmx = moment(alldta[.,aidx] - meanc(alldta[.,aidx])',0)/rows(alldta);
$addon';
cmx;

gvkey = alldta[.,1];
gvs = unique(gvkey,1);
newdta = zeros(1,cols(alldta)+3);
for ii(1,rows(gvs),1);
   dtai = selif(alldta,gvkey.==gvs[ii,1]);
   if rows(dtai)>5 and sumc(dtai[.,2].==seqa(dtai[1,2],1,rows(dtai)))==rows(dtai);
   lvs = dtai[.,aidx];      @lvs means lagged variables@
   lvs = (miss(1,1)*ones(1,3))|lvs[1:rows(lvs)-1,.];
   dtai = packr(dtai~lvs);
  @ dtai = dtai~seqa(1,1,rows(dtai))~(rows(dtai)*ones(rows(dtai),1));   @
   newdta = newdta|dtai;
   endif;
endfor;

newdta = newdta[2:rows(newdta),.];



ydum =design(newdta[.,2]-1973);
ars = zeros(6,1);
for yy(1,3,1);               @This is clearly doing the serial correlations@

yt0 = newdta[.,yy+2];
yt1 = newdta[.,yy+5];
des=yt1~ydum;
dd = invpd(moment(des,0));
allb=dd*(des'yt0);
ars[yy,1] = allb[1,1];
endfor;


y = alldta[.,aidx];                            @This is doing the partialled out variables @
for yy(1,2,1);
des=ones(rows(alldta),1)~alldta[.,aidx[3]];
dd = invpd(moment(des,0));

yt1 = alldta[.,aidx[yy]];
allb=dd*(des'yt1);
y = y~(yt1-des*allb);

if yy==1; "Muy  ";; allb[2,1]; else; "Mux  ";; allb[2,1];endif;
endfor;



des=ones(rows(alldta),1)~alldta[.,aidx[2:3]];
dd = invpd(moment(des,0));

yt1 = alldta[.,aidx[1]];
allb=dd*(des'yt1);
"Coefficient on q         ";; allb[2,1];
"Coefficient on cash flow ";; allb[3,1];

y = y~(yt1-des*allb);





"Covariance matrix of raw variables and residuals ";
moment(y-meanc(y)',0)/rows(y); print;
format /rd 12,3;
addon = addon[aidx,1] | "ikresid" | "qresid" | "uhat" ;
ncount = rows(y);
"                        Mean           Variance          Skew         Kurtosis          5ness            Min           Max           6ness           7ness            AR1        ";
for yy(1,6,1);
    Ey      = meanc(y[.,yy]);
    Vy      = (meanc(y[.,yy]^2) - Ey^2)*ncount/(ncount - 1);
    SDy     = sqrt(Vy);
    thirdy  = meanc((y[.,yy] - Ey)^3);
    skewy   = thirdy/(SDy^3);
    fourthy = meanc((y[.,yy] - Ey)^4);
    kurty   = (fourthy/(SDy^4));
    fifthy  = meanc((y[.,yy]- Ey)^5);
    std5y   = (fifthy)/(SDy^5);
    sixthy  = meanc((y[.,yy] - Ey)^6);
    std6y   = (sixthy)/(SDy^6);
    seveny  = meanc((y[.,yy] - Ey)^7);
    std7y   = (seveny)/(SDy^7);

    $addon[yy];; " & ";; Ey;; " & ";; Vy;; " & ";; skewy;; " & ";; kurty;; " & ";; std5y;; " & ";; minc(y[.,yy]);; " & ";; maxc(y[.,yy]);; " & ";; std6y;; " & ";; std7y;; " & ";; ars[yy,1];
endfor;

" ---------------------Ey2x-----------------------------";

zork = (y[.,4]^2).*y[.,5];
meanc(zork)
stdc(zork)/sqrt(rows(zork));

" ---------------------Eyx2-----------------------------";

zork = (y[.,5]^2).*y[.,4];
meanc(zork)
stdc(zork)/sqrt(rows(zork));
end;



print; "year by year"; print;

"                           Mean       Variance       Skew      Kurtosis       5ness       6ness       7ness       N ";
for year(1973,2005,1);

y = selif(alldta[.,3 4 5],alldta[.,2].==year);
ncount = rows(y);
for yy(1,3,1);
    Ey      = meanc(y[.,yy]);
    Vy      = (meanc(y[.,yy]^2) - Ey^2)*ncount/(ncount - 1);
    SDy     = sqrt(Vy);
    thirdy  = meanc((y[.,yy] - Ey)^3);
    skewy   = thirdy/(SDy^3);
    fourthy = meanc((y[.,yy] - Ey)^4);
    kurty   = (fourthy/(SDy^4));
    fifthy  = meanc((y[.,yy]- Ey)^5);
    std5y   = (fifthy)/(SDy^5);
    sixthy  = meanc((y[.,yy] - Ey)^6);
    std6y   = (sixthy)/(SDy^6);
    seveny  = meanc((y[.,yy] - Ey)^7);
    std7y   = (seveny)/(SDy^7);

    year;; $addon[yy+2];; Ey;; Vy;; skewy;; kurty;; std5y;; std6y;; std7y;; rows(y);
endfor;
endfor;
"-------------------------------------------------------------------------------------------------------";
"-------------------------------------------------------------------------------------------------------";
"-----------------------------------Fixed Effects-------------------------------------------------------";
"-------------------------------------------------------------------------------------------------------";
"-------------------------------------------------------------------------------------------------------";
   let indx = 3 4 5;
   newdta = zeros(1,cols(alldta));
   gvs = unique(alldta[.,1],1);
   for ii(1,rows(gvs),1);
       dtai = selif(alldta,alldta[.,1].==gvs[ii,1]);
       if rows(dtai)>1;
         dtai[.,indx] = dtai[.,indx]-meanc(dtai[.,indx])';
         newdta = newdta|dtai;
       endif;
   endfor;
   alldta = newdta[2:rows(newdta),.];



format /rd 12,4;
"Summary stats";
$addon';
"Max    ";; maxc(alldta)';
"Min    ";;minc(alldta)';
"Std    ";;stdc(alldta)';
"Mean   ";;meanc(alldta)';
"Median ";;median(alldta)';

print;
"Correlation Matrix";
cmx = corrx(alldta);
$addon';
cmx;


print;
"Covariance Matrix";
cmx = moment(alldta - meanc(alldta)',0)/rows(alldta);
$addon';
cmx;






y = alldta[.,3 4 5];
for yy(1,2,1);
des=ones(rows(alldta),1)~alldta[.,5];
dd = invpd(moment(des,0));

yt1 = alldta[.,yy+2];
allb=dd*(des'yt1);
y = y~(yt1-des*allb);

if yy==1; "Muy  ";; allb[2,1]; else; "Mux  ";; allb[2,1];endif;
endfor;

ncount = rows(y);
"                   Mean       Variance       Skew      Kurtosis       5ness       6ness       7ness          AR1 ";
for yy(1,5,1);
    Ey      = meanc(y[.,yy]);
    Vy      = (meanc(y[.,yy]^2) - Ey^2)*ncount/(ncount - 1);
    SDy     = sqrt(Vy);
    thirdy  = meanc((y[.,yy] - Ey)^3);
    skewy   = thirdy/(SDy^3);
    fourthy = meanc((y[.,yy] - Ey)^4);
    kurty   = (fourthy/(SDy^4));
    fifthy  = meanc((y[.,yy]- Ey)^5);
    std5y   = (fifthy)/(SDy^5);
    sixthy  = meanc((y[.,yy] - Ey)^6);
    std6y   = (sixthy)/(SDy^6);
    seveny  = meanc((y[.,yy] - Ey)^7);
    std7y   = (seveny)/(SDy^7);

    $addon[yy+2];; Ey;; Vy;; skewy;; kurty;; std5y;; std6y;; std7y;; ars[yy,1];
endfor;

print; "year by year"; print;

"                           Mean       Variance       Skew      Kurtosis       5ness       6ness       7ness       N ";
for year(1973,2005,1);

y = selif(alldta[.,3 4 5],alldta[.,2].==year);
ncount = rows(y);
for yy(1,3,1);
    Ey      = meanc(y[.,yy]);
    Vy      = (meanc(y[.,yy]^2) - Ey^2)*ncount/(ncount - 1);
    SDy     = sqrt(Vy);
    thirdy  = meanc((y[.,yy] - Ey)^3);
    skewy   = thirdy/(SDy^3);
    fourthy = meanc((y[.,yy] - Ey)^4);
    kurty   = (fourthy/(SDy^4));
    fifthy  = meanc((y[.,yy]- Ey)^5);
    std5y   = (fifthy)/(SDy^5);
    sixthy  = meanc((y[.,yy] - Ey)^6);
    std6y   = (sixthy)/(SDy^6);
    seveny  = meanc((y[.,yy] - Ey)^7);
    std7y   = (seveny)/(SDy^7);

    year;; $addon[yy+2];; Ey;; Vy;; skewy;; kurty;; std5y;; std6y;; std7y;; rows(y);
endfor;
endfor;

end;
