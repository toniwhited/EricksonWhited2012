new ,60000;
ttt = time;
ddd = date;
call sysstate(14,1e-128);
format /rd 12,5;
dataset = "wdset";
assets = 1;



output file = simpleiv.out reset;

if assets == 0;
let keepnames = gvkey fyear ik q cfk numobs nyr;
else;
let keepnames = gvkey fyear ia mb cfa numobs nyr;
endif;

outwidth  256;
closeall f1;
open f1=^dataset;

vnames = getname(dataset);

dta = readr(f1,rowsf(f1));
dta = selif(dta,dta[.,indcv("nyr",vnames)].>=4);

dta = dta[.,indcv(keepnames,vnames)];

year = dta[.,2];
ally = dta[.,3];
allx = dta[.,4];
allz = dta[.,5];
yr   = dta[.,6];
nyr  = dta[.,7];

ydum = year.==seqa(minc(year),1,maxc(year)-minc(year)+1)';
des = allx~allz~ones(rows(allx),1)~ydum[.,2:cols(ydum)];
dd = invpd(moment(des,0));
allb = dd*des'ally;
uhat = ally - des*allb;
wh = uhat.*des;
seb = sqrt(diag(dd*(wh'wh)*dd));

r2 = 1 - moment(uhat,0)/(moment((ally - meanc(ally)),0));
allb[1:2,1]~seb[1:2,1];
r2;








yt0 = selif(ally,yr.>3);
xt0 = selif(allx,yr.>3);
zt0 = selif(allz,yr.>3);

yt1 = selif(ally,yr.<(nyr-0) .and yr.>2);
xt1 = selif(allx,yr.<(nyr-0) .and yr.>2);
zt1 = selif(allz,yr.<(nyr-0) .and yr.>2);

yt2 = selif(ally,yr.<(nyr-1) .and yr.>1);
xt2 = selif(allx,yr.<(nyr-1) .and yr.>1);
zt2 = selif(allz,yr.<(nyr-1) .and yr.>1);

yt3 = selif(ally,yr.<(nyr-2));
xt3 = selif(allx,yr.<(nyr-2));
zt3 = selif(allz,yr.<(nyr-2));

year = selif(year,yr.>3);
ydum = year.==seqa(minc(year),1,maxc(year)-minc(year)+1)';







des = (xt0-xt1)~(zt0-zt1)~ones(rows(xt0),1)~ydum[.,2:cols(ydum)];
ivs = (zt0-zt1)~xt2~zt2~(xt3)~(zt3)~ones(rows(xt0),1)~ydum[.,2:cols(ydum)];
depvar = yt0-yt1;
"------------------------------------------------------------------------------";
" Differenced regressors, time dummies, 2 and 3 lags q and cf for instruments. ";
"------------------------------------------------------------------------------";
call do2sls(des,ivs,depvar);


proc do2sls(des,ivs,depvar);
local inezz,inexz,biv,ehat,wh,seiv,eff,w,eeff,jtest;
inezz = invpd(moment(ivs,0));
      inexz = invpd(des'ivs*inezz*ivs'des);
biv = inexz*des'ivs*inezz*ivs'depvar;

ehat = depvar - des*biv;
wh = ehat.*ivs;
seiv = sqrt(diag(inexz*des'ivs*inezz*(wh'wh)*inezz*ivs'des*inexz));
eff = ivs.*ehat;
w = invpd(moment(eff,0)/rows(eff));
Eeff = meanc(eff);
jtest = rows(eff)*Eeff'w*Eeff;
@
jtest = ehat'ivs*inezz*ivs'ehat/meanc(ehat^2);
@

biv[1:2,1]~seiv[1:2,1];
jtest;
cdfchic(jtest,cols(ivs)-cols(des));


retp(1);
endp;



"Start time:  ";; timestr(ttt);
"End time:    ";; timestr(0);
"Start date:  ";; datestr(ddd);
"End date:    ";; datestr(0);

output off;
end;
