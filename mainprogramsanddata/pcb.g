new ,60000;
ttt = time;
ddd = date;
call sysstate(14,1e-128);
output file=pcb.out reset;
outwidth 256;
cls;
indata = "wdset";

numfirm = 0;
@--------Open GAUSS data set. -------------------------------------@

    closeall f1;
    open f1=^indata;


vnames = getname(indata);
let addon =  gvkey fyear ia mb cfa ik q cfk;

x = readr(f1,rowsf(f1));
for ddd(1,2,1);
if ddd==1;
let addon =  gvkey fyear ik q cfk;
else;
let addon =  gvkey fyear ia mb cfa;
endif;

alldta = x[.,indcv(addon,vnames)];
alldta = selif(alldta,alldta[.,2].>1966 .and alldta[.,2].<2009);

format /rd 12,4;
let aidx = 3 4 5;

gvkey = alldta[.,1];
gvs = unique(gvkey,1);
newdta = zeros(1,cols(alldta)+6);
for ii(1,rows(gvs),1);
   dtai = selif(alldta,gvkey.==gvs[ii,1]);
   if rows(dtai)>2 and sumc(dtai[.,2].==seqa(dtai[1,2],1,rows(dtai)))==rows(dtai);
   lvs = dtai[.,aidx];      @lvs means lagged variables@
   lvs = (miss(1,1)*ones(1,3))|lvs[1:rows(lvs)-1,.];
   dvs = dtai[.,aidx]-lvs;
   ldvs = (miss(1,1)*ones(1,3))|dvs[1:rows(dvs)-1,.];
   dtai = packr(dtai~dvs~ldvs);
  @ dtai = dtai~seqa(1,1,rows(dtai))~(rows(dtai)*ones(rows(dtai),1));   @
   newdta = newdta|dtai;
   endif;
endfor;

newdta = newdta[2:rows(newdta),.];


ars = zeros(3,1);
for yy(1,3,1);               @This is clearly doing the serial correlations@

dv0 = newdta[.,yy+5];
dv1 = newdta[.,yy+8];
des = dv1;
y = 2*dv0 + dv1;
allb = invpd(moment(des,0))*des'y;
ars[yy,1] = allb[1,1];
endfor;
if ddd==1; "IK    q     CF";
else;      "IA    MB    CF";
endif;
ars';

endfor;

end;
