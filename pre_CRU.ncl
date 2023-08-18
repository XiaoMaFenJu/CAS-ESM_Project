begin

f1 = addfile("/mnt/hgfs/iap_down0711/cru_ts4.07.1901.2022.tmp.dat.nc","r")
f2 = addfile("/mnt/hgfs/iap_down0711/cru_ts4.07.1901.2022.pre.dat.nc","r")

tmp = f1->tmp(588:1367,:,:)
pre = f2->pre(588:1367,:,:)

printVarSummary(tmp)
out1 = addfile("/mnt/hgfs/iap_down0711/pre_CRU.nc","c")
out1->pre = pre
out2 = addfile("/mnt/hgfs/iap_down0711/t2m_CRU.nc","c")
out2->t2m = tmp

end