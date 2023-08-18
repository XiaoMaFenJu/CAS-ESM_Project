begin

f11 = addfile("/mnt/hgfs/iap_down0711/prc.nc","r")
f12 = addfile("/mnt/hgfs/iap_down0711/prl.nc","r")
print(f11)
f2 = addfile("/mnt/hgfs/iap_down0711/cru_ts4.07.1901.2022.pre.dat.nc","r")

prc = f11->prc
prl = f12->prl
low_lat = f11->lat
low_lon1 = f11->lon
low_lon = low_lon1 - 179.75
printVarSummary(low_lon)

pre1 = prc+prl
pre = pre1
do i = 0,64,1
    do j = 0,11,1
        if j .eq. 1 then
            pre(i*12+j,:,:) = pre1(i*12+j,:,:) * 3600 * 24 * 28
        else if (j .eq. 0 .or. j .eq. 2 .or.j .eq. 4 .or.j .eq. 6 .or.j .eq. 7 .or.j .eq. 9 .or.j .eq. 11 ) then 
            pre(i*12+j,:,:) = pre1(i*12+j,:,:) * 3600 * 24 * 31
            else
        pre(i*12+j,:,:) = pre1(i*12+j,:,:) * 3600 * 24 * 30
            end if
        end if
    end do
end do

copy_VarCoords(prc,pre)

tmp = f2->pre
printVarSummary(tmp)

high_lat = f2->lat
high_lon = f2->lon

t2m_new = linint2_Wrap(low_lon, low_lat, pre, True, high_lon, high_lat, 0)
printVarSummary(t2m_new)


out = addfile("/mnt/hgfs/iap_down0711/pre_new.nc","c")
out->pre = t2m_new

end