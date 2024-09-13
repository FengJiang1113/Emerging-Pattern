load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
resultdir="/data0/fjiang/scripts/contrast/"


;;;;;ORAs5 wspeed pattern
opendir1="/data0/fjiang/data/ocean/oras5/postprocess/"
a = addfile(opendir1+"zs.remap.nc", "r")
b = addfile(opendir1+"ms.remap.nc", "r")

syear1=1980;58;;58;80
syear2=1958;58

latr=30

eyear    = 2022
ymdStrt1 = syear1*100+1
ymdStrt2 = syear2*100+1

ymdLast1 = eyear*100+12
mm1=(eyear-syear1+1)*12
mm2=(eyear-syear2+1)*12


time1      = a->time_counter
ymd1      = cd_calendar(time1, -1)
iStrt1    = ind(ymd1.eq.ymdStrt1)
iStrt2    = ind(ymd1.eq.ymdStrt2)
iLast1    = ind(ymd1.eq.ymdLast1)

delete(time1)
delete(ymd1)

oras5_taux      =a->sozotaux
oras5_tauy      =b->sometauy
lat=oras5_taux&lat
lon=oras5_taux&lon

Curl = uv2dv_cfd(oras5_taux,oras5_tauy,lat,lon,2)
copy_VarMeta(oras5_taux,Curl)

oras5_wspeed=sqrt(oras5_taux*oras5_taux+oras5_tauy*oras5_tauy)
copy_VarMeta(oras5_taux,oras5_wspeed)
Curl_annual1 =month_to_annual(Curl(iStrt1:iLast1,:,:),1)
Curl_annual2 =month_to_annual(Curl(iStrt2:iLast1,:,:),1)

oras5_wspeed_annual1 =month_to_annual(oras5_wspeed(iStrt1:iLast1,:,:),1)
oras5_wspeed_annual2 =month_to_annual(oras5_wspeed(iStrt2:iLast1,:,:),1)

oras5_taux_annual1 =month_to_annual(oras5_taux(iStrt1:iLast1,:,:),1)
oras5_taux_annual2 =month_to_annual(oras5_taux(iStrt2:iLast1,:,:),1)


oras5_tauy_annual1 =month_to_annual(oras5_tauy(iStrt1:iLast1,:,:),1)
oras5_tauy_annual2 =month_to_annual(oras5_tauy(iStrt2:iLast1,:,:),1)

year1    = ispan(syear1,eyear,1)
year2    = ispan(syear2,eyear,1)

dims=dimsizes(oras5_wspeed_annual1)
rc1    = regCoef_n(year1,oras5_wspeed_annual1,0,0)
rc2    = regCoef_n(year2,oras5_wspeed_annual2,0,0)

rc1_taux    = regCoef_n(year1,oras5_taux_annual1,0,0)
rc2_taux    = regCoef_n(year2,oras5_taux_annual2,0,0)

rc1_tauy    = regCoef_n(year1,oras5_tauy_annual1,0,0)
rc2_tauy    = regCoef_n(year2,oras5_tauy_annual2,0,0)

rc1_Curl    = regCoef_n(year1,Curl_annual1,0,0)
rc2_Curl    = regCoef_n(year2,Curl_annual2,0,0)

rc1_taux:= rc1_taux*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
rc2_taux:= rc2_taux*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
rc1_tauy:= rc1_tauy*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
rc2_tauy:= rc2_tauy*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
rc1_Curl:= rc1_Curl*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
rc2_Curl:= rc2_Curl*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
copy_VarMeta(oras5_wspeed(0,:,:),rc1_taux)
copy_VarMeta(oras5_wspeed(0,:,:),rc2_taux)
copy_VarMeta(oras5_wspeed(0,:,:),rc1_tauy)
copy_VarMeta(oras5_wspeed(0,:,:),rc2_tauy)
copy_VarMeta(oras5_wspeed(0,:,:),rc1_Curl)
copy_VarMeta(oras5_wspeed(0,:,:),rc2_Curl)
;ddddddddddddrc1_taux:=dim_avg_n_Wrap(oras5_taux,0)



confi1=student_t(rc1@tval,rc1@nptxy-2)
confi1:=onedtond(confi1,(/dims(1),dims(2)/))
rc1@long_name = "Trend"
rc1:= rc1*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
copy_VarMeta(oras5_wspeed(0,:,:),rc1)
copy_VarMeta(rc1,confi1)




;confi2=student_t(rc2@tval,rc2@nptxy-2)
;confi2:=onedtond(confi2,(/dims(1),dims(2)/))
;rc2@long_name = "Trend"
;rc2:= rc2*10;(eyear-syear2+1);/10                      ; (C/year)*(nyrs)
;copy_VarMeta(oras5_wspeed(0,:,:),rc2)
;copy_VarMeta(rc2,confi2)
lonL=120
lonR=310;270

d = addfile("/data0/fjiang/data/ocean/oras5/oras5.d20.all.nc", "r")
d20       =d->D20(:,{-latr:latr},{120:lonR})
d20_he1=month_to_annual(d20(iStrt1:iLast1,:,:),1)
d20_he2=month_to_annual(d20(iStrt2:iLast1,:,:),1)


d20_he1:=d20_he1(:,{-latr:latr},{lonR})
d20_he2:=d20_he2(:,{-latr:latr},{lonR})

delete(d20_he1@_FillValue)
delete(d20_he2@_FillValue)

d20_he1=where(d20_he1.gt.10e3.or.d20_he1.lt.-10e3,0.,d20_he1)
d20_he2=where(d20_he2.gt.10e3.or.d20_he2.lt.-10e3,0.,d20_he2)


rc1_he= 10*regCoef_n(year1,d20_he1,0,0)
rc2_he= 10*regCoef_n(year2,d20_he2,0,0)

d20:=dim_avg_n_Wrap(d20,0)
HH=d20



;d20:=dim_avg_n_Wrap(d20,1)
dims_d20=dimsizes(d20)
he=0;d20(:,dims_d20(1)-1)
he1=d20(:,dims_d20(1)-2)
;rc1_taux:=1
;rc1_taux({-30:30},{120:270})=-1

HH_ori=HH
delete(HH_ori@_FillValue)
HH_ori=where(HH_ori.gt.1000,99999,HH_ori)
;Tau1_value=dim_avg_n_Wrap(Tau_ref1,1)

;;;;;;;;Simplest Case :Uniform Wind Stress in the zonal direction and no meridonal component


Tau_ref1=rc1_taux({-latr:latr},{lonL:lonR})
Tau_ref2=rc2_taux({-latr:latr},{lonL:lonR})
Tau_ref =rc1_taux({-latr:latr},{lonL:lonR})

Tauy_ref1=rc1_tauy({-latr:latr},{lonL:lonR})
Tauy_ref2=rc2_tauy({-latr:latr},{lonL:lonR})

delete(dims)
lat:=doubletofloat(Tau_ref1&lat)
lon:=doubletofloat(Tau_ref1&lon)
lamdae=20
lamda=lonR
delete(a)
a= 6378388.                  ; m
rho1=1025.0                  ; density of seawater (kg/m^3)
rho2=1027.7                  ; density of seawater (kg/m^3)
g=9.8                        ; m/s2
delta_rho=rho2-rho1       ;
g1=g*delta_rho/rho2          ;0.019 m/s2
print(g1)
pi=3.1415926
phi=lat;*0.0174533););*lat(nl) (pi/180.0));;;*6378388.
print(lat)
print(phi)
dims=dimsizes(Tau_ref1)
nlat=dims(0)
nlon=dims(1)
;taux_phi=rc1_taux
;deltalamda=(lon(1)-lon(0))* (pi/180.0)*6378388.
omega=2*pi/86400
f=2*omega*sin(phi*pi/180)
f=where(f.eq.0,f@_FillValue,f)

Tau1=Tau_ref1
Tau2=Tau_ref2

Tau1_y=Tauy_ref1
Tau2_y=Tauy_ref2

;Tau1_value=dim_avg_n_Wrap(Tau_ref1,1)
dims_Tau=dimsizes(Tau_ref1)


;end do
;print(Tau1_value)
delete(Tau1@_FillValue)
delete(Tau2@_FillValue)
delete(Tau1_y@_FillValue)
delete(Tau2_y@_FillValue)
delete(HH@_FillValue)

Tau1_rm=where(Tau1.gt.10e1.or.Tau1.lt.-10e1,0,Tau1)
Tau2_rm=where(Tau2.gt.10e1.or.Tau2.lt.-10e1,0,Tau2)
Tau1_y_rm=where(Tau1_y.gt.10e1.or.Tau1_y.lt.-10e1,0,Tau1_y)
Tau2_y_rm=where(Tau2_y.gt.10e1.or.Tau2_y.lt.-10e1,0,Tau2_y)
HH_rm=where(HH.gt.10e3.or.HH.lt.-10e3,99999,HH)

;HH_rm=where(HH.gt.10e8.or.HH.lt.-10e8,999,HH)
copy_VarMeta(Tau1,Tau1_rm)
copy_VarMeta(Tau2,Tau2_rm)
copy_VarMeta(Tau1_y,Tau1_y_rm)
copy_VarMeta(Tau2_y,Tau2_y_rm)
copy_VarMeta(HH,HH_rm)
HH=150.;HH_rm
Tau1=Tau1_rm
Tau2=Tau2_rm
Tau1_y=Tau1_y_rm
Tau2_y=Tau2_y_rm
;copy_VarMeta(rc1_taux,rc1_taux_rm)
print(HH)


;;;;;Get a voundary!
he_lon=HH(:,0)
do j=0,dims_Tau(0)-1
 he_lon(j)=min(ind(HH_ori(j,{240:lonR}).eq.99999))+240
end do
;;;;;;get rid of some equtorial ilands
do j=18,22;dims_Tau(0)-1
 he_lon(j)=min(ind(HH_ori(j,{275:lonR}).eq.99999))+275
end do

;;;;;;;;Process the boundary!
do j=0,dims_Tau(0)-1
  aa=floattoint(he_lon(j)-120)
  print(aa)
  HH(j,aa::)=99999
end do
print(HH)

delete(HH@_FillValue)
HH=where(HH.gt.1000.or.HH.lt.10,99999,HH)



taux_phi1 =Tau1
tauy_phi1 =Tau1_y

h_square1 =Tau1
tau_phi1=0.
he_square1=new(nlat,float)
he_square1=rc1_he;*he;10000;2500.

taux_phi2 =Tau2
tauy_phi2 =Tau2_y

h_square2 =Tau2
tau_phi2  =0.
he_square2=new(nlat,float)
he_square2=rc2_he;0.;he*he;10000;2500.
print(rc2_he)

Tau1=Tau1/1025
Tau2=Tau2/1025

Tau1_y=Tau1_y/1025
Tau2_y=Tau2_y/1025

dlat  = (lat(2)-lat(1))*0.0174533 ; convert to radians0.25*0.0174533 ; convert to radians
dY    = doubletofloat(dlat);;;6378388.*

taux_phi1=center_finite_diff_n(Tau1, pi/180.0 ,False,0,0);(numerator(j,i)-numerator(j-1,i))/(phi(j)-phi(j-1))
taux_phi2=center_finite_diff_n(Tau2, pi/180.0 ,False,0,0);(numerator(j,i)-numerator(j-1,i))/(phi(j)-phi(j-1))

tauy_phi1=center_finite_diff_n(Tau1_y, pi/180.0 ,False,0,1);(numerator(j,i)-numerator(j-1,i))/(phi(j)-phi(j-1))
tauy_phi2=center_finite_diff_n(Tau2_y, pi/180.0 ,False,0,1);(numerator(j,i)-numerator(j-1,i))/(phi(j)-phi(j-1))


lamdae=lonR

;lon:=where(lon.gt.180,long-180,lon)


h_square1=0.
h_square2=0.

;;;;;simple hypothesis
;;;;;  do j=0,dims_Tau(0)-1
;;;;;    do i=0,dims_Tau(1)-1
;;;;;   ;(j,i)=-2*3.2*10e7*(Tau1(j,i)*cos(lat(j)))*(lamdae-lon(i))* (pi/180.0);*6378388.
;;;;;;h_square(j,i)=a/(g1)*Tau1(j,i)*(lamdae-lon(i))* (pi/180.0)/HH(j,i);*6378388.
;;;;;   h_square1(j,i)=he_square1(j)-a/(g1)*(sin(lat(j)*pi/180.0)*sin(lat(j)*pi/180.0)/cos(lat(j)*pi/180.0)*Tau1(j,i)+Tau1(j,i)*cos(lat(j)*pi/180.0)-sin(lat(j)*pi/180.0)*taux_phi1(j,i))*(lamdae-lon(i))* (pi/180.0)/HH(j,i);*6378388.
;;;;;   h_square2(j,i)=he_square2(j)-a/(g1)*(sin(lat(j)*pi/180.0)*sin(lat(j)*pi/180.0)/cos(lat(j)*pi/180.0)*Tau2(j,i)+Tau2(j,i)*cos(lat(j)*pi/180.0)-sin(lat(j)*pi/180.0)*taux_phi2(j,i))*(lamdae-lon(i))* (pi/180.0)/HH(j,i);*6378388.
;;;;;   ;hh(j,i)=h_square(j,i)+he_square(j)
;;;;;  end do
;;;;;end do

f1=h_square1
f2=h_square1
rr=1./(5.5*365*86400);;;-s
A=HH
A=0.
;HH=100;
HH_phix=center_finite_diff_n(HH, pi/180.0 ,False,0,1);(numerator(j,i)-numerator(j-1,i))/(phi(j)-phi(j-1))

HH_phix=where(HH_phix.gt.10e2.or.HH_phix.lt.-10e2,0.,HH_phix)
print(HH_phix)


;;;zonally-varying hypothesis
beta=2.3*1e-11;l;(/(m*s);2*omega*cos(lat(nl) * (pi/180.0))

do j=0,dims_Tau(0)-1
  do i=0,dims_Tau(1)-1
    A(j,i)=1/HH(j,i)*(-a*a*a/g1*beta*(lat(j)*pi/180.0)*(lat(j)*pi/180.0)*cos(lat(j)*pi/180.0)*rr);;;HH_phix(j,i)
  end do
end do


;A=0.001;where(A.gt.10e1.or.A.lt.-10e1,0.,A)
print(A)
;print(aaa)

do j=0,dims_Tau(0)-1
  do i=0,dims_Tau(1)-1
 ;(j,i)=-2*3.2*10e7*(Tau1(j,i)*cos(lat(j)))*(lamdae-lon(i))* (pi/180.0);*6378388.
;h_square(j,i)=a/(g1)*Tau1(j,i)*(lamdae-lon(i))* (pi/180.0)/HH(j,i);*6378388.

 f1(j,i)=-a*exp(A(j,i)*lon(i)*pi/180.0)/(2*g1)*(cos(lat(j)*pi/180.0)*Tau1(j,i)+Tau1(j,i)*sin(lat(j)*pi/180.0)*(lat(j)*pi/180.0)-cos(lat(j)*pi/180.0)*(lat(j)*pi/180.0)*taux_phi1(j,i)+tauy_phi1(j,i)*(lat(j)*pi/180.0))* (pi/180.0)/HH(j,i);*6378388.
 f2(j,i)=-a*exp(A(j,i)*lon(i)*pi/180.0)/(2*g1)*(cos(lat(j)*pi/180.0)*Tau2(j,i)+Tau2(j,i)*sin(lat(j)*pi/180.0)*(lat(j)*pi/180.0)-cos(lat(j)*pi/180.0)*(lat(j)*pi/180.0)*taux_phi2(j,i)+tauy_phi2(j,i)*(lat(j)*pi/180.0))* (pi/180.0)/HH(j,i);*6378388.
;;
;f1(j,i)=-a/(2*g1)*(sin(lat(j)*pi/180.0)*sin(lat(j)*pi/180.0)/cos(lat(j)*pi/180.0)*Tau1(j,i)+Tau1(j,i)*cos(lat(j)*pi/180.0)-sin(lat(j)*pi/180.0)*taux_phi1(j,i))* (pi/180.0)/HH(j,i);*6378388.
;f2(j,i)=-a/(2*g1)*(sin(lat(j)*pi/180.0)*sin(lat(j)*pi/180.0)/cos(lat(j)*pi/180.0)*Tau1(j,i)+Tau2(j,i)*cos(lat(j)*pi/180.0)-sin(lat(j)*pi/180.0)*taux_phi2(j,i))* (pi/180.0)/HH(j,i);*6378388.

end do
end do
f1=where(f1.gt.1e1.or.f1.lt.-1e1,0.,f1)
f2=where(f2.gt.1e1.or.f2.lt.-1e1,0.,f2)



int1= h_square1
int2= h_square1

do j=0,dims_Tau(0)-1
do i=1,dims_Tau(1)-1
  kk=dims_Tau(1)-1-i
  int1(j,kk)=(int1(j,kk+1)+f1(j,kk+1)+f1(j,kk))
  int2(j,kk)=(int2(j,kk+1)+f2(j,kk+1)+f2(j,kk))
 h_square1(j,kk)=exp(-A(j,kk)*lon(kk)*pi/180.0)*int1(j,kk)
 h_square2(j,kk)=exp(-A(j,kk)*lon(kk)*pi/180.0)*int2(j,kk)
 ;he_square2(j)-a/(g1)*(sin(lat(j)*pi/180.0)*sin(lat(j)*pi/180.0)/cos(lat(j)*pi/180.0)*Tau2(j,i)+Tau2(j,i)*cos(lat(j)*pi/180.0)-sin(lat(j)*pi/180.0)*taux_phi2(j,i))*(lamdae-lon(i))* (pi/180.0)/HH(j,i);*6378388.
 ;hh(j,i)=h_square(j,i)+he_square(j)
end do
end do

ff1=f1
ff2=f2

do j=0,dims_Tau(0)-1
  do i=0,dims_Tau(1)-1
 ;(j,i)=-2*3.2*10e7*(Tau1(j,i)*cos(lat(j)))*(lamdae-lon(i))* (pi/180.0);*6378388.
;h_square(j,i)=a/(g1)*Tau1(j,i)*(lamdae-lon(i))* (pi/180.0)/HH(j,i);*6378388.

 ff1(j,i)=exp(A(j,i)*lon(i)*pi/180.0)* (pi/180.0)/2.;/HH(j,i);*6378388.
 ff2(j,i)=exp(A(j,i)*lon(i)*pi/180.0)* (pi/180.0)/2.
;;
;f1(j,i)=-a/(2*g1)*(sin(lat(j)*pi/180.0)*sin(lat(j)*pi/180.0)/cos(lat(j)*pi/180.0)*Tau1(j,i)+Tau1(j,i)*cos(lat(j)*pi/180.0)-sin(lat(j)*pi/180.0)*taux_phi1(j,i))* (pi/180.0)/HH(j,i);*6378388.
;f2(j,i)=-a/(2*g1)*(sin(lat(j)*pi/180.0)*sin(lat(j)*pi/180.0)/cos(lat(j)*pi/180.0)*Tau1(j,i)+Tau2(j,i)*cos(lat(j)*pi/180.0)-sin(lat(j)*pi/180.0)*taux_phi2(j,i))* (pi/180.0)/HH(j,i);*6378388.

end do
end do


intt1= h_square1
intt2= h_square1

scale1=2*f1;h_square1
scale2=2*f2;h_square1
do j=0,dims_Tau(0)-1
do i=1,dims_Tau(1)-1
  kk=dims_Tau(1)-1-i
 ;scale1(j,kk)=exp(-A(j,kk)*lon(kk)*pi/180.0)*2*f1(j,kk)/(pi/180.0)
 ;scale2(j,kk)=exp(-A(j,kk)*lon(kk)*pi/180.0)*2*f2(j,kk)/(pi/180.0)
 scale1(j,kk)=2*f1(j,kk)/(pi/180.0)
 scale2(j,kk)=2*f2(j,kk)/(pi/180.0)



 ;intt1(j,kk)=(intt1(j,kk+1)+ff1(j,kk+1)+ff1(j,kk))
 ;intt2(j,kk)=(intt2(j,kk+1)+ff2(j,kk+1)+ff2(j,kk))
 ;scale1(j,kk)=exp(-A(j,kk)*lon(kk)*pi/180.0)*intt1(j,kk)
 ;scale2(j,kk)=exp(-A(j,kk)*lon(kk)*pi/180.0)*intt2(j,kk)

 ;f1(j,i)=-a*exp(A(j,i)*lon(i)*pi/180.0)/(2*g1)*(cos(lat(j)*pi/180.0)*Tau1(j,i)+Tau1(j,i)*sin(lat(j)*pi/180.0)*(lat(j)*pi/180.0)-cos(lat(j)*pi/180.0)*(lat(j)*pi/180.0)*taux_phi1(j,i)+tauy_phi1(j,i)*(lat(j)*pi/180.0))* (pi/180.0)/HH(j,i);*6378388.
 ;f2(j,i)=-a*exp(A(j,i)*lon(i)*pi/180.0)/(2*g1)*(cos(lat(j)*pi/180.0)*Tau2(j,i)+Tau2(j,i)*sin(lat(j)*pi/180.0)*(lat(j)*pi/180.0)-cos(lat(j)*pi/180.0)*(lat(j)*pi/180.0)*taux_phi2(j,i)+tauy_phi2(j,i)*(lat(j)*pi/180.0))* (pi/180.0)/HH(j,i);*6378388.
 ;;;
 ;he_square2(j)-a/(g1)*(sin(lat(j)*pi/180.0)*sin(lat(j)*pi/180.0)/cos(lat(j)*pi/180.0)*Tau2(j,i)+Tau2(j,i)*cos(lat(j)*pi/180.0)-sin(lat(j)*pi/180.0)*taux_phi2(j,i))*(lamdae-lon(i))* (pi/180.0)/HH(j,i);*6378388.
 ;hh(j,i)=h_square(j,i)+he_square(j)
end do
end do

print(scale1)
print(scale2)

;end do
;end do




copy_VarMeta(Tau_ref(:,:),scale1)
copy_VarMeta(Tau_ref(:,:),scale2)


h_square1:=h_square1/1027.7*2.7
h_square2:=h_square2/1027.7*2.7

copy_VarMeta(Tau_ref(:,:),h_square1)
copy_VarMeta(Tau_ref(:,:),h_square2)

opendir1="/data0/fjiang/data/ocean/oras5/postprocess/"
ss = addfile(opendir1+"ssh.remap.nc", "r")
syear1=1980
syear2=1958
eyear    = 2022
ymdStrt1 = syear1*100+1
ymdStrt2 = syear2*100+1

ymdLast1 = eyear*100+12
mm1=(eyear-syear1+1)*12
mm2=(eyear-syear2+1)*12


time1      = ss->time_counter
ymd1      = cd_calendar(time1, -1)
iStrt1    = ind(ymd1.eq.ymdStrt1)
iStrt2    = ind(ymd1.eq.ymdStrt2)
iLast1    = ind(ymd1.eq.ymdLast1)

delete(time1)
delete(ymd1)

ssh       =ss->sossheig(:,:,:)
ssh_annual1 =month_to_annual(ssh(iStrt1:iLast1,:,:),1)
ssh_annual2 =month_to_annual(ssh(iStrt2:iLast1,:,:),1)


year1    = ispan(syear1,eyear,1)
year2    = ispan(syear2,eyear,1)

rc1_ssh    = regCoef_n(year1,ssh_annual1,0,0)
rc2_ssh    = regCoef_n(year2,ssh_annual2,0,0)

dims_ssh=dimsizes(ssh)
confi1_ssh=student_t(rc1_ssh@tval,rc1_ssh@nptxy-2)
confi1_ssh:=onedtond(confi1_ssh,(/dims_ssh(1),dims_ssh(2)/))
copy_VarMeta(ssh(0,:,:),confi1_ssh)
confi2_ssh=student_t(rc2_ssh@tval,rc2_ssh@nptxy-2)
confi2_ssh:=onedtond(confi2_ssh,(/dims_ssh(1),dims_ssh(2)/))
copy_VarMeta(ssh(0,:,:),confi2_ssh)


rc1_ssh:= rc1_ssh*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
copy_VarMeta(ssh(0,:,:),rc1_ssh)
rc2_ssh:= rc2_ssh*10;(eyear-syear2+1);/10                      ; (C/year)*(nyrs)
copy_VarMeta(ssh(0,:,:),rc2_ssh)


printVarSummary(h_square1({-latr:latr},{120:270}))
printVarSummary(rc1_ssh({-latr:latr},{120:270}))

pc1=pattern_cor(h_square1({-15:15},{120:270}),rc1_ssh({-15:15},{120:270}),1.0,0)
pc2=pattern_cor(h_square2({-15:15},{120:270}),rc2_ssh({-15:15},{120:270}),1.0,0)
print(pc1)
print(pc2)

guess     = 1                ; use zonal means
is_cyclic = True             ; cyclic [global]
nscan     = 1500             ; usually much less than this
eps       = 1.e-2            ; variable dependent
relc      = 0.6              ; relaxation coefficient
opt       = 0                ; not used

poisson_grid_fill(rc1_ssh, is_cyclic, guess, nscan, eps, relc, opt)
poisson_grid_fill(rc2_ssh, is_cyclic, guess, nscan, eps, relc, opt)

poisson_grid_fill(scale1, is_cyclic, guess, nscan, eps, relc, opt)
poisson_grid_fill(scale2, is_cyclic, guess, nscan, eps, relc, opt)


wks = gsn_open_wks("eps",resultdir+"v6-4");1979-2018")
gsn_define_colormap(wks,"MPL_RdBu")
;gsn_define_colormap(wks,"MPL_RdBu");"CBR_drywet");
gsn_reverse_colormap(wks);
gsn_define_colormap(wks,"MPL_RdBu")
;gsn_define_colormap(wks,"MPL_RdBu");"CBR_drywet");
gsn_reverse_colormap(wks);
res = True
res@gsnDraw   = False
res@gsnFrame  = False
res@gsnAddCyclic      = False
res@mpShapeMode="FreeAspect"
;res@vpWidthF =0.9
res@mpFillOn           = False
res@mpMinLonF  = 120 ;lon2d(0,0)
res@mpCenterLonF=180
res@mpMaxLonF = 290 ;lon2d(0,ilon-1)
;res@lbBoxEndCapStyle="";TriangleBothEnds"

res@vpHeightF = 0.15
res@vpWidthF  = 0.34

 res@lbBoxLinesOn=True
res@mpMinLatF            = -20;latr      ; zoom in on map
res@mpMaxLatF            =  20;latr
;
res@pmTickMarkDisplayMode="Always"
res@mpFillOn=True
res@mpLandFillColor="gray60"
res@mpOceanFillColor=-1
res@mpFillDrawOrder="PostDraw"
res@tmXBMinorOn=False
;res@tmYLLabels =(/"30S","20S","10S","0","10N","20N","30N"/)
res@tmYLMinorOn=False
res@cnInfoLabelOn=False
res@mpGeophysicalLineThicknessF=1
res@mpGeophysicalLineColor="gray70"
res@pmTickMarkDisplayMode = "Always"
res@tmXBLabelFontHeightF  = 0.008
res@tmYLLabelFontHeightF  = 0.008

res@cnFillOn    = True
res@cnLinesOn = False

res@tmXBMajorLengthF=0.004
res@tmYLMajorLengthF=0.004
res@tmXBMajorOutwardLengthF=0.004
res@tmYLMajorOutwardLengthF=0.004
;
res@cnFillOn             = True
res@cnMaxLevelValF  =1.5
res@cnMaxLevelValF  =-1.5

res@cnLinesOn =False
res@cnLevelSelectionMode = "ExplicitLevels"  ; use explicit levels
;res@cnFillColors=(/111,100,96,87,77,67,64,54,44,34,30,20/);
res@cnFillColors=(/10,15,20,25,30,35,40,45,50,60,0,0,70,80,85,90,95,100,105,110,115,120/)

res@gsnLeftStringFontHeightF=0.012
res@gsnLeftStringOrthogonalPosF=0.02;0.005
res@gsnRightStringFontHeightF=0.012
res@gsnRightStringOrthogonalPosF=-0.003;5;0.005
res@cnLineLabelsOn=False
res@gsnRightString="";(~S~o~N~C/decade)"
res@gsnLeftString="";"(a) Wind-driven SSH Trend 1980-2022  "
res@tiMainString="(a) Observed SSH Trend 1980-2022"
res@tiMainFont=22
res@tiMainFontHeightF=0.0112
res@tiMainOffsetYF=-0.014
plot1=new(4,graphic)

;dum1=new(4,graphic)
res@vpXF=0.06
res@vpYF=0.9

res@vpHeightF= 0.15            ; change aspect ratio of plot

res@cnLevels             = 0.004*ispan(-10,10,1 )

res@lbLabelBarOn=False
res@gsnRightStringOrthogonalPosF=0.02;0.005
res@gsnRightString="";sprintf("%5.2f", pc1);(~S~o~N~C/decade)"
;plot1(0) = gsn_csm_contour_map(wks,rc1_ssh,res)
res@tiMainString="(b) Observed SSH Trend 1958-2022"

;plot1(1) = gsn_csm_contour_map(wks,rc2_ssh,res)

res@lbLabelBarOn=False
res@pmLabelBarOrthogonalPosF=0.2
res@pmLabelBarParallelPosF=-0.11
res@tiMainString="(a) Wind-driven SSH Trend 1980-2022 (R="+sprintf("%5.2f", pc1)+")"

res@cnLevels             = 0.004*ispan(-10,10,1 )
res@gsnRightStringOrthogonalPosF=0.02;0.005
res@gsnRightString="";sprintf("%5.2f", pc2);(~S~o~N~C/decade)"

plot1(0) = gsn_csm_contour_map(wks,h_square1,res)
res@tiMainString="(b) Wind-driven SSH Trend 1958-2022 (R="+sprintf("%5.2f", pc2)+")";(~S~o~N~C/decade)"
res@lbLabelBarOn=True
res@pmLabelBarHeightF=0.037
res@pmLabelBarWidthF=0.6
res@lbLabelStride=2.
res@lbLabelFontHeightF=0.0085
res@pmLabelBarParallelPosF=-0.12;OrthogonalPosF

plot1(1) = gsn_csm_contour_map(wks,h_square2,res)

delete(res@cnLevels)
delete(res@cnFillColors)

res@cnLevels             = 1.*ispan(-1,10,1 )
res@cnFillPalette="cmocean_dense"
res@cnFillColors=(/-1,2,15,30,40,50,70,80,100,120,130,140,150/)

print(scale1({-5:5},{150:200}))
res@lbLabelBarOn=False
res@tiMainString="(c) Wind Effect 1980-2022 ";(~S~o~N~C/decade)"

plot1(2) = gsn_csm_contour_map(wks,scale1,res)
res@tiMainString="(d) Wind Effect 1958-2022 ";(~S~o~N~C/decade)"

res@lbLabelBarOn=True
plot1(3) = gsn_csm_contour_map(wks,scale2,res)



res1=True
res1@gsnDraw=False
res1@gsnFrame=False
res1@cnFillOn=True
res1@cnLinesOn=False
res1@cnMonoFillPattern=False
res1@gsnLeftString=""
res1@gsnRightString=""
res1@cnLevelSelectionMode="ExplicitLevels"
res1@cnLevels=(/0,0.05/)
res1@cnFillPatterns=(/-1,17,-1/)
res1@cnFillColors=(/"black","black","black"/)
res1@cnLineLabelsOn=False
res1@cnInfoLabelOn=False
res1@lbLabelBarOn=False
  plot2=new(2,graphic)

  plot2(0)=gsn_csm_contour(wks,confi1_ssh,res1)
  plot2(1)=gsn_csm_contour(wks,confi2_ssh,res1)
;  plot2(2)=gsn_csm_contour(wks,confi3,res1)
;;  plot2(3)=gsn_csm_contour(wks,confi4,res1)
;;
;;;  plot2(1)=gsn_csm_contour(wks,confislp,res1)
;;
  ;overlay(plot1(0),plot2(0))
  ;overlay(plot1(1),plot2(1))
;;  overlay(plot1(2),plot2(2))
;;  overlay(plot1(3),plot2(3))

;  ;draw(plot1(2))
;  ;draw(plot1(3))
;  ;draw(plot1(4))
;  ;draw(plot1(5))
;
;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
delete(res)
res                 = True
res@gsnDraw         = False
res@gsnFrame        = False
res@vpHeightF = 0.15
res@vpWidthF  = 0.34
res@trXMinF   =  1870.
res@trXMaxF   =  2022.
res@trYMaxF   = 1.
res@trYMinF   = -1.
res@tiYAxisString = " "
res@tiXAxisString = " "
res@tiMainString  = " "
res@tmXBMajorThicknessF     = 1.5
res@tmYLMajorThicknessF     = 1.5
res@tmXBMinorThicknessF     = 1.0
res@tmYLMinorThicknessF     = 1.0
res@tmYLMinorOn=False
res@tmYLValues    = 0.5*ispan(-5,5,1)


txres=True
txres@txFontHeightF=0.012
;gsn_text_ndc(wks, "[m/decade]",0.503,0.483, txres)


pgres=True
pgres@gsMarkerIndex      = 12

pgres@gsMarkerSizeF = 4.
pgres@tfPolyDrawOrder ="PostDraw"
pgres@gsMarkerColor="orangered3"
;dum_gon_1 = gsn_add_polymarker (wks, plot1(2),151.2+1870,patterncor1_oras5, pgres)
;pgres@gsMarkerColor="gray40"
;pgres@gsMarkerIndex      = 16
;dum_gon_2 = gsn_add_polymarker (wks, plot1(2),151.2+1870,patterncor1_a_oras5, pgres)
;pgres@gsMarkerColor="deepskyblue3"
;pgres@gsMarkerIndex      = 14
;dum_gon_3 = gsn_add_polymarker (wks, plot1(2),151.2+1870,patterncor1_b_oras5, pgres)
;pgres@gsMarkerColor="orangered3"
;pgres@gsMarkerIndex      = 12
;dum_gon_4 = gsn_add_polymarker (wks, plot1(3),151.2+1870,patterncor2_oras5, pgres)
;pgres@gsMarkerColor="gray40"
;pgres@gsMarkerIndex      = 16
;dum_gon_5 = gsn_add_polymarker (wks, plot1(3),151.2+1870,patterncor2_a_oras5, pgres)
;pgres@gsMarkerColor="deepskyblue3"
;pgres@gsMarkerIndex      = 14
;dum_gon_6 = gsn_add_polymarker (wks, plot1(3),151.2+1870,patterncor2_b_oras5, pgres)

pgres@gsMarkerColor="orangered3"
pgres@gsMarkerIndex      = 16
pgres@gsMarkerSizeF = 8.


txres=True
txres@txFontHeightF=0.012
gsn_text_ndc(wks, "[m/decade]",0.5,0.489, txres)
gsn_text_ndc(wks, "[m/decade]",0.5,0.135, txres)


pres=True
pres@gsnPanelCenter=False
pres@gsnPanelYWhiteSpacePercent=25
pres@gsnPanelLabelBar =False
pres@pmLabelBarHeightF=0.1/2.
pres@pmLabelBarWidthF=0.9
pres@lbLabelStride=2.
pres@lbLabelFontHeightF=0.012

;pres@gsnPanelFigureStrings =(/"a-spring","b-spring","c-summer","d-summer","e-autumn","f-autumn","g-winter","h-winter"/);,"c","d"
;pres@gsnPanelFigureStrings =(/"a","b","c"/)
;pres@amJust="TopLeft"
gsn_panel(wks,plot1,(/2,2/),pres)
;delete(wks)



  end
