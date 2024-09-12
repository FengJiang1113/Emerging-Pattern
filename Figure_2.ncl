load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin

syear    = 1870
eyear    = 2022
ymdStrt1 = syear*100+1
ymdLast1 = eyear*100+12
ymdbetween=2000*100+1
mm=(eyear-syear+1)*12
mm1=(2018-2000)*12;(2000-1980)*12;

opendir="/data0/fjiang/data/ocean/"
resultdir="/data0/fjiang/scripts/contrast/"


a= addfile(opendir+"HadISST_sst.nc","r");sst.mnmean.nc
time1      = a->time
ymd1      = cd_calendar(time1, -1)
iStrt3    = ind(ymd1.eq.ymdStrt1)
iLast3    = ind(ymd1.eq.ymdLast1)
ibe=ind(ymd1.eq.ymdbetween)
delete(time1)
delete(ymd1)
ha   = a->sst(iStrt3:iLast3,::-1,:)
ha:=lonPivot(ha,0.5)
;ha:=lonPivot(ha,180.5)

hor= where(abs(ha).ge.50,ha@_FillValue,ha)
hor= where(ismissing(hor),0,hor)
copy_VarMeta(ha, hor)
h=hor

sstannual=month_to_annual(h,1)
printVarSummary(sstannual)


syear1=1980
syear2=1958
year1    = ispan(syear1,eyear,1)
year2    = ispan(syear2,eyear,1)


rc1    = regCoef_n(year1,sstannual(syear1-syear:eyear-syear,:,:),0,0)
rc2    = regCoef_n(year2,sstannual(syear2-syear:eyear-syear,:,:),0,0)


  rc1:= rc1*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
  rc2:= rc2*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)

copy_VarMeta(sstannual(0,:,:),rc1)
copy_VarMeta(sstannual(0,:,:),rc2)

rc1_all=new((/111,180,360/),float)
rc2_all=new((/90,180,360/),float)
patterncor1=new((/153/),float)
patterncor1_a=new((/153/),float)
patterncor1_b=new((/153/),float)
EW1=new((/153/),float)
NS1=new((/153/),float)

patterncor2=new((/153/),float)
patterncor2_a=new((/153/),float)
patterncor2_b=new((/153/),float)
EW2=new((/153/),float)
NS2=new((/153/),float)

lats=-60
late=60
lons=120
lone=290
xx1=year1
printVarSummary(xx1)
do i=0,110
  rc1_all(i,:,:)=regCoef_n(year1,sstannual(0+i:42+i,:,:),0,0)
  copy_VarMeta(sstannual(0,:,:),rc1_all(i,:,:))
  patterncor1(i+42)=pattern_cor(rc1({lats:late},{lons:lone}),rc1_all(i,{lats:late},{lons:lone}), 1.0, 0)
  patterncor1_a(i+42)=pattern_cor(rc1({lats+30:late-30},{lons:lone}),rc1_all(i,{lats+30:late-30},{lons:lone}), 1.0, 0)
  patterncor1_b(i+42)=pattern_cor(rc1({lats+50:late-50},{lons:lone}),rc1_all(i,{lats+50:late-50},{lons:lone}), 1.0, 0)
  EW1(i+42)=10*(wgt_areaave_Wrap(rc1_all(i,{-5:5},{140:170}), 1.0,1.0,0)-wgt_areaave_Wrap(rc1_all(i,{-5:5},{190:270}), 1.0,1.0,0))
  NS1(i+42)=10*(0.5*(wgt_areaave_Wrap(rc1_all(i,{3:9},{190:240}), 1.0,1.0,0)+wgt_areaave_Wrap(rc1_all(i,{-9:-3},{190:240}), 1.0,1.0,0))-wgt_areaave_Wrap(rc1_all(i,{-5:5},{190:240}), 1.0,1.0,0))
end do

xx2=year2
do i=0,88
  rc2_all(i,:,:)=regCoef_n(year2,sstannual(0+i:64+i,:,:),0,0)
  copy_VarMeta(sstannual(0,:,:),rc2_all(i,:,:))
  patterncor2(i+64)=pattern_cor(rc2({lats:late},{lons:lone}),rc2_all(i,{lats:late},{lons:lone}), 1.0, 0)
  patterncor2_a(i+64)=pattern_cor(rc2({lats+30:late-30},{lons:lone}),rc2_all(i,{lats+30:late-30},{lons:lone}), 1.0, 0)
  patterncor2_b(i+64)=pattern_cor(rc2({lats+50:late-50},{lons:lone}),rc2_all(i,{lats+50:late-50},{lons:lone}), 1.0, 0)
  EW2(i+64)=10*(wgt_areaave_Wrap(rc2_all(i,{-5:5},{140:170}), 1.0,1.0,0)-wgt_areaave_Wrap(rc2_all(i,{-5:5},{190:270}), 1.0,1.0,0))
  NS2(i+64)=10*(0.5*(wgt_areaave_Wrap(rc2_all(i,{3:9},{190:240}), 1.0,1.0,0)+wgt_areaave_Wrap(rc2_all(i,{-9:-3},{190:240}), 1.0,1.0,0))-wgt_areaave_Wrap(rc2_all(i,{-5:5},{190:240}), 1.0,1.0,0))
end do


opendir="/data0/fjiang/data/ocean/oras5/postprocess/"
a = addfile(opendir+"sst.remap.nc", "r")


syear1=1980
syear2=1958

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

sst       =a->sosstsst(:,:,:)

sstannual1 =month_to_annual(sst(iStrt1:iLast1,:,:),1)
sstannual2 =month_to_annual(sst(iStrt2:iLast1,:,:),1)

printVarSummary(sst)


dims=dimsizes(sstannual1)
rc1_oras5    = regCoef_n(year1,sstannual1,0,0)
rc2_oras5    = regCoef_n(year2,sstannual2,0,0)
rc1_oras5:= rc1_oras5*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
rc2_oras5:= rc2_oras5*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)

copy_VarMeta(sstannual1(0,:,:),rc1_oras5)
copy_VarMeta(sstannual1(0,:,:),rc2_oras5)

printVarSummary(rc1({lats:late},{lons:lone}))
printVarSummary(rc1_oras5({lats:late},{lons:lone}))
patterncor1_oras5=pattern_cor(rc1({lats:late},{lons:lone}),rc1_oras5({lats+1:late},{lons+1.:lone}), 1.0, 0)
patterncor1_a_oras5=pattern_cor(rc1({lats+30:late-30},{lons:lone}),rc1_oras5({lats+31:late-30},{lons+0.5:lone}), 1.0, 0)
patterncor1_b_oras5=pattern_cor(rc1({lats+50:late-50},{lons:lone}),rc1_oras5({lats+51:late-50},{lons+0.5:lone}), 1.0, 0)

patterncor2_oras5=pattern_cor(rc2({lats:late},{lons:lone}),rc2_oras5({lats+1:late},{lons+1.:lone}), 1.0, 0)
patterncor2_a_oras5=pattern_cor(rc2({lats+30:late-30},{lons:lone}),rc2_oras5({lats+31:late-30},{lons+0.5:lone}), 1.0, 0)
patterncor2_b_oras5=pattern_cor(rc2({lats+50:late-50},{lons:lone}),rc2_oras5({lats+51:late-50},{lons+0.5:lone}), 1.0, 0)



;;  ;confi=rc
;;;;;;



ipo_obs  =wgt_areaave_Wrap(h(:,{-10:10},{170:270}), 1.0,1.0,0)-0.5*(wgt_areaave_Wrap(h(:,{25:45},{140:215}), 1.0,1.0,0)+wgt_areaave_Wrap(h(:,{-50:-15},{150:200}), 1.0,1.0,0))
fcbt =1/(12.*13)
cft1 =ezfftf_n(ipo_obs, 0)
fcft =fspan(0, 0.5,mm/2)
ifcft1=ind(fcft.gt.fcbt)  ;;decadaal
cft1(:,ifcft1)=0.0
ipo_obs:=ezfftb_n(cft1,cft1@xbar,1);



ihp       = 0
sigma     = 1
nWgt      = 157;;;78                          ; loose 156 months each end
fca       = 1./156                       ; decadal- 13*12
wgtt      = filwgts_lanczos (nWgt, ihp, fca, -999., sigma )
;nino34_obs=wgt_runave (    nino34_obs, wgtt, 0 )   ; 10 year
  ;ipo_obs=wgt_runave (       ipo_obs, wgtt, 0 )   ; 10 year


rc1_ipo   = regCoef_n(ipo_obs,h,0,0)
copy_VarMeta(h(0,:,:),rc1_ipo)

rc1_avg=wgt_areaave_Wrap(rc1_ipo({-5:5},{140:170}), 1.0,1.0,0)-wgt_areaave_Wrap(rc1_ipo({-5:5},{190:270}), 1.0,1.0,0)
ew_avg =wgt_areaave_Wrap(h(:,{-5:5},{140:170}), 1.0,1.0,0)-wgt_areaave_Wrap(h(:,{-5:5},{190:270}), 1.0,1.0,0)



dims=dimsizes(h)
confi1=student_t(rc1_ipo@tval,rc1_ipo@nptxy-2)
confi1:=onedtond(confi1,(/dims(1),dims(2)/))
rc1_ipo@long_name = "Trend"
rc1_ipo=-rc1_ipo
copy_VarMeta(h(0,:,:),rc1_ipo)
copy_VarMeta(rc1_ipo,confi1)

ipo_ew=rc1_avg*month_to_annual(ipo_obs,1)
ew=month_to_annual(ew_avg,1)

rc1_all_ipo =new((/153/),float)
rc2_all_ipo =new((/153/),float)
rc1_raw=new((/153/),float)
rc2_raw=new((/153/),float)
syear1=1980
syear2=1958
year1:=ispan(syear1,eyear,1)
year2:=ispan(syear2,eyear,1)

do i=0,110
  rc1_all_ipo(i+42)=10*regCoef_n(year1,ipo_ew(0+i:42+i),0,0)
  rc1_raw(i+42)=10*regCoef_n(year1,ew(0+i:42+i),0,0)

end do

do i=0,88
  rc2_all_ipo(i+64)=10*regCoef_n(year2,ipo_ew(0+i:64+i),0,0)
  rc2_raw(i+64)=10*regCoef_n(year2,ew(0+i:64+i),0,0)

end do





;;;;;
sstannual=month_to_annual(h,1)


syear1=1980
syear2=1958
year1    = ispan(syear1,eyear,1)
year2    = ispan(syear2,eyear,1)


rc1_map    = regCoef_n(year1,sstannual(syear1-syear:eyear-syear,:,:),0,0)

copy_VarMeta(sstannual(0,:,:),rc1_map)
dims=dimsizes(sstannual)
confi1_map=student_t(rc1_map@tval,rc1_map@nptxy-2)
confi1_map:=onedtond(confi1,(/dims(1),dims(2)/))
rc1_map@long_name = "Trend"
rc1_map:= rc1*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
copy_VarMeta(sstannual(0,:,:),rc1_map)
copy_VarMeta(rc1_map,confi1_map)


rc1_all_map=new((/111,180,360/),float)
confi1_all_map=new((/111,180,360/),float)


lats=-60
late=60
lons=120
lone=290
xx1=year1
printVarSummary(xx1)
do i=0,110
  rc1_temp=regCoef_n(year1,sstannual(0+i:42+i,:,:),0,0)
  copy_VarMeta(sstannual(0,:,:),rc1_temp)
  rc1_all_map(i,:,:)=rc1_temp
  confi1_temp=student_t(rc1_temp@tval,rc1_temp@nptxy-2)
  confi1_all_map(i,:,:)=onedtond(confi1_temp,(/dims(1),dims(2)/))
  copy_VarMeta(rc1,confi1_all_map(i,:,:))
end do




;
  wks = gsn_open_wks("eps",resultdir+"revise-01/r01-F2");1979-2018")
   gsn_define_colormap(wks,"MPL_RdBu")
   ;gsn_define_colormap(wks,"MPL_RdBu");"CBR_drywet");
  gsn_reverse_colormap(wks);
  cmap = gsn_retrieve_colormap(wks)
res = True

res@tiMainFont=22
res@tiXAxisString = "Pattern Correlation"
res@tiYAxisString = ""
res@tmXBLabelStride=2.
;---Just create plots for now, so we can get info from each one.
res@gsnDraw  = False
res@gsnFrame = False
res@tmXBMajorLengthF    =0.003
res@tmYLMajorLengthF    =0.003
res@tmXBMinorLengthF    =0.0;05
res@tmXBMajorOutwardLengthF=0.003
res@tmYLMajorOutwardLengthF=0.003
res@tmYROn=False
res@tmYLOn=True
res@tmYLLabelsOn=True
res@tmYLMode              = "Explicit"         ; explicit labels
res@tmYRMode              = "Explicit"         ; explicit labels
res@trYMaxF=25;
res@tmYLMinorOn=False
res@tmYRMinorOn=False
res@tmYLValues=5.*ispan(0,20,1)
res@tmYLLabels=5.*ispan(0,20,1)
res@tmYUseLeft=False
res@tmYRLabelsOn=False
res@tmYRValues=10.*ispan(0,20,1)
res@tmYRLabels=1.*ispan(0,20,1)
res@gsnHistogramBinIntervals = 0.2*ispan(-5,5,1);(/-1.2,-0.8,-0.4,0.,0.4,0.8,1.2/)
res@tiMainFontHeightF=0.012
res@tiMainOffsetYF=-0.01;1
res@gsnHistogramBinWidth=0.05
res@tiYAxisFontHeightF=0.01
res@tiXAxisFontHeightF=0.01
res@tmXBLabelFontHeightF=0.01
res@tmYLLabelFontHeightF=0.01
;---Generate the plots so we can get the maximum number in the bins.
res@gsFillColor=20
plot1=new(4,graphic)
patterncor1:=patterncor1_a
patterncor2:=patterncor2_a

patterncor1=where(patterncor1.ge.1.,0.99,patterncor1)
patterncor2=where(patterncor2.ge.1.,0.99,patterncor2)


res@tiMainString  = "(a) 43-year trend distribution for pattern correlation"


res@vpXF          = 0.06    ; Move to left a bit.
res@vpYF          = 0.85
res@vpWidthF      = 0.4
res@vpHeightF     = 0.18
plot1(0) = gsn_histogram(wks,patterncor1(42::),res)
res@tiMainString  = "(b) 65-year trend distribution for pattern correlation"
res@vpXF=0.54

plot1(1) = gsn_histogram(wks,patterncor2(42::),res)



res@tiMainString  = ""
res1=True
res1@tiMainString  = ""
plot1_a=new(111,graphic)
plot2_a=new(111,graphic)

res1@tmXTBorderOn=False
res1@tmYRBorderOn=False
res1@tmXBBorderOn=False
res1@tmYLBorderOn=False
res1@tmYLOn  =False
res1@tmYROn  =False
res1@tmXBOn  =False
res1@tmXTOn  =False
res@tmXTBorderOn=False
res@tmYRBorderOn=False
res@tmXBBorderOn=False
res@tmYLBorderOn=False
res@tmYLOn  =False
res@tmYROn  =False
res@tmXBOn  =False
res@tmXTOn  =False
do i=0,110
res@gsFillColor=5+i
res@gsEdgeColor=5+i
res1@gsFillColor=5+i
res1@gsEdgeColor=5+i
plot1_a(i) = gsn_histogram(wks,patterncor1(i+42::),res)
plot2_a(i) = gsn_histogram(wks,patterncor2(i+42::),res)
end do

lnres=True
lnres@gsLineDashPattern=0
lnres@gsLineThicknessF=1.
lnres@gsLineColor="gray50"
lnres@gsLineDashPattern=2
lnres@gsLineThicknessF=1.


line1=gsn_add_polyline(wks, plot1(0), (/0.5,0.5/),(/-100,100/), lnres)
line2=gsn_add_polyline(wks, plot1(1), (/0.5,0.5/),(/-100,100/), lnres)



x1_1=dim_avg_n_Wrap(patterncor1,0)
x1_1_a=dim_avg_n_Wrap(patterncor1(42:104),0)
x1_1_b=dim_avg_n_Wrap(patterncor1(105::),0)

x1_2=dim_avg_n_Wrap(patterncor2,0)
x1_2_a=dim_avg_n_Wrap(patterncor2(64:104),0)
x1_2_b=dim_avg_n_Wrap(patterncor2(105::),0)

v1_1=sqrt(dim_variance_n_Wrap(patterncor1,0))
v1_1_a=sqrt(dim_variance_n_Wrap(patterncor1(42:104),0))
v1_1_b=sqrt(dim_variance_n_Wrap(patterncor1(105::),0))

v1_2=sqrt(dim_variance_n_Wrap(patterncor2,0))
v1_2_a=sqrt(dim_variance_n_Wrap(patterncor2(64:104),0))
v1_2_b=sqrt(dim_variance_n_Wrap(patterncor2(105::),0))


x2_1=dim_avg_n_Wrap(EW1,0)
x2_1_a=dim_avg_n_Wrap(EW1(42:104),0)
x2_1_b=dim_avg_n_Wrap(EW1(105::),0)

x2_2  =dim_avg_n_Wrap(EW2,0)
x2_2_a=dim_avg_n_Wrap(EW2(64:104),0)
x2_2_b=dim_avg_n_Wrap(EW2(105::),0)

v2_1  =sqrt(dim_variance_n_Wrap(EW1,0))
v2_1_a=sqrt(dim_variance_n_Wrap(EW1(42:104),0))
v2_1_b=sqrt(dim_variance_n_Wrap(EW1(105::),0))

v2_2  =sqrt(dim_variance_n_Wrap(EW2,0))
v2_2_a=sqrt(dim_variance_n_Wrap(EW2(64:104),0))
v2_2_b=sqrt(dim_variance_n_Wrap(EW2(105::),0))


opt=True
opt@N=100
nd1_1   = genNormalDist(x1_1,v1_1,     opt)
nd1_1_a = genNormalDist(x1_1_a,v1_1_a, opt)
nd1_1_b = genNormalDist(x1_1_b,v1_1_b, opt)

nd1_2   = genNormalDist(x1_2,v1_2,     opt)
nd1_2_a = genNormalDist(x1_2_a,v1_2_a, opt)
nd1_2_b = genNormalDist(x1_2_b,v1_2_b, opt)

pdf_nd1_1   =nd1_1   /sum(nd1_1   )*100
pdf_nd1_1_a =nd1_1_a /sum(nd1_1_a )*100
pdf_nd1_1_b =nd1_1_b /sum(nd1_1_b )*100
pdf_nd1_2   =nd1_2   /sum(nd1_2   )*100
pdf_nd1_2_a =nd1_2_a /sum(nd1_2_a )*100
pdf_nd1_2_b =nd1_2_b /sum(nd1_2_b )*100

nd2_1   = genNormalDist(x2_1,  v2_1,     opt)
nd2_1_a = genNormalDist(x2_1_a,v2_1_a, opt)
nd2_1_b = genNormalDist(x2_1_b,v2_1_b, opt)

nd2_2   = genNormalDist(x2_2,  v2_2,     opt)
nd2_2_a = genNormalDist(x2_2_a,v2_2_a, opt)
nd2_2_b = genNormalDist(x2_2_b,v2_2_b, opt)

pdf_nd2_1   =nd2_1   /sum(nd2_1   )*100
pdf_nd2_1_a =nd2_1_a /sum(nd2_1_a )*100
pdf_nd2_1_b =nd2_1_b /sum(nd2_1_b )*100
pdf_nd2_2   =nd2_2   /sum(nd2_2   )*100
pdf_nd2_2_a =nd2_2_a /sum(nd2_2_a )*100
pdf_nd2_2_b =nd2_2_b /sum(nd2_2_b )*100



  resh                         = True
  resh@gsnDraw                 = False
  resh@gsnFrame                = False
  ;resh@tmXBLabelStride          = 2
  resh@tiYAxisFontHeightF=0.012
  resh@tiXAxisFontHeightF=0.012
  ;resh@tmXBLabelStride          = 2
  ;resh@tmYLLabelStride          = 2

  resh@trXMaxF=1.
  resh@trXMinF=-1.
  resh@trYMaxF=4.
  resh@tmXTBorderOn=False
  resh@tmYRBorderOn=False
  resh@tmXBBorderOn=False
  resh@tmYLBorderOn=False
;  resh@tmXBOn  =True
;  resh@tmYLLabelFontHeightF=0.01
;  resh@tmXBLabelFontHeightF=0.01
;
;
;  resh@tmXBMinorOn=False
;  resh@tmYLMinorOn=False
;
  resh@tmYLOn  =False
  resh@tmYROn  =False
  resh@tmXBOn  =False
  resh@tmXTOn  =False
;  resh@tmYLValues   = ispan(0,4,1)   ;
;  resh@tmXBValues   = 0.4*ispan(-3,3,1)   ;
;  resh@tmXBLabels   = (/"-1.2","-0.8","-0.4"," 0.0"," 0.4"," 0.8"," 1.2"/)
;  resh@tmYLLabels   = (/"0.0","1.0","2.0","3.0","4.0"/)
;
;  resh@tmXBMajorLengthF    =0.003
;  resh@tmYLMajorLengthF    =0.003
;  resh@tmXBMinorLengthF    =0.0;05
;  resh@tmXBMajorOutwardLengthF=0.003
;  resh@tmYLMajorOutwardLengthF=0.003
  resh@xyCurveDrawOrder="PostDraw"
;  resh@tiMainString    = "";"CO2: N="+nx+"  nBands="+resh@gsnHistogramNumberOfBins
;  ;plot2 = gsn_histogram(wks,x,resh)         ; create histogram with 20 bins
;  resh@tiYAxisString="Probability Density"

  resh@xyLineThicknessF=4.
  ;resh@gsnXRefLine=0.
  ;resh@gsnXRefLineColor="gray"
  plot2=new(4,graphic)
  resh@xyLineColor="deepskyblue3"
  plot2(0) = gsn_csm_xy (wks,(nd1_1_a@x+1.2)/2.4,pdf_nd1_1_a*10,resh) ; create plot
  plot2(1) = gsn_csm_xy (wks,(nd1_2_a@x+1.2)/2.4,pdf_nd1_2_a*10,resh) ; create plot
  ;plot2(2) = gsn_csm_xy (wks,(nd2_1_a@x+0.18)/0.36,pdf_nd2_1_a*20,resh) ; create plot
  ;plot2(3) = gsn_csm_xy (wks,(nd2_2_a@x+0.18)/0.36,pdf_nd2_2_a*20,resh) ; create plot


  plot3=new(4,graphic)
  resh@xyLineColor="orangered3"
  plot3(0) = gsn_csm_xy (wks,(nd1_1_b@x+1.2)/2.4,pdf_nd1_1_b*10,resh) ; create plot
  plot3(1) = gsn_csm_xy (wks,(nd1_2_b@x+1.2)/2.4,pdf_nd1_2_b*10,resh) ; create plot
  ;plot3(2) = gsn_csm_xy (wks,(nd2_1_b@x+0.18)/0.36,pdf_nd2_1_b*20,resh) ; create plot
  ;plot3(3) = gsn_csm_xy (wks,(nd2_2_b@x+0.18)/0.36,pdf_nd2_2_b*20,resh) ; create plot




  do i=0,110
  overlay(plot1(0),plot1_a(i))
  overlay(plot1(1),plot2_a(i))
  ;overlay(plot1(2),plot3_a(i))
  ;overlay(plot1(3),plot4_a(i))

  end do





labels=ispan(1870+42,2022,1)
levels=ispan(1870+42,2022,1)

firstdx = 5                                    ; starting index
stride  = 1
nlevels = dimsizes(levels)                                ; interval between colors
colors  = ispan(firstdx,firstdx+stride*nlevels,stride)
ncolors = dimsizes(colors)
labelbarcolors        = new((/ncolors,3/),float)
ncolors=111
do n = 0, ncolors-1
  labelbarcolors(n,:) = cmap(colors(n),:)
end do
;labels:=ispan(1912,2022,10)+8
;levels:=ispan(1912,2022,10)+8


    lbres = True
    lbres@lbPerimOn            = False               ; no label bar box
    lbres@lbOrientation        = "Horizontal"        ; orientation
    lbres@vpWidthF             = 0.77                ; size
    lbres@vpHeightF            = 0.07
    lbres@lbLabelFontHeightF   = 0.008            ; label font height
    lbres@lbLabelAlignment     = "InteriorEdges"     ; where to label
    lbres@lbBoxLinesOn=False
    lbres@lbMonoFillPattern    = True                ; fill sold
    lbres@lbFillColors         = labelbarcolors     ; must be RGB triplets
    lbres@lbLabelStride=4.
    gsn_labelbar_ndc (wks,nlevels,labels,0.1,0.63,lbres)

    txres                       = True
    ;txres@txFont                = "helvetica-bold"
    txres@txFontHeightF         = 0.011
    gsn_text_ndc (wks,"End year of the trend",0.5,0.56,txres)
    txres@txAngleF=90
    ;gsn_text_ndc (wks,"Probability Density (%)",0.98,0.76,txres)
    ;gsn_text_ndc (wks,"Probability Density (%)",0.96,0.376,txres)
    gsn_text_ndc (wks,"Frequency",0.01,0.76,txres)
    ;gsn_text_ndc (wks,"Frequecy",0.04,0.376,txres)
    txres@txAngleF=0
    gsn_text_ndc (wks,"[~S~o~N~C/decade]",0.96,0.3,txres)


    ;;;;;;;;;;;;;;;;;;;;legend
    lgres=True
    amres=True
    lgres@lgLineColors=(/"orangered3","deepskyblue3"/)
    lgres@lgLeftMarginF=0.1
    lgres@lgRightMarginF=0.1

    lgres@lgDashIndexes=(/0,0/)
    lgres@lgLineThicknesses=(/4.,4./)
    lgres@lgMonoLineThickness=False
    lgres@vpWidthF=0.162
    lgres@vpHeightF=0.05
    lgres@lgLabelFontHeightF=0.12
    lgres@lgLineDashSegLenF=0.04
    lgres@lgPerimThicknessF=2.
    lgres@lgPerimOn=False
    lgres@lgPerimFill="SolidFill"
    lgres@lgPerimFillColor="white"
    ;lgres@lgLabelOffsetF=-0.1
    lgres@lbBoxMajorExtentF=0.3
    lbid1=gsn_create_legend(wks,2,(/"   End year in 1975-2022 (48)","   End year in 1942-1974 (63)"/),lgres)
    lbid3=gsn_create_legend(wks,2,(/"   End year in 1975-2022 (48)","   End year in 1942-1974 (63)"/),lgres)

    amres@amParallelPosF=-0.29
    amres@amOrthogonalPosF=-0.3
    ;annoid1=gsn_add_annotation(plot1(0),lbid1,amres)
    ;annoid3=gsn_add_annotation(plot1(2),lbid3,amres)

    lbid2=gsn_create_legend(wks,2,(/"  End year in 1975-2022 (48)","  End year in 1934-1974 (41)"/),lgres)
    lbid4=gsn_create_legend(wks,2,(/"  End year in 1975-2022 (48)","  End year in 1934-1974 (41)"/),lgres)

    amres@amParallelPosF=-0.29
    amres@amOrthogonalPosF=-0.3
    ;annoid2=gsn_add_annotation(plot1(1),lbid2,amres)
    ;annoid4=gsn_add_annotation(plot1(3),lbid4,amres)




draw(plot1(0))
draw(plot1(1))


delete(res1)
res1                 = True
res1@gsnDraw         = False
res1@gsnFrame        = False
res1@vpWidthF      = 0.4
res1@vpHeightF     = 0.18
res1@trXMinF   =  0.;1870.
res1@trXMaxF   =  mm-1;2022.
res1@trYMaxF   = 1.
res1@trYMinF   = -1.
res1@tiYAxisString = " "
res1@tiXAxisString = " "
res1@tiMainString  = " "
res1@tmXBMajorThicknessF     = 1.5
res1@tmYLMajorThicknessF     = 1.5
res1@tmXBMinorThicknessF     = 1.0
res1@tmYLMinorThicknessF     = 1.0
res1@tmYLMinorOn=False
res1@tmYLValues    = 0.5*ispan(-5,5,1)
res1@tmYLLabels    = 0.5*ispan(-5,5,1);
res1@tmXBValues    = ispan(10*12,mm-1,20*12); (/"1960","1965","1970","1975","1980","1985","1985","1990","1995","2000","2005","2010","2015","2020"/)
res1@tmXBLabels    = ispan(1880,2022,20); (/"1960","1965","1970","1975","1980","1985","1985","1990","1995","2000","2005","2010","2015","2020"/)

res1@gsnYRefLine=0.
res1@tmXBMajorLengthF=0.006
res1@tmYLMajorLengthF=0.004
res1@tmXBMajorOutwardLengthF=0.006
res1@tmYLMajorOutwardLengthF=0.004
res1@tmXBMinorLengthF        = 0.004
res1@tmYLMinorLengthF        = 0.004
res1@tmXBMinorOn=True
;res@tmXBMinorValues=ispan(1958,2022,1)

res1@gsnLeftString=" "
res1@gsnLeftStringFontHeightF= 0.012
res1@gsnLeftStringOrthogonalPosF=0.00;5
res1@tmXBOn = True
;res@tmYLOn = True
;res@tmYLLabelsOn = True
;res@tmYLMinorOn = True
;res@tmYLLabelFontColor = "black"
res1@tiXAxisString         = "";End year of the trend"
res1@tiXAxisFontHeightF=0.01
res1@tiYAxisFontColor      = "black"
res1@tmYLOn    = True
res1@tmYLMode          = "Explicit"         ; explicit labels
res1@tmXBMode          = "Explicit"         ; explicit labels
res1@tmXBOn    = True

res1@tmYLLabelsOn=True

res1@tmXBLabelFontHeightF  = 0.008
res1@tmYLLabelFontHeightF  = 0.008
res1@tiYAxisSide           = "Left"
res1@tiYAxisFontHeightF    = 0.01
res1@xyLineColor      = "orangered3"
res1@tmXBBorderOn=True
res1@tiYAxisFontHeightF=0.012

res1@tmXTOn=False
res1@tmYROn=False
res1@gsnLeftString=""
res1@tiMainOffsetYF=-0.01
res1@tiMainFontHeightF=0.012
res1@gsnRightString=""
res1@xyLineThicknessF = 2.

res1@tiMainFont=22
res1@gsnXYBarChart = True            ; Create bar plot
res1@gsnXYBarWidth=0.2
res1@gsnBelowYRefLineBarColors="gray"
res1@gsnAboveYRefLineBarColors="gray"
res1@xyLineColor=-1
res1@tiXAxisFontHeightF=0.01
res1@tiXAxisString="Year"
res1@tiYAxisFontHeightF=0.01
res1@tiYAxisString="[~S~o~N~C]"

delete(res1@gsnXYBarChart)

res1@tmXBValues    = ispan(10,mm/12-1,20); (/"1960","1965","1970","1975","1980","1985","1985","1990","1995","2000","2005","2010","2015","2020"/)
res1@trXMinF:= 1870-1870
res1@trXMaxF:= 2022-1870
res1@trYMaxF   = 0.18
res1@trYMinF   = -0.18
res1@tmYLValues    = 0.05*ispan(-5,5,1)
res1@tmYLLabels    = 0.05*ispan(-5,5,1);
res1@xyLineColor      = "orangered3"
res1@tiMainString="(c) Zonal SST Gradient Trend and IPO Contribution"
res1@tiYAxisString="[~S~o~N~C/decade]"
res1@tiMainOffsetYF=-0.0066
res1@xyDashPattern=2.
res1@xyLineThicknessF=2.;3-s1.5
res1@tiXAxisString="End year of the trend"
res1@vpXF=0.06
res1@vpYF=0.5
plot1(2) = gsn_csm_xy(wks,ispan(0,mm/12-1,1),rc1_all_ipo,res1)
res1@tiMainString=""
res1@xyDashPattern=0.
res1@xyLineThicknessF=1.5

print(rc1_all_ipo)

plot3a = gsn_csm_xy(wks,ispan(0,mm/12-1,1),rc1_raw,res1)
res1@xyDashPattern=2.
res1@xyLineThicknessF=1.5
res1@xyLineThicknessF=2.;3-s1.5

res1@xyLineColor      = "deepskyblue3"
plot3b   = gsn_csm_xy(wks,ispan(0,mm/12-1,1),rc2_all_ipo,res1)
res1@xyDashPattern=0.
res1@xyLineThicknessF=1.5

plot3c   = gsn_csm_xy(wks,ispan(0,mm/12-1,1),rc2_raw,res1)

overlay(plot1(2),plot3a)
overlay(plot1(2),plot3b)
overlay(plot1(2),plot3c)

delete(lgres)
delete(amres)
lgres=True
amres=True
lgres@lgLineColors=(/"deepskyblue3","deepskyblue3","orangered3","orangered3"/)
lgres@lgLeftMarginF=0.1
lgres@lgRightMarginF=0.1

lgres@lgDashIndexes=(/2.,0.,2.,0./)
lgres@lgLineThicknesses=(/1.,1.,1.,1./)
lgres@lgMonoLineThickness=False
lgres@vpWidthF=0.162
lgres@vpHeightF=0.07
lgres@lgLabelFontHeightF=0.12
lgres@lgLineDashSegLenF=0.04
lgres@lgPerimThicknessF=2.
lgres@lgPerimOn=False
lgres@lgPerimFill="SolidFill"
lgres@lgPerimFillColor="white"
;lgres@lgLabelOffsetF=-0.1
lgres@lbBoxMajorExtentF=0.3
lbid1=gsn_create_legend(wks,4,(/"   IPO Contribution to 65-year trend","   65-year trend","   IPO Contribution to 43-year trend","   43-year trend"/),lgres)

amres@amParallelPosF=-0.28
amres@amOrthogonalPosF=-0.28
annoid1=gsn_add_annotation(plot1(2),lbid1,amres)
;annoid2=gsn_add_annotation(plot1(3),lbid2,amres)

draw(plot1(2))


delete(res)
res = True
res@gsnDraw   = False
res@gsnFrame  = False
res@gsnAddCyclic      = True
res@mpMinLonF  = 0 ;lon2d(0,0)
res@mpCenterLonF=180
res@mpMaxLonF = 360 ;lon2d(0,ilon-1)
res@mpFillOn=True
res@mpLandFillColor="gray70"
res@mpOceanFillColor=-1
res@mpFillDrawOrder="PostDraw"
res@tmXBMinorOn=False
;res@tmYLLabels =(/"30S","20S","10S","0","10N","20N","30N"/)
res@tmYLMinorOn=False
res@cnInfoLabelOn=False
res@mpGeophysicalLineThicknessF=1
res@mpGeophysicalLineColor="gray70"
;res@pmTickMarkDisplayMode = "Always"
res@tmXBLabelFontHeightF  = 0.007
res@tmYLLabelFontHeightF  = 0.007

res@cnFillOn    = True
res@cnLinesOn = False

res@tmXBMajorLengthF=0.004
res@tmYLMajorLengthF=0.004
res@tmXBMajorOutwardLengthF=0.004
res@tmYLMajorOutwardLengthF=0.004
res@mpProjection="Robinson"
res@mpPerimOn=False
res@cnLevelSelectionMode = "ExplicitLevels"  ; use explicit levels
;res@cnFillColors=(/111,100,96,87,77,67,64,54,44,34,30,20/);
;res@cnFillColors=(/10,15,20,25,30,35,40,45,50,57,0,0,73,80,85,90,95,100,105,110,115,120/)
res@cnFillColors=(/12,20,30,40,50,57,0,0,73,80,90,96,105,115/)

res@gsnLeftStringFontHeightF=0.01
res@gsnLeftStringOrthogonalPosF=0.02;0.005
res@gsnRightStringFontHeightF=0.012
res@gsnRightStringOrthogonalPosF=-0.003;5;0.005
res@cnLineLabelsOn=False
res@gsnRightString="";(~S~o~N~C/decade)"
res@gsnLeftString=""
res@vpHeightF= 0.2            ; change aspect ratio of plot
res@cnLevels             = 0.02*ispan(-12,12,2 )

res@lbLabelBarOn=False
res@lbLabelBarOn=True
res@lbOrientation="Vertical"
;res@lbBoxEndCapStyle="TriangleBothEnds"
res@pmLabelBarHeightF=0.16
res@pmLabelBarWidthF=0.03
res@lbLabelStride=4.
res@lbBoxLinesOn=True
res@pmLabelBarOrthogonalPosF=0.05
res@lbLabelFontHeightF=0.01
;res@pmLabelBarOrthogonalPosF=0.2
;res@pmLabelBarParallelPosF=-0.11
res@tiMainOffsetYF=-0.007
res@tiMainFontHeightF=0.012
res@vpXF=0.36
res@vpYF=0.5
res@gsnLeftStringOrthogonalPosF=-0.02;0.005
rc1_all_map_10=rc1_all_map*10
copy_VarMeta(rc1_all_map,rc1_all_map_10)
res@tiMainString="(d)  SST Trend 1970-2012 (near-zero IPO Contribution)  "
plot1(3) = gsn_csm_contour_map(wks,rc1_all_map_10(142-42,:,:),res)

sres1=True
sres1@gsnDraw=False
sres1@gsnFrame=False
sres1@cnFillOn=True
sres1@cnLinesOn=False
sres1@cnMonoFillPattern=False
sres1@gsnLeftString=""
sres1@gsnRightString=""
sres1@cnLevelSelectionMode="ExplicitLevels"
sres1@cnLevels=(/0,0.05/)
sres1@cnFillPatterns=(/-1,17,-1/)
sres1@cnFillColors=(/"black","black","black"/)
sres1@cnLineLabelsOn=False
sres1@cnInfoLabelOn=False
sres1@lbLabelBarOn=False
;
plots=gsn_csm_contour(wks,confi1_all_map(141-42,:,:),sres1)

overlay(plot1(3),plots)
draw(plot1(3))




   end
