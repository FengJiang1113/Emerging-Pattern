load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
;;;ipo_obs
syear    = 1958
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

hor= where(abs(ha).ge.50,ha@_FillValue,ha)
hor= where(ismissing(hor),0,hor)
hor=rmMonAnnCycTLL(hor)
copy_VarMeta(ha, hor)

h=hor

ipo_obs  =wgt_areaave_Wrap(h(:,{-10:10},{170:270}), 1.0,1.0,0)-0.5*(wgt_areaave_Wrap(h(:,{25:45},{140:215}), 1.0,1.0,0)+wgt_areaave_Wrap(h(:,{-50:-15},{150:200}), 1.0,1.0,0))
fcbt =1/(12.*13)
cft1 =ezfftf_n(ipo_obs, 0)
fcft =fspan(0, 0.5,mm/2)
ifcft1=ind(fcft.gt.fcbt)  ;;decadaal
cft1(:,ifcft1)=0.0
ipo_obs:=-ezfftb_n(cft1,cft1@xbar,1);
printVarSummary(ipo_obs)
ipo_obs:=dim_standardize_n_Wrap(ipo_obs,1,0)
rc_ipo=120*regCoef_n(ispan(0,43*12-1,1),ipo_obs(22*12::),0,0)
print(rc_ipo)

rc1_ipo= regCoef_n(ipo_obs,h,0,0)
copy_VarMeta(h(0,:,:),rc1_ipo)
;;;;;


resultdir = "/data0/fjiang/scripts/contrast/figuredata/hb/"
fout1 = addfile(resultdir+"uatc.nc","r")

latS   = -21
latN   = 21
lonL   = 120
lonR   = 290
time1      = fout1->time_counter
ymd1      = cd_calendar(time1, -1)
syear1=1958
eyear1=1977
ymdStrt1 = syear1*100+1
ymdLast1 = eyear1*100+12

syear2=1980
eyear2=1999
ymdStrt2 = syear2*100+1
ymdLast2 = eyear2*100+12

syear3=2003
eyear3=2022
ymdStrt3 = syear3*100+1
ymdLast3 = eyear3*100+12

iStrt1    = ind(ymd1.eq.ymdStrt1)
iStrt2    = ind(ymd1.eq.ymdStrt2)
iStrt3    = ind(ymd1.eq.ymdStrt2)
iLast1    = ind(ymd1.eq.ymdLast1)
iLast2    = ind(ymd1.eq.ymdLast2)
iLast3    = ind(ymd1.eq.ymdLast3)


resultdir = "/data0/fjiang/scripts/contrast/figuredata/hb/"
fout1 = addfile(resultdir+"uatc.nc","r")
uatc1 = dim_avg_n_Wrap(fout1->UATC(iStrt1:iLast1,:,:),0)
uatc2 = dim_avg_n_Wrap(fout1->UATC(iStrt2:iLast2,:,:),0)
uatc3 = dim_avg_n_Wrap(fout1->UATC(iStrt3:iLast3,:,:),0)
uatc4 = regCoef_n(ipo_obs,fout1->UATC(iStrt1:iLast3,:,:),0,0)


fout2 = addfile(resultdir+"ucta.nc","r")
ucta1 = dim_avg_n_Wrap(fout2->UCTA(iStrt1:iLast1,:,:),0)
ucta2 = dim_avg_n_Wrap(fout2->UCTA(iStrt2:iLast2,:,:),0)
ucta3 = dim_avg_n_Wrap(fout2->UCTA(iStrt3:iLast3,:,:),0)
ucta4 = regCoef_n(ipo_obs,fout2->UCTA(iStrt1:iLast3,:,:),0,0)


fout3 = addfile(resultdir+"uata.nc","r")
uata1 = dim_avg_n_Wrap(fout3->UATA(iStrt1:iLast1,:,:),0)
uata2 = dim_avg_n_Wrap(fout3->UATA(iStrt2:iLast2,:,:),0)
uata3 = dim_avg_n_Wrap(fout3->UATA(iStrt3:iLast3,:,:),0)
uata4 = regCoef_n(ipo_obs,fout3->UATA(iStrt1:iLast3,:,:),0,0)


ut1=(uatc1+ucta1+uata1)/3.
ut2=(uatc2+ucta2+uata2)/3.
ut3=(uatc3+ucta3+uata3)/3.
copy_VarMeta(uatc1,ut1)
copy_VarMeta(uatc1,ut2)
copy_VarMeta(uatc1,ut3)



fout4 = addfile(resultdir+"vatc.nc","r")
vatc1 = dim_avg_n_Wrap(fout4->VATC(iStrt1:iLast1,:,:),0)
vatc2 = dim_avg_n_Wrap(fout4->VATC(iStrt2:iLast2,:,:),0)
vatc3 = dim_avg_n_Wrap(fout4->VATC(iStrt3:iLast3,:,:),0)
vatc4 = regCoef_n(ipo_obs,fout4->VATC(iStrt1:iLast3,:,:),0,0)


fout5 = addfile(resultdir+"vcta.nc","r")
vcta1 = dim_avg_n_Wrap(fout5->VCTA(iStrt1:iLast1,:,:),0)
vcta2 = dim_avg_n_Wrap(fout5->VCTA(iStrt2:iLast2,:,:),0)
vcta3 = dim_avg_n_Wrap(fout5->VCTA(iStrt3:iLast3,:,:),0)
vcta4 = regCoef_n(ipo_obs,fout5->VCTA(iStrt1:iLast3,:,:),0,0)


fout6 = addfile(resultdir+"vata.nc","r")
vata1 = dim_avg_n_Wrap(fout6->VATA(iStrt1:iLast1,:,:),0)
vata2 = dim_avg_n_Wrap(fout6->VATA(iStrt2:iLast2,:,:),0)
vata3 = dim_avg_n_Wrap(fout6->VATA(iStrt3:iLast3,:,:),0)
vata4 = regCoef_n(ipo_obs,fout6->VATA(iStrt1:iLast3,:,:),0,0)

vt1=(vatc1+vcta1+vata1)/3.
vt2=(vatc2+vcta2+vata2)/3.
vt3=(vatc3+vcta3+vata3)/3.
copy_VarMeta(vatc1,vt1)
copy_VarMeta(vatc1,vt2)
copy_VarMeta(vatc1,vt3)


fout7 = addfile(resultdir+"watc.nc","r")
watc1 = dim_avg_n_Wrap(fout7->WATC(iStrt1:iLast1,:,:),0)
watc2 = dim_avg_n_Wrap(fout7->WATC(iStrt2:iLast2,:,:),0)
watc3 = dim_avg_n_Wrap(fout7->WATC(iStrt3:iLast3,:,:),0)
watc4 = regCoef_n(ipo_obs,fout7->WATC(iStrt1:iLast3,:,:),0,0)

fout8 = addfile(resultdir+"wcta.nc","r")
wcta1 = dim_avg_n_Wrap(fout8->WCTA(iStrt1:iLast1,:,:),0)
wcta2 = dim_avg_n_Wrap(fout8->WCTA(iStrt2:iLast2,:,:),0)
wcta3 = dim_avg_n_Wrap(fout8->WCTA(iStrt3:iLast3,:,:),0)
wcta4 = regCoef_n(ipo_obs,fout8->WCTA(iStrt1:iLast3,:,:),0,0)

fout9 = addfile(resultdir+"wata.nc","r")
wata1 = dim_avg_n_Wrap(fout9->WATA(iStrt1:iLast1,:,:),0)
wata2 = dim_avg_n_Wrap(fout9->WATA(iStrt2:iLast2,:,:),0)
wata3 = dim_avg_n_Wrap(fout9->WATA(iStrt3:iLast3,:,:),0)
wata4 = regCoef_n(ipo_obs,fout9->WATA(iStrt1:iLast3,:,:),0,0)

wt1=(watc1+wcta1+wata1)/3.
wt2=(watc2+wcta2+wata2)/3.
wt3=(watc3+wcta3+wata3)/3.
copy_VarMeta(watc1,wt1)
copy_VarMeta(watc1,wt2)
copy_VarMeta(watc1,wt3)


fout10 = addfile(resultdir+"flux.nc","r")
flux1 = dim_avg_n_Wrap(fout10->flux(iStrt1:iLast1,:,:),0)
flux2 = dim_avg_n_Wrap(fout10->flux(iStrt2:iLast2,:,:),0)
flux3 = dim_avg_n_Wrap(fout10->flux(iStrt3:iLast3,:,:),0)
flux4 = regCoef_n(ipo_obs,fout10->flux(iStrt1:iLast3,:,:),0,0)

dims=dimsizes(flux1)
data1=new((/4,dims(0), dims(1)/),float)
data2=new((/4,dims(0), dims(1)/),float)
data1(0,:,:)=86400*12*(ut3-ut1)
data1(1,:,:)=86400*12*(vt3-vt1)
data1(2,:,:)=86400*12*(wt3-wt1)
data1(3,:,:)=86400*12*(flux3-flux1)

data2(0,:,:)=86400*12*(ut3-ut2)
data2(1,:,:)=86400*12*(vt3-vt2)
data2(2,:,:)=86400*12*(wt3-wt2)
data2(3,:,:)=86400*12*(flux3-flux2)

data1_all=new((/10,dims(0), dims(1)/),float)
data2_all=new((/10,dims(0), dims(1)/),float)
data1_all(0,:,:)=86400*30*(uatc3-uatc1)/6.5
data1_all(1,:,:)=86400*30*(ucta3-ucta1)/6.5
data1_all(2,:,:)=86400*30*(uata3-uata1)/6.5
data1_all(3,:,:)=86400*30*(vatc3-vatc1)/6.5
data1_all(4,:,:)=86400*30*(vcta3-vcta1)/6.5
data1_all(5,:,:)=86400*30*(vata3-vata1)/6.5
data1_all(6,:,:)=86400*30*(watc3-watc1)/6.5
data1_all(7,:,:)=86400*30*(wcta3-wcta1)/6.5
data1_all(8,:,:)=86400*30*(wata3-wata1)/6.5
data1_all(9,:,:)=86400*30*(flux3-flux1)/6.5

data2_all(0,:,:)=86400*30*(uatc3-uatc2)
data2_all(1,:,:)=86400*30*(ucta3-ucta2)
data2_all(2,:,:)=86400*30*(uata3-uata2)
data2_all(3,:,:)=86400*30*(vatc3-vatc2)
data2_all(4,:,:)=86400*30*(vcta3-vcta2)
data2_all(5,:,:)=86400*30*(vata3-vata2)
data2_all(6,:,:)=86400*30*(watc3-watc2)
data2_all(7,:,:)=86400*30*(wcta3-wcta2)
data2_all(8,:,:)=86400*30*(wata3-wata2)
data2_all(9,:,:)=86400*30*(flux3-flux2)

copy_VarMeta(ut1,data1(0,:,:))
copy_VarMeta(ut1,data1(1,:,:))
copy_VarMeta(ut1,data1(2,:,:))
copy_VarMeta(ut1,data1(3,:,:))
copy_VarMeta(ut1,data2(0,:,:))
copy_VarMeta(ut1,data2(1,:,:))
copy_VarMeta(ut1,data2(2,:,:))
copy_VarMeta(ut1,data2(3,:,:))


do i=0,9
copy_VarMeta(ut1,data1_all(i,:,:))
copy_VarMeta(ut1,data2_all(i,:,:))
end do


data1_test1=new((/10,dims(0), dims(1)/),float)
data1_test1(0,:,:)=86400*30*uatc4*rc_ipo
data1_test1(1,:,:)=86400*30*ucta4*rc_ipo
data1_test1(2,:,:)=86400*30*uata4*rc_ipo
data1_test1(3,:,:)=86400*30*vatc4*rc_ipo
data1_test1(4,:,:)=86400*30*vcta4*rc_ipo
data1_test1(5,:,:)=86400*30*vata4*rc_ipo
data1_test1(6,:,:)=86400*30*watc4*rc_ipo
data1_test1(7,:,:)=86400*30*wcta4*rc_ipo
data1_test1(8,:,:)=86400*30*wata4*rc_ipo
data1_test1(9,:,:)=86400*30*flux4*rc_ipo

do i=0,9
copy_VarMeta(ut1,data1_test1(i,:,:))
end do



;if (any(isnan_ieee(rc_wt2))) then

do i=0,9
copy_VarMeta(ut1,data1_all(i,:,:))
copy_VarMeta(ut1,data2_all(i,:,:))
end do


delete(data1_all@_FillValue)
delete(data1_test1@_FillValue)
data1_all:=where(data1_all.gt.100.or.data1_all.lt.-100,0,data1_all)
copy_VarMeta(data1_test1,data1_all)
data1_test1:=where(data1_test1.gt.100.or.data1_test1.lt.-100,0,data1_test1)
copy_VarMeta(data1_all,data1_test1)

data=new((/2,3,10/),float) ;;;0:d&c 1:E,C,W  2:Component

data(0,0,1:9)=wgt_areaave_Wrap(data1_all(0:8,{-5:5},{150:180}), 1.0,1.0,0)
data(0,1,1:9)=wgt_areaave_Wrap(data1_all(0:8,{-5:5},{180:240}), 1.0,1.0,0)
data(0,2,1:9)=wgt_areaave_Wrap(data1_all(0:8,{-3:3},{190:270}), 1.0,1.0,0)
data(0,2,0)=sum(data(0,2,1:9))

data(1,0,1:9)=wgt_areaave_Wrap(data1_test1(0:8,{-5:5},{150:180}), 1.0,1.0,0)
data(1,1,1:9)=wgt_areaave_Wrap(data1_test1(0:8,{-5:5},{180:240}), 1.0,1.0,0)
data(1,2,1:9)=wgt_areaave_Wrap(data1_test1(0:8,{-3:3},{190:270}), 1.0,1.0,0)
data(1,2,0)=sum(data(1,2,1:9))




  wks = gsn_open_wks("eps","/data0/fjiang/scripts/contrast/v6-6");2000_original");1979-2018")
  ; gsn_define_colormap(wks,"MPL_RdBu");L_RdBu")
   ;gsn_define_colormap(wks,"MPL_RdBu");"CBR_drywet");
   ;gsn_reverse_colormap(wks);


   res=True
   res@tmXBMajorLengthF=0.015
   res@tmYLMajorLengthF=0.004
   res@tmXBMajorOutwardLengthF=0.015
   res@tmYLMajorOutwardLengthF=0.004
   res@gsnDraw=False
   res@gsnFrame=False
   res@tmYROn=False
   res@tmXTOn=False
   res@vpHeightF= 0.18          ; change aspect ratio of plot
   res@vpWidthF=0.6

   res@tmXBMode="Explicit"

   ;res@trYReverse=True
   res@trXMaxF=10+.75-1
   res@trXMinF=-0.75
   res@tmYLMode="Explicit"
   res@tmYLValues=0.012*ispan(-15,15,1)
   res@tmYLLabels=0.012*ispan(-15,15,1)
   res@tmYLLabelFontHeightF=0.008
   res@gsnYRefLineDashPattern=0
   res@gsnYRefLineColor="gray60"
   res@gsnYRefLine=0.
   res@tmYLMinorOn=False
   res@tmXBLabelFontHeightF=0.0078
   res@tmXBLabelAngleF=-35
   res@tmXBLabelJust="TopLeft"
   res@gsnXYBarChart=True
   res@gsnXYBarChartFillLineThicknessF=0.1
   res@gsnXYBarChartBarWidth=0.42
   res@gsnLeftStringOrthogonalPosF=0.005
   res@gsnLeftStringFontHeightF=0.012
   res@tiYAxisString = "[~S~o~N~C/month/decade]"                  ; yaxis
   res@tiYAxisFontHeightF=0.012
   res@tiXAxisFontHeightF=0.012

   res@tiXAxisString =  "";West Pacific                               Central Pacific                               East Pacific"
   res@gsnLeftString="";(j) Decomposition of dynamical contribution"
   ;res@xyLineColor=(/"black"/)
   ;res@gsnYRefLine=0
   md=0.28
   res@tmXBLabels=(/"Sum","UaTc","UcTa","UaTa","VaTc","VcTa","VaTa","WaTc","WcTa","WaTa"/)
   res@tmXBValues=ispan(0,9,1);(/0+md,1,2-md,3+md,4,5-md,6+md,7,8-md/)



   res@vpXF=0.1
   res@vpYF=0.27

   plot=new(2,graphic)
   print(data)
   res@tiMainFont=22
   res@tiMainString="(a) IPO-related Central-to-Eastern SST Change"
   res@tiMainFontHeightF=0.014
   res@tiMainOffsetYF=-0.005
   res@trYMaxF=0.036;5
   res@trYMinF=-0.036
   res@tiYAxisString = "[~S~o~N~C/month]"                  ; yaxis
   cmap = read_colormap_file("amwg")
   res@gsnXYBarChartColors = "salmon2";(/"SlateBlue","SkyBlue","blue","navy","green","orange","red","firebrick","orangered3"/)             ; choose colors
   ;res@gsnXYBarChartColors=(/cmap(1,:),cmap(3,:),cmap(5,:),cmap(6,:),cmap(7,:),cmap(8,:),cmap(10,:),cmap(13,:),cmap(15,:)/);(/"indianred2","salmon2","tomato3","darkseagreen","pink","purple","lightblue3","brown2","yellow"/);

   plot(0)= gsn_csm_xy(wks, ispan(0,9,1), data(1,2,:),res)
   res@gsnXYBarChartColors = "lightblue3";(/"SlateBlue","SkyBlue","blue","navy","green","orange","red","firebrick","orangered3"/)             ; choose colors


   res@tiMainString="(b) Emerging Central-to-Eastern Pacific SST Change"
   plot(1)= gsn_csm_xy(wks, ispan(0,9,1), data(0,2,:),res)


   res@gsnXYBarChartPatterns=17
   res@gsnXYBarChartColors = "black";(/"SlateBlue","SkyBlue","blue","navy","green","orange","red","firebrick","orangered3"/)             ; choose colors
   plota= gsn_csm_xy(wks, 1, data(1,2,1),res)
   overlay(plot(0),plota)
;   res@tiMainString="(e) IPO-related East Pacific SST Change"
;
;   plot(4)= gsn_csm_xy(wks, ispan(0,8,1), data(1,2,:),res)
;   res@tiYAxisString = "[~S~o~N~C/month/decade]"                  ; yaxis
;
;   res@trYMaxF=0.15;5
;   res@trYMinF=-0.15
;   res@tiMainString="(b) Emerging West Pacific SST Change"
;
;   plot(1)= gsn_csm_xy(wks, ispan(0,8,1), data(0,0,:),res)
;   res@tiMainString="(d) Emerging Central Pacific SST Change"
;
;   plot(3)= gsn_csm_xy(wks, ispan(0,8,1), data(0,1,:),res)
;   res@tiMainString="(f) Emerging East Pacific SST Change"
;
;   plot(5)= gsn_csm_xy(wks, ispan(0,8,1), data(0,2,:),res)
;
;
;   res@gsnXYBarChartColors="lightblue3"
;
;   res@gsnXYBarChartColors="gray60"
;
;   res@gsnXYBarChartColors="indianred2"











   lnres=True
   lnres@gsLineColor="gray40"
   lnres@gsLineDashPattern=2
   lnres@gsLineThicknessF=1.

   ;;dum1=gsn_add_polyline(wks,plot(0),fspan(20.5,379.5,360)-20,t_g_Height,lnres)
   dumx1=gsn_add_polyline(wks,plot(0),(/.5,.5/),(/-4,4/),lnres)
   dumx2=gsn_add_polyline(wks,plot(1),(/.5,.5/),(/-4,4/),lnres)
   dumx3=gsn_add_polyline(wks,plot(0),(/3.5,3.5/),(/-4,4/),lnres)
   dumx4=gsn_add_polyline(wks,plot(1),(/3.5,3.5/),(/-4,4/),lnres)
   dumx5=gsn_add_polyline(wks,plot(0),(/6.5,6.5/),(/-4,4/),lnres)
   dumx6=gsn_add_polyline(wks,plot(1),(/6.5,6.5/),(/-4,4/),lnres)
   ;dumx2=gsn_add_polyline(wks,plot(0),(/2.5,2.5/)+3,(/-4,4/),lnres)
   ;dumx3=gsn_add_polyline(wks,plot(1),(/2.5,2.5/),(/-4,4/),lnres)
   ;dumx4=gsn_add_polyline(wks,plot(1),(/2.5,2.5/)+3,(/-4,4/),lnres)
   ;lnres@gsLineDashPattern=0
   ;lnres@gsLineThicknessF=2.
   ;lnres@gsLineColor="black"

   ;dumx6=gsn_add_polyline(wks,plot1(9),(/2.5,2.5/)+6,(/-4,4/),lnres)


   ;txres@txFontHeightF=0.0085
   ;tx1=gsn_add_text(wks,plot1(9),"Zonal Advection",     .8,   -0.22,txres)
   ;tx2=gsn_add_text(wks,plot1(9),"Meridional Advection",1.+3, -0.22,txres)
   ;tx3=gsn_add_text(wks,plot1(9),"Vertical Advection",  .9+6,  0.06,txres)
   ;tx4=gsn_add_text(wks,plot1(9),"Zonal Advection",     .85+9, -0.22,txres)
   ;tx5=gsn_add_text(wks,plot1(9),"Meridional Advection",1.+12,-0.22,txres)
   ;tx6=gsn_add_text(wks,plot1(9),"Vertical Advection",  1.05+15,-0.22,txres)


   pres=True
   pres@gsnPanelCenter=False
   pres@gsnPanelYWhiteSpacePercent=5
   pres@gsnPanelXWhiteSpacePercent=5

   ;pres@gsnPanelLeft=0.2
   pres@gsnPanelLabelBar=True
   pres@pmLabelBarWidthF=0.6
   pres@pmLabelBarHeightF=0.06
   pres@lbLabelFontHeightF=0.012
   pres@lbLabelBarOn=True;False
   pres@lbBoxEndCapStyle="TriangleBothEnds"
   pres@lbBoxLinesOn=True
   ;lpres@gsnPanelFigureStrings =(/"a-spring","b-spring","c-summer","d-summer","e-autumn","f-autumn","g-winter","h-winter"/);,"c","d"
   ;pres@gsnPanelFigureStrings =(/"a","b","c"/)
   pres@amJust="TopLeft"
   ;maximize_output(wks,pres)
   gsn_panel(wks,plot,(/2,1/),pres)


end
