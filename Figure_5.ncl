load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin


opendir="/data0/fjiang/data/ocean/oras5/postprocess/"
resultdir="/data0/fjiang/scripts/contrast/"

fname1 = resultdir + "figuredata/current-trend-surface.nc";;eq.
fout1 = addfile(fname1, "r")
data1=fout1->data1
data2=fout1->data2

fname2 = resultdir + "figuredata/current-IPO.nc";;eq.
fout2 = addfile(fname2, "r")
data3=fout2->data1

aa=150
bb=270


data1:=runave_n_Wrap(dim_avg_n_Wrap(data1(:,:,{aa:bb}),2)  ,5,0,1)
data2:=runave_n_Wrap(dim_avg_n_Wrap(data2(:,:,{aa:bb}),2)  ,5,0,1)
data3:=runave_n_Wrap(dim_avg_n_Wrap(data3(:,:,{aa:bb}),2)  ,5,0,1)
;data4:=runave_n_Wrap(dim_avg_n_Wrap(data4(:,:,{aa:bb}),2)/2,5,0,1)

data3(0,:)=data3(0,:);/1.5
;data4(0,:)=data4(0,:)/1.5


  wks = gsn_open_wks("eps",resultdir+"v6-5");1979-2018")
   gsn_define_colormap(wks,"MPL_RdBu")
   ;gsn_define_colormap(wks,"MPL_RdBu");"CBR_drywet");
  gsn_reverse_colormap(wks);

   resh = True
   resh@gsnDraw   = False
   resh@gsnFrame  = False
   resh@vpWidthF=0.2
   resh@vpHeightF=0.36
   resh@gsnXRefLine=0.
   resh@gsnYRefLine=21.
   resh@gsnXRefLineColor="gray"
   resh@gsnYRefLineColor="gray"
    resh@gsnXRefLineDashPattern=2.
    resh@gsnYRefLineDashPattern=2.

   ;resh@trYMaxF=.9
   ;resh@trYMinF=-.9
   ;resh@trXMaxF=2.
   ;resh@trXMinF=0.

   resh@tmXTMinorOn=False
   resh@tmXBOn=False
   resh@tmYROn=False
   resh@tmXTBorderOn=True

   ;resh@tmXBBorderOn=False
   ;resh@tmYRBorderOn=False

   resh@tmXTOn=True

   resh@tmXTMode              = "Explicit"         ; explicit labels
   resh@tmYLLabelFontHeightF=0.008
   resh@tmXTLabelFontHeightF=0.008
   resh@tmXBLabelFontHeightF=0.008
   resh@gsnLeftStringFontHeightF=0.012
   resh@gsnLeftStringOrthogonalPosF=0.03


   resh@tmXTMinorOn=False
   resh@tmYLMinorOn=False

   resh@tmYLOn  =True
   resh@tmYLMode              = "Explicit"         ; explicit labels
   resh@tmXBOn  =True
   resh@tmXBMode              = "Explicit"         ; explicit l
   ;resh@tmXTValues   = ispan(0,4,1)   ;
   ;resh@tmXBValues:= 0.5*ispan(0,4,1)   ;
   ;resh@tmXBLabels:= (/"0.0","0.5","1.0","1.5","2.0"/)
   resh@trYMaxF=41-5
   resh@trYMinF=1+5

   resh@tmYLValues:= (/0,5,10,15,20,25,30,35,40/)+1;ispan(0,42,1)
   resh@tmXTLabelsOn=True ;
   resh@tmXBLabelsOn=True ;

   resh@tmYLLabels:= (/"20~S~o~N~S","15~S~o~N~S","10~S~o~N~S","5~S~o~N~S","0~S~o~N~","5~S~o~N~N","10~S~o~N~N","15~S~o~N~N","20~S~o~N~N"/);(/"-1.2","-0.8","-0.4"," 0.0"," 0.4"," 0.8"," 1.2"/)
   ;resh@tmXTLabels   = (/"0.0","1.0","2.0","3.0","4.0"/)

   resh@tmXTMajorLengthF    =0.002
   resh@tmYLMajorLengthF    =0.004
   resh@tmXTMinorLengthF    =0.0;05
   resh@tmXTMajorOutwardLengthF=0.002
   resh@tmYLMajorOutwardLengthF=0.004
   resh@tmXBMajorLengthF    =0.002
   resh@tmXBMajorOutwardLengthF=0.002

   resh@tiMainString    =""; "(a) Zonal current Trend  (b) Meridional current Trend  (c) Vertical current Trend";"CO2: N="+nx+"  nBands="+resh@gsnHistogramNumberOfBins
   ;plot2 = gsn_histogram(wks,x,resh)         ; create histogram with 20 bins
   resh@tiYAxisString=""
   resh@tiXAxisSide="Bottom"
   resh@xyLineColor="deepskyblue3"
   resh@xyLineThicknessF=3.
   plot=new(7,graphic)
   resh@trXMaxF=0.04
   resh@trXMinF=-0.04
   resh@tmXUseBottom=False
   resh@tmXTValues   = -0.02*ispan(-4,4,1)   ;
   resh@tmXTLabels   = -0.02*ispan(-4,4,1)   ;(/"0.0","0.5","1.0","1.5","2.0"/)
   resh@tmXBValues   = -0.02*ispan(-4,4,1)    ;
   resh@tmXBLabels   = -0.02*ispan(-4,4,1)     ;(/"0.0","0.5","1.0","1.5","2.0"/)
   resh@tiXAxisString="[m/s/decade]"
   resh@tiXAxisFontColor="black"
   resh@tiXAxisFontHeightF=0.01
   resh@gsnLeftString="(a) Zonal Current Trend"
   resh@gsnStringFont=22
   resh@vpXF=0.1
   resh@vpYF=0.95
   plot(0) = gsn_csm_xy (wks,data1(0,:),ispan(0,42,1),resh) ; create plot
   resh@trXMaxF=0.006
   resh@trXMinF=-0.006
   resh@tmXTValues   = -0.003*ispan(-4,4,1)  ;   ;
   resh@tmXTLabels   = -0.003*ispan(-4,4,1)  ;(/"0.0","0.5","1.0","1.5","2.0"/)
   resh@tmXBValues   = -0.003*ispan(-4,4,1)  ;   ;
   resh@tmXBLabels   = -0.003*ispan(-4,4,1)  ;(/"0.0","0.5","1.0","1.5","2.0"/)
   resh@gsnLeftString="(b) Meridional Current Trend"
   resh@vpXF=0.375

   plot(1) = gsn_csm_xy (wks,data1(1,:),ispan(0,42,1),resh) ; create plot
   resh@trXMaxF=0.06
   resh@trXMinF=-0.06
   resh@tmXTValues   = -0.03*ispan(-4,4,1)   ;
   resh@tmXTLabels   = -0.03*ispan(-4,4,1)   ;(/"0.0","0.5","1.0","1.5","2.0"/)
   resh@tmXBValues   = -0.03*ispan(-4,4,1)    ;
   resh@tmXBLabels   = -0.03*ispan(-4,4,1)     ;(/"0.0","0.5","1.0","1.5","2.0"/)
   resh@gsnLeftString="(c) Vertical current Trend"
   resh@vpXF=0.65
   resh@tiXAxisString="[m/day/decade]"

   plot(2) = gsn_csm_xy (wks,data1(2,:),ispan(0,42,1),resh) ; create plot
   resh@gsnLeftString=""

   resh@xyLineColor="orangered3"


   plot_1=new(3,graphic)
   plot_1(0) = gsn_csm_xy (wks,data2(0,:),ispan(0,42,1),resh) ; create plot
   plot_1(1) = gsn_csm_xy (wks,data2(1,:),ispan(0,42,1),resh) ; create plot
   plot_1(2) = gsn_csm_xy (wks,data2(2,:),ispan(0,42,1),resh) ; create plot

   resh@xyLineColor="gray40"
   resh@xyLineThicknessF=1.66

   plot_2=new(3,graphic)
   plot_2(0) = gsn_csm_xy (wks,data3(0,:),ispan(0,42,1),resh) ; create plot
   plot_2(1) = gsn_csm_xy (wks,data3(1,:),ispan(0,42,1),resh) ; create plot
   plot_2(2) = gsn_csm_xy (wks,data3(2,:),ispan(0,42,1),resh) ; create plot
   resh@xyDashPattern=16

   ;plot3=new(3,graphic)
   ;plot3(0) = gsn_csm_xy (wks,data4(0,:),ispan(0,42,1),resh) ; create plot
   ;plot3(1) = gsn_csm_xy (wks,data4(1,:),ispan(0,42,1),resh) ; create plot
   ;plot3(2) = gsn_csm_xy (wks,data4(2,:),ispan(0,42,1),resh) ; create plot



overlay(plot(0),plot_1(0))
overlay(plot(1),plot_1(1))
overlay(plot(2),plot_1(2))


overlay(plot(0),plot_2(0))
overlay(plot(1),plot_2(1))
overlay(plot(2),plot_2(2))
;
;overlay(plot(0),plot3(0))
;overlay(plot(1),plot3(1))
;overlay(plot(2),plot3(2))



lgres=True
amres=True
lgres@lgLineColors=(/"gray40","deepskyblue3","orangered3"/)
lgres@lgLeftMarginF=0.1
lgres@lgRightMarginF=0.1

lgres@lgDashIndexes=(/0,0,0/)
lgres@lgLineThicknesses=(/2,4,4/)
lgres@lgMonoLineThickness=False
lgres@vpWidthF=0.085
lgres@vpHeightF=0.05
lgres@lgLabelFontHeightF=0.13
lgres@lgLineDashSegLenF=0.04
lgres@lgPerimThicknessF=2.
lgres@lgPerimOn=False
lgres@lgPerimFill="SolidFill"
lgres@lgPerimFillColor="white"
;lgres@lgLabelOffsetF=-0.1
lgres@lbBoxMajorExtentF=0.3
lbid1=gsn_create_legend(wks,3,(/"   IPO-related","  1958-2022","  1980-2022"/),lgres)
amres@amParallelPosF=-0.26
amres@amOrthogonalPosF=-0.41
annoid1=gsn_add_annotation(plot(0),lbid1,amres)

draw(plot(0))
draw(plot(1))
draw(plot(2))

;;;ipo_obs
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

hor= where(abs(ha).ge.50,ha@_FillValue,ha)
hor= where(ismissing(hor),0,hor)
hor=rmMonAnnCycTLL(hor)
copy_VarMeta(ha, hor)
hor:=dtrend_n(hor,False,0)
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
ipo_obs:=-ezfftb_n(cft1,cft1@xbar,1);
printVarSummary(ipo_obs)
ipo_obs:=dim_standardize_n_Wrap(ipo_obs,1,0)
rc_ipo=120*regCoef_n(ispan(0,43*12-1,1),ipo_obs(110*12::),0,0)
print(rc_ipo)
ipo_obs:=ipo_obs(88*12::)

opendir="/data0/fjiang/data/ocean/oras5/postprocess/"
resultdir="/data0/fjiang/scripts/contrast/"
a = addfile(opendir+"temp.remap.nc", "r")
b = addfile(opendir+"uo.remap.nc", "r")
c = addfile(opendir+"vo.remap.nc", "r")
d = addfile(opendir+"w.remap.cal.nc", "r")

syear1=1980
syear2=1958

eyear    = 2022
ymdStrt1 = syear1*100+1
ymdStrt2 = syear2*100+1

ymdLast1 = eyear*100+12
mm1=(eyear-syear1+1)*12
mm2=(eyear-syear2+1)*12


time1      = b->time_counter
ymd1      = cd_calendar(time1, -1)
iStrt1    = ind(ymd1.eq.ymdStrt1)
iStrt2    = ind(ymd1.eq.ymdStrt2)
iLast1    = ind(ymd1.eq.ymdLast1)

delete(time1)
delete(ymd1)

uo        =b->vozocrtx(:,0,{-21:21},{120:290});dim_avg_n_Wrap(b->vozocrtx(iStrt1:iLast1,0:17,{-21:21},{120:290}),1)
uoannual1  =month_to_annual(uo(iStrt1:iLast1,:,:),1)
uoannual2  =month_to_annual(uo(iStrt2:iLast1,:,:),1)

year1    = ispan(syear1,eyear,1)
year2    = ispan(syear2,eyear,1)

dims=dimsizes(uoannual1)
rc1    = regCoef_n(ipo_obs,dtrend_n(uo(iStrt2:iLast1,:,:),False,0),0,0);year1,uoannual1,0,0)
rc2    = regCoef_n(year2,uoannual2,0,0)


lat=uo&lat
lon=uo&lon



  confi1=student_t(rc1@tval,rc1@nptxy-2)
  confi1:=onedtond(confi1,(/dims(1),dims(2)/))
  rc1@long_name = "Trend"
  rc1:= rc1*rc_ipo                     ; (C/year)*(nyrs)
  copy_VarMeta(uo(0,:,:),rc1)
  copy_VarMeta(rc1,confi1)



    confi2=student_t(rc2@tval,rc2@nptxy-2)
    confi2:=onedtond(confi2,(/dims(1),dims(2)/))
    rc2@long_name = "Trend"
    rc2:= rc2*10                      ; (C/year)*(nyrs)
    copy_VarMeta(uo(0,:,:),rc2)
    copy_VarMeta(rc2,confi2)


  ;;sshannual

  opendir1="/data0/fjiang/data/ocean/oras5/postprocess/"
  a = addfile(opendir1+"ssh.remap.nc", "r")

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

  oras5_d20       =a->sossheig(:,:,:)

  oras5_d20_annual1 =month_to_annual(oras5_d20(iStrt1:iLast1,:,:),1)
  oras5_d20_annual2 =month_to_annual(oras5_d20(iStrt2:iLast1,:,:),1)

  year1    = ispan(syear1,eyear,1)
  year2    = ispan(syear2,eyear,1)

  dims=dimsizes(oras5_d20_annual1)
  rc1_ssh    = regCoef_n(ipo_obs,dtrend_n(oras5_d20(iStrt2:iLast1,:,:),False,0),0,0);year1,oras5_d20_annual1,0,0)
  rc2_ssh    = regCoef_n(year2,oras5_d20_annual2,0,0)
  ;confi1=student_t(rc1@tval,rc1@nptxy-2)
  ;confi1:=onedtond(confi1,(/dims(1),dims(2)/))
  ;rc1@long_name = "Trend"
  rc1_ssh:= rc1_ssh*rc_ipo;*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
  copy_VarMeta(oras5_d20(0,:,:),rc1_ssh)
  ;copy_VarMeta(rc1,confi1)

  ;confi2=student_t(rc2@tval,rc2@nptxy-2)
  ;confi2:=onedtond(confi2,(/dims(1),dims(2)/))
  ;rc2@long_name = "Trend"
  rc2_ssh:= rc2_ssh*10;(eyear-syear2+1);/10                      ; (C/year)*(nyrs)
  copy_VarMeta(oras5_d20(0,:,:),rc2_ssh)
  ;copy_VarMeta(rc2,confi2)



    uu1=rc1_ssh
    uu2=rc2_ssh
    ;vv1=rc1
    ;vv2=rc2
   delete(a)
    a= 6378388.                  ; m
    g=9.8                        ; m/s2
    pi=3.1415926
    omega=2*pi/86400
    phi=doubletofloat(rc1_ssh&lat)
    f=2*omega*sin(phi*pi/180)
    f=where(f.eq.0,f@_FillValue,f)
    pi=3.1415926
    h1_y=center_finite_diff_n(rc1_ssh, pi/180.0 ,False,0,0);(numerator(j,i)-numerator(j-1,i))/(phi(j)-phi(j-1))
    h2_y=center_finite_diff_n(rc2_ssh, pi/180.0 ,False,0,0);(numerator(j,i)-numerator(j-1,i))/(phi(j)-phi(j-1))

    h1_x=center_finite_diff_n(rc1_ssh, pi/180.0 ,False,0,1);(numerator(j,i)-numerator(j-1,i))/(phi(j)-phi(j-1))
    h2_x=center_finite_diff_n(rc2_ssh, pi/180.0 ,False,0,1);(numerator(j,i)-numerator(j-1,i))/(phi(j)-phi(j-1))


    h1_yy=center_finite_diff_n(h1_y, pi/180.0 ,False,0,0);(numerator(j,i)-numerator(j-1,i))/(phi(j)-phi(j-1))
    h2_yy=center_finite_diff_n(h2_y, pi/180.0 ,False,0,0);(numerator(j,i)-numerator(j-1,i))/(phi(j)-phi(j-1))

    h1_xx=center_finite_diff_n(h1_x, pi/180.0 ,False,0,1);(numerator(j,i)-numerator(j-1,i))/(phi(j)-phi(j-1))
    h2_xx=center_finite_diff_n(h2_x, pi/180.0 ,False,0,1);(numerator(j,i)-numerator(j-1,i))/(phi(j)-phi(j-1))


    dims_Tau=dimsizes(rc1_ssh)
    do j=0,dims_Tau(0)-1
    do i=0,dims_Tau(1)-1
     uu1(j,i)=-g/(a*f(j))*h1_y(j,i)
     uu2(j,i)=-g/(a*f(j))*h2_y(j,i)
     ;vv1(j,i)= g/(f(j)*a*cos(phi(j)*pi/180.0))*h1_x(j,i)
    ;vv2(j,i)= g/(f(j)*a*cos(phi(j)*pi/180.0))*h2_x(j,i)
    end do
    end do
    ;uu1=-g/(a*f)*h1_y
    ;uu2=-g/(a*f)*h2_y
    beta=2.3*1e-11;l;(/(m*s);2*omega*cos(lat(nl) * (pi/180.0))

    do i=0,dims_Tau(1)-1
     uu1(dims_Tau(0)/2,i)=-g/(a*a*beta)*h1_yy(dims_Tau(0)/2,i)
     uu2(dims_Tau(0)/2,i)=-g/(a*a*beta)*h2_yy(dims_Tau(0)/2,i)
     ;vv1(dims_Tau(0)/2,i)=g/(a*a*beta)*h1_xx(dims_Tau(0)/2,i)
     ;vv2(dims_Tau(0)/2,i)=g/(a*a*beta)*h2_xx(dims_Tau(0)/2,i)
    end do




    copy_VarMeta(rc1_ssh,uu1)
    copy_VarMeta(rc1_ssh,uu2)
    ;copy_VarMeta(rc1,vv1)
    ;copy_VarMeta(rc1,vv2)

    ;print(vv1({-5:5},{120:290}))

    guess     = 1                ; use zonal means
    is_cyclic = True             ; cyclic [global]
    nscan     = 1500             ; usually much less than this
    eps       = 1.e-2            ; variable dependent
    relc      = 0.6              ; relaxation coefficient
    opt       = 0                ; not used

    poisson_grid_fill(rc1, is_cyclic, guess, nscan, eps, relc, opt)
    poisson_grid_fill(rc2, is_cyclic, guess, nscan, eps, relc, opt)
    poisson_grid_fill(uu1, is_cyclic, guess, nscan, eps, relc, opt)
    poisson_grid_fill(uu2, is_cyclic, guess, nscan, eps, relc, opt)



    delete(a)
    delete(b)
    delete(lat)
    delete(lon)

    opendir="/data0/fjiang/data/ocean/oras5/postprocess/"
    resultdir="/data0/fjiang/scripts/contrast/"
    a = addfile(opendir+"zs.remap.nc", "r")
    b = addfile(opendir+"ms.remap.nc", "r")
    lat=a->lat
    lon=a->lon

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


    ;;
    ws        =a->sozotaux(:,:,:)
    wsannual1  =month_to_annual(ws(iStrt1:iLast1,:,:),1)
    wsannual2  =month_to_annual(ws(iStrt2:iLast1,:,:),1)

    ms        =b->sometauy(:,:,:)
    msannual1  =month_to_annual(ms(iStrt1:iLast1,:,:),1)
    msannual2  =month_to_annual(ms(iStrt2:iLast1,:,:),1)


    year1    = ispan(syear1,eyear,1)
    year2    = ispan(syear2,eyear,1)

    dims=dimsizes(wsannual1)
    rc1_ws    = regCoef_n(ipo_obs,dtrend_n(ws(iStrt2:iLast1,:,:),False,0),0,0);regCoef_n(year1,wsannual1,0,0)
    rc2_ws    = regCoef_n(year2,wsannual2,0,0)
    rc1_ms    = regCoef_n(ipo_obs,dtrend_n(ms(iStrt2:iLast1,:,:),False,0),0,0);regCoef_n(year1,msannual1,0,0)
    rc2_ms    = regCoef_n(year2,msannual2,0,0)
    ;confi1=student_t(rc1@tval,rc1@nptxy-2)
    ;confi1:=onedtond(confi1,(/dims(1),dims(2)/))
    ;rc1@long_name = "Trend"
    rc1_ws:= rc1_ws*rc_ipo;10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
    rc2_ws:= rc2_ws*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
    rc1_ms:= rc1_ms*rc_ipo;*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
    rc2_ms:= rc2_ms*10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
    copy_VarMeta(wsannual1(0,:,:),rc1_ws)
    copy_VarMeta(wsannual1(0,:,:),rc2_ws)
    copy_VarMeta(msannual1(0,:,:),rc1_ms)
    copy_VarMeta(msannual1(0,:,:),rc2_ms)


    delete(a)
    a= 6378388.                  ; m
    g=9.8                        ; m/s2
    pi=3.1415926
    omega=2*pi/86400
    pi=3.1415926

    dims=dimsizes(ws)
    nl=dims(1)
    rs=1/(2*86400.);;/s
    ;;ws: N/m2
    rho = 1037.0    ; density of seawater (kg/m^3)
    omega = 7.2921e-5
    pi=3.14
    beta=2.26*1e-11;l;(/(m*s);2*omega*cos(lat(nl) * (pi/180.0))

    ekman_u1=rc1_ws
    ekman_v1=rc1_ms
    ekman_u2=rc2_ws
    ekman_v2=rc2_ms

    do i=0,nl-1
    y=lat(i)* (pi/180.0)*6378388.
    ff=doubletofloat(beta*y)
    ekman_u1(i,:)=rs*rc1_ws(i,:)+ff*rc1_ms(i,:)/(rho*(ff^2+rs^2))
    ekman_v1(i,:)=(rs*rc1_ms(i,:)-ff*rc1_ws(i,:)/(rho*(ff^2+rs^2)))
    ekman_u2(i,:)=rs*rc2_ws(i,:)+ff*rc2_ms(i,:)/(rho*(ff^2+rs^2))
    ekman_v2(i,:)=(rs*rc2_ms(i,:)-ff*rc2_ws(i,:)/(rho*(ff^2+rs^2)))
    end do


    ekman_u=ws
    ekman_v=ms


    do i=0,nl-1
    y=lat(i)* (pi/180.0)*6378388.
    ff=doubletofloat(beta*y)
    ekman_u(:,i,:)=rs*ws(:,i,:)+ff*ms(:,i,:)/(rho*(ff^2+rs^2))
    ekman_v(:,i,:)=(rs*ms(:,i,:)-ff*ws(:,i,:)/(rho*(ff^2+rs^2)))
    end do
    copy_VarMeta(ws,ekman_u)
    copy_VarMeta(ws,ekman_v)


    Ekman_pumping = 86400*uv2dv_cfd(ekman_u,ekman_v,lat,lon,2)
    if (any(isnan_ieee(Ekman_pumping))) then
        value = 1.e20
        replace_ieeenan (Ekman_pumping, value, 0)
        Ekman_pumping@_FillValue = value
    end if
    copy_VarMeta(ws,Ekman_pumping)

    Ekman_pumping_annual  =month_to_annual(Ekman_pumping(:,:,:),1)
    dims=dimsizes(wsannual1)

    year1    = ispan(syear1,eyear,1)
    year2    = ispan(syear2,eyear,1)

    rc1_ekman    = regCoef_n(ipo_obs,Ekman_pumping(iStrt2:iLast1,:,:),0,0);regCoef_n(year1,wsannual1,0,0)
    rc2_ekman    = regCoef_n(year2,Ekman_pumping_annual,0,0)
    copy_VarMeta(wsannual1(0,:,:),rc1_ekman)
    copy_VarMeta(wsannual1(0,:,:),rc2_ekman)

    confi_rc1_ekman=student_t(rc1_ekman@tval,rc1_ekman@nptxy-2)
    confi_rc1_ekman:=onedtond(confi_rc1_ekman,(/dims(1),dims(2)/))
    copy_VarMeta(rc1_ekman,confi_rc1_ekman)


    confi_rc2_ekman=student_t(rc2_ekman@tval,rc2_ekman@nptxy-2)
    confi_rc2_ekman:=onedtond(confi_rc2_ekman,(/dims(1),dims(2)/))
    copy_VarMeta(rc1_ekman,confi_rc2_ekman)


    rc1_ekman:= rc1_ekman*rc_ipo;10;(eyear-syear1+1);/10                      ; (C/year)*(nyrs)
    rc2_ekman:= rc2_ekman*10;(eyear-syear1+1);/10

    copy_VarMeta(wsannual1(0,:,:),rc1_ekman)
    copy_VarMeta(wsannual1(0,:,:),rc2_ekman)





    rc1_ekman:=where(rc1_ekman.gt.10e2.or.rc1_ekman.lt.-10e2,0,rc1_ekman)
    rc2_ekman:=where(rc2_ekman.gt.10e2.or.rc2_ekman.lt.-10e2,0,rc2_ekman)


    copy_VarMeta(wsannual1(0,:,:),rc1_ekman)
    copy_VarMeta(wsannual2(0,:,:),rc2_ekman)


    printMinMax(rc1_ekman,1)
    printMinMax(rc2_ekman,1)


    guess     = 1                ; use zonal means
    is_cyclic = True             ; cyclic [global]
    nscan     = 1500             ; usually much less than this
    eps       = 1.e-2            ; variable dependent
    relc      = 0.6              ; relaxation coefficient
    opt       = 0                ; not used


    poisson_grid_fill(rc1_ekman, is_cyclic, guess, nscan, eps, relc, opt)
    poisson_grid_fill(rc2_ekman, is_cyclic, guess, nscan, eps, relc, opt)



     res = True
     res@gsnDraw   = False
     res@gsnFrame  = False
     res@gsnAddCyclic      = False
     res@lbLabelBarOn=True
     ;res@vpWidthF =0.9
     res@mpFillOn           = False
     res@mpMinLonF  = 120 ;lon2d(0,0)
     res@mpCenterLonF=180
     res@mpMaxLonF = 290 ;lon2d(0,ilon-1)
     ;res@lbBoxEndCapStyle="TriangleBothEnds"

     res@lbLabelStride=4
     res@lbBoxLinesOn=True;False
     res@mpMinLatF            = -20       ; zoom in on map
     res@mpMaxLatF            = 20
     res@vpHeightF = 0.15
     res@vpWidthF  = 0.34
     res@mpShapeMode="FreeAspect"
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
     res@mpGeophysicalLineColor="gray60"

     res@pmTickMarkDisplayMode = "Always"
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
     res@cnLineLabelsOn=False
     res@cnLinesOn =False
     res@tmXBLabelFontHeightF  = 0.008
     res@tmYLLabelFontHeightF  = 0.008
     res@cnLevelSelectionMode = "ExplicitLevels"  ; use explicit levels
     ;res@cnFillColors=(/111,100,96,87,77,67,64,54,44,34,30,20/);
     res@cnFillColors=(/10,15,20,25,30,35,48,60,0,0,70,82,95,100,105,110,115,120/)

  ;   res@cnLineLabelsOn =True
  ;   res@cnLineLabelBackgroundColor="-1"
  ;   res@cnLineLabelInterval =0.5
  ;   res@cnLineThicknessF=1.5
  ;   res@cnInfoLabelOn=False
  ;
  ;   res@gsnContourNegLineDashPattern=2
  ;   res@gsnContourPosLineDashPattern=0
  ;   res@gsnContourZeroLineThicknessF=0
  ;   res@lbBoxEndCapStyle="TriangleBothEnds"
  ;   res@lbBoxLinesOn=False
  ;   res@pmLabelBarHeightF=0.1
     res@gsnLeftStringFontHeightF=0.01
     res@gsnLeftStringOrthogonalPosF=0.01;0.005
     res@gsnRightStringFontHeightF=0.01
     res@gsnRightStringOrthogonalPosF=0.005;0.005

     res@gsnRightString="";(~S~o~N~C/65 yrs)"
     res@tiMainOffsetYF=-0.012
     res@tiMainFontHeightF=0.0115
  ;
  ;   res@tmXBMajorLengthF=0.004
  ;   res@tmYLMajorLengthF=0.004
  ;   res@tmXBMajorOutwardLengthF=0.004
  ;   res@tmYLMajorOutwardLengthF=0.004
     ;dum1=new(4,graphic)


     res@gsnRightString=""

     res@gsnLeftString=""
     res@tiMainString="(d) IPO-related Geostrophic Zonal Current Trend"

     res@cnLevels             = 0.005*ispan(-8,8,1 )
     res@vpXF=0.1
     res@vpYF=0.497
     res@pmLabelBarOrthogonalPosF=0.2
    ; res@pmLabelBarParallelPosF=-0.11
     res@lbLabelBarOn=True
     res@pmLabelBarHeightF=0.03
     res@pmLabelBarWidthF=0.3
     plot(3) = gsn_csm_contour_map(wks,rc1,res)
     res@tiMainString="(e) Emerging Geostrophic Zonal Current Trend "
     res@lbLabelBarOn=True
     res@vpXF=0.51
     res@vpYF=0.497
     plot(4) = gsn_csm_contour_map(wks,rc2,res)

     ;res@cnLevels             = 0.005*ispan(-8,8,1 )
     res@cnLinesOn=True
     res@cnFillOn=False
     res@cnLineThicknessF=1.6
     res@cnLineColor="black"
     res@gsnContourNegLineDashPattern=1
     res@gsnContourZeroLineThicknessF=0.
     res@cnLineLabelsOn=True
     res@cnLineLabelBackgroundColor=-1
     res@cnLineLabelPlacementMode="Constant"
     res@cnLevels             = 0.01*ispan(-8,8,1 )

     res@gsnLeftString=""
     res@cnConstFLabelConstantSpacingF=1.
     res@tiMainString="";(c) Geostrophic U Trend 1980-2022  "
     plot1a = gsn_csm_contour(wks,uu1,res)
     res@tiMainString="";(d) Geostrophic U Trend 1958-2022  "

     plot1b = gsn_csm_contour(wks,uu2,res)

     overlay(plot(3),plot1a)
     overlay(plot(4),plot1b)

     res@cnLevels             = .03*ispan(-8,8,1 )


    delete(res)
    res = True
    res@gsnDraw   = False
    res@gsnFrame  = False
    res@gsnAddCyclic      = True
    res@mpShapeMode="FreeAspect"
    ;res@vpWidthF =0.9
    res@mpFillOn           = False
    res@mpMinLonF  = 120 ;lon2d(0,0)
    res@mpCenterLonF=180
    res@mpMaxLonF = 290 ;lon2d(0,ilon-1)
    ;res@lbBoxEndCapStyle="TriangleBothEnds"

    res@vpHeightF = 0.15
    res@vpWidthF  = 0.34
    res@lbLabelStride=4.
     res@lbBoxLinesOn=True
    res@mpMinLatF            = -20       ; zoom in on map
    res@mpMaxLatF            = 20
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
    res@mpGeophysicalLineColor="gray60"
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

    res@gsnLeftStringFontHeightF=0.012
    res@gsnLeftStringOrthogonalPosF=0.02;0.005
    res@gsnRightStringFontHeightF=0.012
    res@gsnRightStringOrthogonalPosF=-0.003;5;0.005
    res@cnLineLabelsOn=False
    res@gsnRightString="";(~S~o~N~C/decade)"
    res@gsnLeftString="";
    res@tiMainOffsetYF=-0.012
    res@tiMainFontHeightF=0.0115
    ;dum1=new(4,graphic)


    res@vpHeightF= 0.15            ; change aspect ratio of plot
    ;res@cnFillColors=(/10,15,20,25,30,35,40,45,50,57,0,0,73,80,85,90,95,100,105,110,115,120/)

    ;res@cnLevels             = 0.003*ispan(-10,10,1 )
    ;res@cnLevels             = 0.004*ispan(-10,10,1 )
    res@cnFillColors=(/10,15,20,25,30,35,40,50,0,0,80,90,95,100,105,110,115,120/)
    res@cnLevels             = 0.01*ispan(-8,8,1 )




    res@tiMainString="(f) IPO-related Ekman Pumping Trend"
    res@vpXF=0.1
    res@vpYF=0.23
    res@pmLabelBarOrthogonalPosF=0.2
   ; res@pmLabelBarParallelPosF=-0.11
    res@lbLabelBarOn=True
    res@pmLabelBarHeightF=0.03
    res@pmLabelBarWidthF=0.3
    plot(5) = gsn_csm_contour_map(wks,rc1_ekman,res)
    res@vpXF=0.51
    res@vpYF=0.23
    res@tiMainString="(g) Emerging Ekman Pumping Trend "

    plot(6) = gsn_csm_contour_map(wks,rc2_ekman,res)



    res1=True
    res1@gsnDraw=False
    res1@gsnFrame=False
    res1@cnFillOn=True
    res1@cnLinesOn=False
    res1@cnMonoFillPattern=False
    res1@gsnLeftString=""
    res1@gsnRightString=""
    res1@cnLevelSelectionMode="ExplicitLevels"
    res1@cnLevels=(/0,0.1/)
    res1@cnFillPatterns=(/-1,17,-1/)
    res1@cnFillColors=(/"black","black","black"/)
    res1@cnLineLabelsOn=False
    res1@cnInfoLabelOn=False
    res1@lbLabelBarOn=False

    plot2=new(6,graphic)

    plot2(0)=gsn_csm_contour(wks,confi1,res1)
    plot2(1)=gsn_csm_contour(wks,confi2,res1)
    plot2(2)=gsn_csm_contour(wks,confi_rc1_ekman,res1)
    plot2(3)=gsn_csm_contour(wks,confi_rc2_ekman,res1)
    ;plot2(4)=gsn_csm_contour(wks,confi4_a,res1)
    ;plot2(5)=gsn_csm_contour(wks,confi4,res1)

    ;plot2(3)=gsn_csm_contour(wks,confi4,res1)

  ;  plot2(1)=gsn_csm_contour(wks,confislp,res1)

    overlay(plot(3),plot2(0))
    overlay(plot(4),plot2(1))
    overlay(plot(5),plot2(2))
    overlay(plot(6),plot2(3))
;  overlay(plot1(2),plot2(2))
  ;  overlay(plot1(3),plot2(3))
  ;  overlay(plot1(4),plot2(4))
  ;  overlay(plot1(5),plot2(5))
  ;  ;overlay(plot1(3),plot2(3))


  txres=True
  txres@txFontHeightF=0.012

  gsn_text_ndc(wks, "[m/s/decade]",0.27,0.277, txres)
  ;gsn_text_ndc(wks, "[10~S~-6~N~ m/s/decade",0.5,0.135, txres)
  gsn_text_ndc(wks, "[m/s/decade]",0.69,0.277, txres)
  gsn_text_ndc(wks, "[m/day]",0.27,0.01, txres)
  gsn_text_ndc(wks, "[m/day/decade]",0.69,0.01, txres)

  draw(plot(3))
  draw(plot(4))
  draw(plot(5))
  draw(plot(6))


  ;pres=True
  ;pres@gsnPanelCenter=False
  ;pres@gsnPanelYWhiteSpacePercent=3
  ;;pres@gsnPanelLabelBar =True
  ;;pres@gsnPanelFigureStrings =(/"a-spring","b-spring","c-summer","d-summer","e-autumn","f-autumn","g-winter","h-winter"/);,"c","d"
  ;;pres@gsnPanelFigureStrings =(/"a","b","c"/)
  ;pres@amJust="TopLeft"
  ;gsn_panel(wks,plot,(/1,3/),pres)
  frame(wks)

   end
