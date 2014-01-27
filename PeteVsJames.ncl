load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin

;
; Load interpolated data files and set plot title
;



plot = new(6,graphic)

panel_res = True
panel_res@cnInfoLabelOn = False
panel_res@gsnPaperOrientation = "auto"          
panel_res@gsnMaximize         = True 
;
; Read data
;	

dataDir="/Volumes/Warehouse/IUTAMPaper/gaussVort/contourPlots/"

wks = gsn_open_wks("pdf",dataDir+"gaussVort_contours3b")

; James's data
nLat = 90
nLon = 180
deg2file0=dataDir+"g2x2initRelVort.dat"
deg2file1=dataDir+"g2x2t1relVort.dat"
deg2file2=dataDir+"g2x2t2relVort.dat"

deg2t0 = asciiread(deg2file0,(/ nLat, nLon /), "double")
deg2t1 = asciiread(deg2file1,(/ nLat, nLon /), "double")
deg2t2 = asciiread(deg2file2,(/ nLat, nLon /), "double")


; Pete's data
nLat = 181
nLon = 360

tri5file0=dataDir+"smooth_gvTri5_0__relVort.dat"
tri5file1=dataDir+"smooth_gvTri5_1__relVort.dat"
tri5file2=dataDir+"smooth_gvTri5_2__relVort.dat"

tri5t0 = asciiread(tri5file0,(/ nLat, nLon /), "double")
tri5t1 = asciiread(tri5file1,(/ nLat, nLon /), "double")
tri5t2 = asciiread(tri5file2,(/ nLat, nLon /), "double")




; set contour levels
	mnlvl = -12.0
	mxlvl = 12.0
	lvlSpcng = 0.5
	nCntrs = (mxlvl-mnlvl)/lvlSpcng + 1
	
	;gsn_define_colormap(wks,"BlWhRe")	
	gsn_define_colormap(wks,"BlueYellowRed")	
	;gsn_define_colormap(wks,"BlueWhiteOrangeRed")
	;gsn_define_colormap(wks,"testcmap")
;
;	Map Resources
;
	map_res = True
	map_res@gsnDraw = False
	map_res@gsnFrame = False
	map_res@gsnMaximize = False
	map_res@mpGridAndLimbOn = True
	map_res@mpGridSpacingF = 90.0
	map_res@mpLimbLineColor = 1
	map_res@mpGridLineColor = 0
	map_res@mpGreatCircleLinesOn = True
	map_res@mpPerimLineThicknessF = 1.0
	map_res@mpCenterLonF = -180.0
	map_res@mpOutlineOn = False ; no continents
	
    map_res@tiMainOn = True
    map_res@tiMainFontHeightF = 0.04
	
	map_res@vpXF = 0.1
	map_res@vpYF = 0.92
	map_res@vpWidthF = 0.8
	map_res@vpHeightF = 0.8

;
; Contour line resources
;	
	map_res@cnLevelSelectionMode = "ManualLevels"
	map_res@cnMinLevelValF = mnlvl
	map_res@cnMaxLevelValF = mxlvl
	map_res@cnLevelSpacingF = lvlSpcng
	
	map_res@cnInfoLabelOn = False
	map_res@cnLineLabelsOn = False
	map_res@cnInfoLabelPerimOn = False
	
	map_res@cnFillOn = True
	map_res@gsnSpreadColors = True
	map_res@gsnSpreadColorStart = 2
	map_res@gsnSpreadColorEnd = -1
	
	map_res@lbLabelBarOn = False
	
	map_res@lbLabelStride = 4
	map_res@lbTitleFontThicknessF = 0.08
	map_res@lbAutoManage = False
	map_res@lbLabelFontHeightF = 0.0175
	
	map_res@pmTickMarkDisplayMode = "Always"
	
	map_res@cnLinesOn = False
	map_res@cnLineDashPattern = 0
	map_res@cnLineThicknessF = 1.0
	map_res@cnMonoLineColor = True
	map_res@cnLineColors = 1
	
	map_res@cnLevelFlag = True
	
	map_res@sfXCStartV = -180.0
	map_res@sfXCEndV = 180.0
	map_res@sfYCStartV = -90.0
	map_res@sfYCEndV = 90.0
	
	;
	; Annotation
	;
	amres = True
	amres@amParallelPosF = -0.48
	amres@amOrthogonalPosF = -0.45
	amres@amJust = "TopLeft"
	
	txres = True
	txres@txPerimOn = True
	txres@txBackgroundFillColor = "White"
	txres@txFontHeightF = 0.045
	
	txid = gsn_create_text(wks,"a)",txres)
	txid2 = gsn_create_text(wks,"b)",txres)
	
	;
	; Build the plots
	;	
	map_res@tiMainString = "day 0"
	
	plot(0) = gsn_csm_contour_map_ce(wks,deg2t0,map_res)
	annid = gsn_add_annotation(plot(0),txid,amres)
	map_res@tiMainOn = False
	
	plot(3) = gsn_csm_contour_map_ce(wks,tri5t0,map_res)
	annoid = gsn_add_annotation(plot(3),txid2,amres)
	
;	mnlvl = -5.0
;	mxlvl = 6.0
;	lvlSpcng = 0.5
;	map_res@cnMinLevelValF = mnlvl
;	map_res@cnMaxLevelValF = mxlvl
;	map_res@cnLevelSpacingF = lvlSpcng
	
	map_res@tiMainString = "day 1"
	map_res@tiMainOn = True
	plot(1) = gsn_csm_contour_map_ce(wks,deg2t1,map_res)
	map_res@tiMainOn = False
	plot(4) = gsn_csm_contour_map_ce(wks,tri5t1,map_res)
	
;	mnlvl = -5.0
;	mxlvl = 6.0
;	lvlSpcng = 0.5
;	map_res@cnMinLevelValF = mnlvl
;	map_res@cnMaxLevelValF = mxlvl
;	map_res@cnLevelSpacingF = lvlSpcng
	
	map_res@tiMainOn = True
	map_res@tiMainString = "day 2"
	plot(2) = gsn_csm_contour_map_ce(wks,deg2t2,map_res)
	map_res@tiMainOn = False
	plot(5) = gsn_csm_contour_map_ce(wks,tri5t2,map_res)
	
	
		
	;
	; Draw the paneled plots
	;	
	;panel_res@txString = "Lin-Rood PPM, 2 x 2 deg"
	;panel_res@txFontHeightF = 0.02
	panel_res@gsnPanelLabelBar = True
	panel_res@lbLabelFontHeightF = 0.01
	panel_res@lbLabelStride = 2
	panel_res@gsnPanelBottom   = 0.05
	gsn_panel(wks,plot,(/2,3/),panel_res)

	
	
end