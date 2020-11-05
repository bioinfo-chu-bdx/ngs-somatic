#!/usr/bin/env python
import os
import sys
import json

results_folder = sys.argv[1]
barcodes_json = sys.argv[2]

##########################################################
with open(barcodes_json, 'r') as g:
	bjson = json.load(g)

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))

sample2barcode = {}
barcode2lib = {}
barcode2runtype = {}
checkContaReport = False

for barcode in bjson.keys():
	sample2barcode[bjson[barcode]['sample']] = barcode
	barcode2lib[barcode] = bjson[barcode]['barcode_description']
	bed_name = bjson[barcode]['target_region_filepath'].split('/unmerged/detail/')[-1]
	for rt in global_param['run_type']:
		if global_param['run_type'][rt]['target_bed'].split('/')[-1] == bed_name:
			barcode2runtype[barcode] = rt
			break

runtypes = barcode2runtype.values()

for lib in list(set(barcode2lib.values())):
	vbs = open("%s/PRINT_Annotation_FinalReport_%s.vbs" % (results_folder,lib), 'w')
	vbs.write('Dim iAnswer\r\n')
	vbs.write('iAnswer = MsgBox("Voulez-vous lancer l\'impression?", vbOKCancel + vbQuestion, "Continue")\r\n')
	vbs.write('if iAnswer = vbCancel Then\r\n')
	vbs.write('WScript.quit()\r\n')
	vbs.write('End If\r\n')
	vbs.write('Set exl = WScript.CreateObject("Excel.Application")\r\n')
	vbs.write('exl.Visible = True \'False\r\n')
	vbs.write('Set sh = WScript.CreateObject("WScript.Shell")\r\n')
	vbs.write('on error resume next\r\n')
	vbs.write('sh.RegWrite "HKEY_CURRENT_USER\Software\Microsoft\Office\\11.0\Excel\Security\\accessVBOM",1,"REG_DWORD"\r\n')
	vbs.write('on error goto 0\r\n')
	vbs.write('dim fso: set fso = CreateObject("Scripting.FileSystemObject")\r\n')
	vbs.write('dim CurrentDirectory\r\n')
	vbs.write('CurrentDirectory = fso.GetAbsolutePathName(".")\r\n')
	vbs.write('dim SBTFinalReportlist\r\n')
	vbs.write('Set SBTFinalReportlist = WScript.CreateObject("System.Collections.ArrayList")\r\n')
	vbs.write('dim HEMATOFinalReportlist\r\n')
	vbs.write('Set HEMATOFinalReportlist = WScript.CreateObject("System.Collections.ArrayList")\r\n')
	vbs.write('dim TP53FinalReportlist\r\n')
	vbs.write('Set TP53FinalReportlist = WScript.CreateObject("System.Collections.ArrayList")\r\n')
	vbs.write('dim OTHERFinalReportlist\r\n')
	vbs.write('Set OTHERFinalReportlist = WScript.CreateObject("System.Collections.ArrayList")\r\n')
	vbs.write('dim CheckContaReportlist\r\n')
	vbs.write('Set CheckContaReportlist = WScript.CreateObject("System.Collections.ArrayList")\r\n')
	
	for sample in sample2barcode.keys():
		if barcode2lib[sample2barcode[sample]] != lib:
			continue
		elif ("HD802-HD748" in sample.upper()) or ("ACROMETRIX" in sample.upper()) or ("HORIZON" in sample.upper()):
			continue
		elif ("H2O" in sample.upper()) or ("H20" in sample.upper()) or ("NTC" in sample.upper()):
			checkContaReport = True
			vbs.write('CheckContaReportlist.Add CurrentDirectory & "\\_checkContamination\\checkContamination_%s.xlsx"\r\n' % sample)
		elif barcode2runtype[sample2barcode[sample]] == "SBT":
			vbs.write('SBTFinalReportlist.Add CurrentDirectory & "\%s" & "\\%s_%s.finalReport.xlsx"\r\n' % (sample,sample,sample2barcode[sample]))
		elif barcode2runtype[sample2barcode[sample]] in ["LAM","Leuc","FLT3","ABL1"]:
			vbs.write('HEMATOFinalReportlist.Add CurrentDirectory & "\%s" & "\\%s_%s.finalReport.xlsx"\r\n' % (sample,sample,sample2barcode[sample]))
		elif barcode2runtype[sample2barcode[sample]] == "TP53":
			vbs.write('TP53FinalReportlist.Add CurrentDirectory & "\%s" & "\\%s_%s.finalReport.xlsx"\r\n' % (sample,sample,sample2barcode[sample]))
		else:
			vbs.write('OTHERFinalReportlist.Add CurrentDirectory & "\%s" & "\\%s_%s.finalReport.xlsx"\r\n' % (sample,sample,sample2barcode[sample]))

	if "SBT" in runtypes:
		vbs.write('For Each SBTFinalReport in SBTFinalReportlist\r\n')
		vbs.write('set fichxl = exl.workbooks.add(SBTFinalReport)\r\n')
		vbs.write('Set mdle = fichxl.VBProject.VBComponents.Add(1)\r\n')
		vbs.write('num=0\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Sub PrintPage1()"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "With Sheets(""Annotation"").PageSetup"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".Orientation = xlLandscape"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".CenterHeader = ""Page &A"""\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".RightHeader = ""Printed &D & &T"""\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".CenterFooter = ""&F"""\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "End With"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Application.PrintCommunication = False"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "With Sheets(""Annotation"").PageSetup"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".FitToPagesWide = 1"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".FitToPagesTall = False"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".LeftMargin = Application.InchesToPoints(0.1)"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".RightMargin = Application.InchesToPoints(0.1)"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "End With"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Application.PrintCommunication = True"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "NbLignes = Sheets(""Annotation"").UsedRange.Rows.Count"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Sheets(""Annotation"").Range(""A1:AZ1"").Interior.Color = vbWhite"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Sheets(""Annotation"").Range(""K:K"").EntireColumn.Hidden = True"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Sheets(""Annotation"").Range(""A1:U"" & NbLignes).printOut"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "End Sub"\r\n')
		vbs.write('exl.Application.Run "PrintPage1"\r\n')
		vbs.write('exl.DisplayAlerts = False\r\n')
		vbs.write('exl.ActiveWorkbook.Close False\r\n')
		vbs.write('Next\r\n')
			
	if ("LAM" in runtypes) or ("Leuc" in runtypes) or ("FLT3" in runtypes) or ("ABL1" in runtypes):
		vbs.write('For Each HEMATOFinalReport in HEMATOFinalReportlist\r\n')
		vbs.write('set fichxl = exl.workbooks.add(HEMATOFinalReport)\r\n')
		vbs.write('Set mdle = fichxl.VBProject.VBComponents.Add(1)\r\n')
		vbs.write('num=0\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Sub PrintPage1()"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "With Sheets(""Annotation"").PageSetup"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".Orientation = xlLandscape"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".CenterHeader = ""Page &A"""\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".RightHeader = ""Printed &D & &T"""\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".CenterFooter = ""&F"""\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "End With"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Application.PrintCommunication = False"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "With Sheets(""Annotation"").PageSetup"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".FitToPagesWide = 1"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".FitToPagesTall = False"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".LeftMargin = Application.InchesToPoints(0.1)"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".RightMargin = Application.InchesToPoints(0.1)"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "End With"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Application.PrintCommunication = True"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "NbLignes = Sheets(""Annotation"").UsedRange.Rows.Count"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Sheets(""Annotation"").Range(""A1:AZ1"").Interior.Color = vbWhite"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Sheets(""Annotation"").Range(""F:H,N:N,T:AF,AH:AI"").EntireColumn.Hidden = True"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Sheets(""Annotation"").Range(""A1:AJ"" & NbLignes).printOut"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "End Sub"\r\n')
		vbs.write('exl.Application.Run "PrintPage1"\r\n')
		vbs.write('exl.DisplayAlerts = False\r\n')
		vbs.write('exl.ActiveWorkbook.Close False\r\n')
		vbs.write('Next\r\n')
	
	if "TP53" in runtypes:
		vbs.write('For Each TP53FinalReport in TP53FinalReportlist\r\n')
		vbs.write('set fichxl = exl.workbooks.add(TP53FinalReport)\r\n')
		vbs.write('Set mdle = fichxl.VBProject.VBComponents.Add(1)\r\n')
		vbs.write('num=0\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Sub PrintPage1()"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "With Sheets(""Annotation"").PageSetup"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".Orientation = xlLandscape"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".CenterHeader = ""Page &A"""\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".RightHeader = ""Printed &D & &T"""\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".CenterFooter = ""&F"""\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "End With"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Application.PrintCommunication = False"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "With Sheets(""Annotation"").PageSetup"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".FitToPagesWide = 1"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".FitToPagesTall = False"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".LeftMargin = Application.InchesToPoints(0.1)"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".RightMargin = Application.InchesToPoints(0.1)"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "End With"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Application.PrintCommunication = True"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "NbLignes = Sheets(""Annotation"").UsedRange.Rows.Count"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Sheets(""Annotation"").Range(""A1:AZ1"").Interior.Color = vbWhite"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Sheets(""Annotation"").Range(""F:H,N:N,V:AH,AJ:AK"").EntireColumn.Hidden = True"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Sheets(""Annotation"").Range(""A1:AL"" & NbLignes).printOut"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "End Sub"\r\n')
		vbs.write('exl.Application.Run "PrintPage1"\r\n')
		vbs.write('exl.DisplayAlerts = False\r\n')
		vbs.write('exl.ActiveWorkbook.Close False\r\n')
		vbs.write('Next\r\n')
		
	for runtype in set(runtypes):
		if runtype not in ["SBT","LAM","Leuc","FLT3","ABL1","TP53"]: # == OTHER
			vbs.write('For Each OTHERFinalReport in OTHERFinalReportlist\r\n')
			vbs.write('set fichxl = exl.workbooks.add(OTHERFinalReport)\r\n')
			vbs.write('Set mdle = fichxl.VBProject.VBComponents.Add(1)\r\n')
			vbs.write('num=0\r\n')
			vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Sub PrintPage1()"\r\n')
			vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "With Sheets(""Annotation"").PageSetup"\r\n')
			vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".Orientation = xlLandscape"\r\n')
			vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".CenterHeader = ""Page &A"""\r\n')
			vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".RightHeader = ""Printed &D & &T"""\r\n')
			vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".CenterFooter = ""&F"""\r\n')
			vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "End With"\r\n')
			vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Application.PrintCommunication = False"\r\n')
			vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "With Sheets(""Annotation"").PageSetup"\r\n')
			vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".FitToPagesWide = 1"\r\n')
			vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".FitToPagesTall = False"\r\n')
			vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".LeftMargin = Application.InchesToPoints(0.1)"\r\n')
			vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".RightMargin = Application.InchesToPoints(0.1)"\r\n')
			vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "End With"\r\n')
			vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Application.PrintCommunication = True"\r\n')
			vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "NbLignes = Sheets(""Annotation"").UsedRange.Rows.Count"\r\n')
			vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Sheets(""Annotation"").Range(""A1:AZ1"").Interior.Color = vbWhite"\r\n')
			vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Sheets(""Annotation"").Range(""A1:S"" & NbLignes).printOut"\r\n')
			vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "End Sub"\r\n')
			vbs.write('exl.Application.Run "PrintPage1"\r\n')
			vbs.write('exl.DisplayAlerts = False\r\n')
			vbs.write('exl.ActiveWorkbook.Close False\r\n')
			vbs.write('Next\r\n')
			break
				
	###
	### CHECK-CONTA REPORTS
	###
	if checkContaReport:
		vbs.write('For Each CheckContaReport in CheckContaReportlist\r\n')
		vbs.write('set fichxl = exl.workbooks.add(CheckContaReport)\r\n')
		vbs.write('Set mdle = fichxl.VBProject.VBComponents.Add(1)\r\n')
		vbs.write('num=0\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Sub PrintPage1()"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "With Sheets(""Annotation"").PageSetup"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".Orientation = xlLandscape"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".CenterHeader = ""Page &A"""\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".RightHeader = ""Printed &D & &T"""\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".CenterFooter = ""&F"""\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "End With"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Application.PrintCommunication = False"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "With Sheets(""Annotation"").PageSetup"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".FitToPagesWide = 1"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".FitToPagesTall = False"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".LeftMargin = Application.InchesToPoints(0.1)"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, ".RightMargin = Application.InchesToPoints(0.1)"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "End With"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Application.PrintCommunication = True"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "NbLignes = Sheets(""Annotation"").UsedRange.Rows.Count"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Sheets(""Annotation"").Range(""A1:AZ1"").Interior.Color = vbWhite"\r\n')
		if ("LAM" in runtypes) or ("TP53" in runtypes) or ("Leuc" in runtypes) or ("FLT3" in runtypes) or ("ABL1" in runtypes):
			vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Sheets(""Annotation"").Range(""A1:S"" & NbLignes).printOut"\r\n')
		else:
			vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Sheets(""Annotation"").Range(""J:J"").EntireColumn.Hidden = True"\r\n')
			vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "Sheets(""Annotation"").Range(""A1:U"" & NbLignes).printOut"\r\n')
		vbs.write('num=num+1:mdle.CodeModule.InsertLines num, "End Sub"\r\n')
		vbs.write('exl.Application.Run "PrintPage1"\r\n')
		vbs.write('exl.DisplayAlerts = False\r\n')
		vbs.write('exl.ActiveWorkbook.Close False\r\n')
		vbs.write('Next\r\n')
		
	vbs.close()
