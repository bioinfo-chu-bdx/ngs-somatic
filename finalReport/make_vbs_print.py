#!/usr/bin/env python
import os
import sys
import json
import sqlite3

def dict_factory(cursor, row):
	d = {}
	for idx, col in enumerate(cursor.description):
		d[col[0]] = row[idx]
	return d

results_folder = sys.argv[1]
barcodes_json = sys.argv[2]

with open(barcodes_json, 'r') as g:
	barcodes_json = json.load(g)

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))

db_path = global_param['VariantBase']
db_con = sqlite3.connect(db_path)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()

panel2project = {}
db_cur.execute("SELECT * FROM Panel")
db_panels = db_cur.fetchall()
for db_panel in db_panels:
	panel2project[db_panel['panelID']] = db_panel['panelProject']
for barcode in barcodes_json:
	barcodes_json[barcode]['project'] = panel2project[barcodes_json[barcode]['panel']]

libs = list(set([barcodes_json[barcode]['library'] for barcode in barcodes_json]))
panels = list(set([barcodes_json[barcode]['panel'] for barcode in barcodes_json]))
projects = list(set([barcodes_json[barcode]['project'] for barcode in barcodes_json]))

for lib in libs:
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

	checkContaReport_to_print = False
	SBT_to_print = False
	HEM_to_print = False
	TP53_to_print = False
	OTHER_to_print = False
	for barcode in barcodes_json.keys():
		if barcodes_json[barcode]['library'] != lib:
			continue
		elif ("HD802-HD748" in barcodes_json[barcode]['sample'].upper()) or ("ACROMETRIX" in barcodes_json[barcode]['sample'].upper()) or ("HORIZON" in barcodes_json[barcode]['sample'].upper()):
			continue
		elif ("-CHECKCONTAMINATION" in barcodes_json[barcode]['sample'].upper()):
			# THIS IS INTERMEDIATE ENTRY FOR H2O ANALYSIS IN BC JSON, to ignore
			continue
		elif ("H2O" in barcodes_json[barcode]['sample'].upper()) or ("H20" in barcodes_json[barcode]['sample'].upper()) or ("NTC" in barcodes_json[barcode]['sample'].upper()):
			checkContaReport_to_print = True
			vbs.write('CheckContaReportlist.Add CurrentDirectory & "\\_checkContamination\\checkContamination_%s.xlsx"\r\n' % barcodes_json[barcode]['sample'])
		elif barcodes_json[barcode]['project'] == "SBT":
			SBT_to_print = True
			vbs.write('SBTFinalReportlist.Add CurrentDirectory & "\%s" & "\\%s_%s.finalReport.xlsx"\r\n' % (barcodes_json[barcode]['sample'],barcodes_json[barcode]['sample'],barcode))
		elif barcodes_json[barcode]['project'] in ["LAM","Leuc","FLT3","ABL1"]:
			HEM_to_print = True
			vbs.write('HEMATOFinalReportlist.Add CurrentDirectory & "\%s" & "\\%s_%s.finalReport.xlsx"\r\n' % (barcodes_json[barcode]['sample'],barcodes_json[barcode]['sample'],barcode))
		elif barcodes_json[barcode]['project'] == "TP53":
			TP53_to_print = True
			vbs.write('TP53FinalReportlist.Add CurrentDirectory & "\%s" & "\\%s_%s.finalReport.xlsx"\r\n' % (barcodes_json[barcode]['sample'],barcodes_json[barcode]['sample'],barcode))
		else:
			OTHER_to_print = True
			vbs.write('OTHERFinalReportlist.Add CurrentDirectory & "\%s" & "\\%s_%s.finalReport.xlsx"\r\n' % (barcodes_json[barcode]['sample'],barcodes_json[barcode]['sample'],barcode))

	if SBT_to_print:
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

	if HEM_to_print:
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

	if TP53_to_print:
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

	if OTHER_to_print: 
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

	### CHECK-CONTA REPORTS
	if checkContaReport_to_print:
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
		# IF "IT IS AN HEMATO RUN". WARNING : if CO-RUN SBT AND HEMATO THIS WILL NOT WORK
		if HEM_to_print:
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
