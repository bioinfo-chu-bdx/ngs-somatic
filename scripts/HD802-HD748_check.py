#!/usr/bin/env python
import sys
import csv
import xlrd
import xlwt
from xlutils.copy import copy
from xlutils.styles import Styles
import subprocess

def representsInt(s): # pour eviter avertissement "nombre ecrit en texte" sous excel
	try: 
		s = int(s)
		return s
	except ValueError:
		return s

################################################################################

xls_path = "/media/n06lbth/sauvegardes_pgm/SBT/Suivi_temoin_HD802-HD748.xls"
HD_tsv_path = sys.argv[1]
sample = sys.argv[2]
run_name = sys.argv[3]

################################################################################

# xls styles
headerStyle = xlwt.easyxf("alignment: horiz centre, wrap on; font: name Calibri, bold on, height 220; borders: left thin, top thin, bottom thin, right thin")
generalStyle = xlwt.easyxf("alignment: horiz left; font: name Calibri, height 220; borders: left thin, top thin, bottom thin, right thin")

# i/o
va_file = open(HD_tsv_path,'r')	
va_reader = csv.reader(va_file,delimiter='\t')

xls_reader = xlrd.open_workbook(xls_path,formatting_info=True)
xls_sheet = xls_reader.sheet_by_index(0)
s = Styles(xls_reader)

final_xls = copy(xls_reader)
final_sheet = final_xls.get_sheet(0)

# header 1
header1 = xls_sheet.row(0)
final_sheet.row(0).write(len(header1),sample+'_'+run_name,headerStyle)

# header 2
header2 = xls_sheet.row(1)
final_sheet.row(1).write(len(header2),'Var.freq',generalStyle)

# variants lines
for i in range(xls_sheet.nrows-2):
	variant_line = xls_sheet.row(i+2)
	variant2check = (xls_sheet.cell(i+2,1).value,xls_sheet.cell(i+2,5).value) # nm c.

	va_file.seek(0)
	va_reader.next()
	for va_line in va_reader:	
		va_variant = (va_line[2],va_line[5]) # nm c.
		va_variant_freq = '?'
		if variant2check == va_variant:
			va_variant_freq = va_line[7]
			break
	final_sheet.row(i+2).write(len(variant_line),representsInt(va_variant_freq),generalStyle)

final_xls.save('suivi_hd_temp.xls') # overwrite previous xls

cmd = subprocess.call(['cp','-f','suivi_hd_temp.xls',xls_path])
cmd = subprocess.call(['rm','-f','suivi_hd_temp.xls'])
