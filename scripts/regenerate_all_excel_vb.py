#!/usr/bin/env python
import os
import sys
import subprocess

#projects = ['SBT','LAM,FLT3','Lymphome_B','Lymphome_T','TP53','Leuc','ABL1']

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
cmd0 = subprocess.Popen(['python','%s/variantBase/find_sample_pathology_SBT.py'])
cmd0.wait()

cmd1 = subprocess.Popen(['python','%s/variantBase/make_excel_VariantBase.py' % pipeline_folder,'--project','SBT','--output','/media/n06lbth/sauvegardes_pgm/SBT/VariantBase_SBT.xlsx'])
cmd2 = subprocess.Popen(['python','%s/variantBase/make_excel_VariantBase.py' % pipeline_folder,'--project','LAM,FLT3','--output','/media/n06lbth/sauvegardes_pgm/LAM/VariantBase_LAM_FLT3.xlsx'])
cmd3 = subprocess.Popen(['python','%s/variantBase/make_excel_VariantBase.py' % pipeline_folder,'--project','Lymphome_B','--output','/media/n06lbth/sauvegardes_pgm/Lymphome_B/VariantBase_Lymphome_B.xlsx'])
cmd4 = subprocess.Popen(['python','%s/variantBase/make_excel_VariantBase.py' % pipeline_folder,'--project','Lymphome_T','--output','/media/n06lbth/sauvegardes_pgm/Lymphome_T/VariantBase_Lymphome_T.xlsx'])
cmd5 = subprocess.Popen(['python','%s/variantBase/make_excel_VariantBase.py' % pipeline_folder,'--project','TP53','--output','/media/n06lbth/sauvegardes_pgm/LAM/VariantBase_TP53.xlsx'])
cmd6 = subprocess.Popen(['python','%s/variantBase/make_excel_VariantBase.py' % pipeline_folder,'--project','Leuc','--output','/media/n06lbth/sauvegardes_pgm/LAM/VariantBase_Leuc.xlsx'])
cmd7 = subprocess.Popen(['python','%s/variantBase/make_excel_VariantBase.py' % pipeline_folder,'--project','ABL1','--output','/media/n06lbth/sauvegardes_pgm/ABL1/VariantBase_ABL1.xlsx'])

cmd1.wait()
cmd2.wait()
cmd3.wait()
cmd4.wait()
cmd5.wait()
cmd6.wait()
cmd7.wait()
