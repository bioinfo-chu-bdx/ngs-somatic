#!/usr/bin/env python
import subprocess
import time

#list = open('/DATA/work/variantBase/temp/run_list_SBT.fullpath.txt','r')

runlist = [
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-130-Run78-LAMv8_2018-70pM-7ech_292',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-134-Run80-LAMv8_2018-adj-70pM-6ech-MM_297',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-139-Run81-2018-04-16-LAMv4-TP53-CM_303',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-141-Run82-LAMv8adj2-TP53-MM_305',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-145-Run84-LAMv8_MM_310',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-149-Run85-LAMv8-TP53-MM-CM_313',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-155-Run87_LAMv8-TP53-MM_321',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-158-RUN87bis-LAMv8-TP53-MM_325',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-159-RUN89-LAMv8-TP53-MM_326',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-163-Run87-ter-LAMv8-TP53-MM_330',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-165-Run90-LAMv8-TP53-MM_332',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-167-Run91-LAMv8-TP53-MM_334',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-170-Run92-LAMv8-TP53-MM_337',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-171-run93-LAMv8-TP53-ABL-MM_338',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-177-run94-2018-07-09-LAMv8-tp53-cm_346',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-187-Run95-LAMv8-TP53-ABL10pM-MM_355',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-189-Run96-LAMv8-MM_358',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-194-Run97-LAMv8-TP53-ABL10pM-FLT310pM-MM_366',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-196-Run98_LAMv8-TP53-ABL10pM-MM_370',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-202-run99-LAMv8-tp53-ABL10pM-CM_377',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-204-RUN99bis-LAMV8-TP53-ABL10Pm-CM-2_380',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-206-Run100-LAMv8-TP53-ABL-FLT3-MM_382',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-208-Run101-LAMv8-TP53-ABL10pM-MM_384',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-210-Run101bis-LAMv8-TP53-ABL10pM-MM_386',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-212-Run103-LAMv8-TP53-FLT3-MM_387',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-215-run104-lamv8-tp53-CM_391',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-217-Run101-nouvelle-chimie-IonC-LAMv8-TP53-ABL-MM_393',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-222-Run106_-LAMv8-TP53-ABL-MM_397',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-225-Run107_LAMv8-TP53-ABL-MM_401',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-228-Run108-LAMv8-TP53-ABL-MM_404',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-230-Run109_LAMv8_ABL_FLT3_MM_406',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-232-Run110_LAMv8-TP53-ABL-FLT3-MM_437',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-234-run111-lamv8-tp53-abl-10-12-18-CM_439',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-236-Run112-LAMv8-ABL-17-12-18-CM_441',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-238-Run113_LAMv8_TP53_ABL_MM_443',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-240-run114-lamv8-tp53-abl-cm_445',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-242-Run115-LAMv8-TP53-ABL-MM_446',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-244-Run116-LAMv8-TP53-ABL-CM_448',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-247-run117-lamv8-tp53-abl-cm_451',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-249-run118-lamv8-tp53-abl-flt3-28-1-19-cm_454',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-253-run119-lamv8-tp53-abl-cm_459',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-255-Run120-LAMv8-TP53-ABL-MM_460',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-257-run121-lamv8-tp53-abl-mm_462',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-259-run123-lamv8-tp53-abl-CF_465',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-260-Run122-LAMv8-ABL-MM-CM_468',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-263-Run124-LAMv8-TP53-ABL-MM_470',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-264-Run125-LAMv8-ABL-MM_471',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-267-Run126-LAMv8-TP53-ABL-FLT3-calS5-MM_474',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-269-Run127_LAMv8-TP53-MM_475',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-270-Run128-LAMv8-ABL-validation-MM_478',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-273-Run129-LAMv8-TP53-MM_479',
'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-274-RUN130--LAMV8-ABL-MM_481'
]

for run in runlist:
	#subprocess.call(['python','/DATA/work/scripts/run_analysis_only_cnv.py','--full-run',run.replace('\n','')])
	#subprocess.call(['python','/DATA/work/scripts/finalReport_replace_cnvsheet.py','--run-folder',run.replace('\n','')])
	
	#runpath = run.replace('\n','')
	#runname = runpath.split('/')[-1]
	
	#subprocess.call(['mkdir',''/media/stuff/rerun_cnv_avec_temoin_EFS/'+runname])
	#subprocess.call(['cp','-r',runpath+'/_CNA',''/media/stuff/rerun_cnv_avec_temoin_EFS/'+runname+'/'])

	subprocess.call(['python','/DATA/work/scripts/mpileup_checkMut.py','--run-folder',run,'--gene','GATA2','--cpos','c.599dup','--run-type','LAM'])
