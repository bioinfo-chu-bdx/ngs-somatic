#!/usr/bin/env python
import json
import glob
import subprocess
from datetime import date

maxDays = 90
exportedReports = '/media/n06lbth/sauvegardes_pgm/archivedReports/'
month2num = {'Jan':1,'Feb':2,'Mar':3,'Apr':4,'May':5,'Jun':6,'Jul':7,'Aug':8,'Sep':9,'Oct':10,'Nov':11,'Dec':12}

toKeep= ['/media/n06lbth/sauvegardes_pgm/archivedReports/Auto_user_PGM-39-Run_16_validation_SBT_colun_lung_v4_85_068', 	# SBT validation
	'/media/n06lbth/sauvegardes_pgm/archivedReports/Auto_user_PGM-40-RUN17_VALIDATION_SBT_colun_lung_v4_86_069', 	# SBT validation
	'/media/n06lbth/sauvegardes_pgm/archivedReports/Auto_user_PGM-42-RUN18-VALIDATION-SBT_colun_lung_v4_88_072', 	# SBT validation
	'/media/n06lbth/sauvegardes_pgm/archivedReports/Auto_user_PGM-43-RUN19-VALIDATION-SBT-colun-lung-v4_89_074', 	# SBT validation
	'/media/n06lbth/sauvegardes_pgm/archivedReports/Auto_user_PGM-45-Run20-Routine-SBT-colun-lungv4_92_078', 	# SBT validation
	'/media/n06lbth/sauvegardes_pgm/archivedReports/Auto_user_PGM-46-RUN21-SBT-colon_Lung_v4_318v2_94_079', 		# SBT validation
	'/media/n06lbth/sauvegardes_pgm/archivedReports/Run_validpanel_name_corrected_039', 		# SBT validation (= Auto_user_S5-0198-2-Run_validpanel_PGM151_Chef_SBT_colon_lung_v5_530_90)
	'/media/n06lbth/sauvegardes_pgm/archivedReports/Auto_user_S5-0198-3-Run_validation_PGM159CQ160_40pM_Chef_SBT-colon-lung_v5_530_99_022',
	'/media/n06lbth/sauvegardes_pgm/archivedReports/Auto_user_PGM-49-run17-LAMv4-25-06-2015_97_084', 		# LAM validation
	'/media/n06lbth/sauvegardes_pgm/archivedReports/Auto_user_PGM-50-RUN18-LAMv4-26-06-2015_98_085', 		# LAM validation
	'/media/n06lbth/sauvegardes_pgm/archivedReports/Auto_user_PGM-52-RUN19-LAMv4-24-07-2015_101_095', 		# LAM validation
	'/media/n06lbth/sauvegardes_pgm/archivedReports/Auto_user_PGM-54-run20-LAMv4-09-10-2015_103_107', 		# LAM validation
	'/media/n06lbth/sauvegardes_pgm/archivedReports/Auto_user_PGM-71-RUN21_LAMv4_25-11-2015_129_127', 		# LAM validation
	'/media/n06lbth/sauvegardes_pgm/archivedReports/Auto_TST_DATA-0_5224_9680c70_161_130',				# TEST mis en service, stockage dans le doute
	'/media/n06lbth/sauvegardes_pgm/archivedReports/Auto_user_S5-0198-15-NGSValidPanel_175_35pM_Chef_SBT-colon-lung_v7_530_118_052'
	]

runs_stored = glob.glob(exportedReports + '*')
print "%s runs stored." % len(runs_stored)
for run_folder in runs_stored:
	try:
		if run_folder in toKeep:
			continue
		#print "- checking run : " + run_folder
		s = run_folder.split('/')[-1].split('_')[0:-1]
		json_name = '_'
		json_name = json_name.join(s)
		json_file = open(run_folder+'/serialized_'+json_name+'.json','r')
		json_data = json.load(json_file)
		json_file.close()
		platform = json_data[1]['fields']['platform']
		log = json_data[1]['fields']['log'].split(',')
		for item in log:
			if 'start_time' in item:
				start_time = item.replace('"','').replace('start_time:','')
				break
	#	print "\t --Start Time: " + start_time
		if platform == 'PGM': # start time look like :  Tue Jan 24 08:54:47 2017
			start_time = start_time.split(' ')
			if '' in start_time:
				start_time.remove('')
			day = int(start_time[2])
			month = month2num[start_time[1]]
			year = int(start_time[4])
			run_date = date(year,month,day)
		elif platform == 'S5': # start time look like :  01/18/2017 10:23:51
			start_time = start_time.split(' ')[0].split('/')
			day = int(start_time[1])
			month = int(start_time[0])
			year = int(start_time[2])
			run_date = date(year,month,day)

		# calcul nb de jours entre date explog et date aujourd'hui, si > x jours, delete dossier
		today = date.today()
		delta = today - run_date
#		print '\t --Run age : ' + str(delta.days) + ' days'

		if delta.days >= maxDays:
			print "- RUN : " + run_folder
			print "\t --Launch Date: %s %s %s" % (str(day),start_time[1],year)
			print "\t --Run age : " + str(delta.days) + " days"
			print "\t -> This run can be deleted"
			# uncomment for auto delete folder (not recommended)
	#		cmd = Popen(['rm','-rf',run_folder], stdout=PIPE)
	#		out, err = cmd.communicate()
	#		print "-Removing previous xls : \nOUT: %s\nERR: %s"%(out, err)

	except:
		print "warning : can't check " + run_folder
