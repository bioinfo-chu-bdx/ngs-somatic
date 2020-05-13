#!/usr/bin/env python
from subprocess import *
import httplib2
import urllib2
import time
import json

Run2PK = {
#	'Auto_S5-0198-0-TF_Run_5224_9680c70_85':16,
#	'Auto_S5-0198-0-TF_Run_5224_9680c70_85_tn':17,
#	'Auto_TST_DATA-0_5224_9680c70_83':10,
	'Auto_user_PGM-255-Run164_35pM_Chef_SBT_colon_lung_v5_318v2_327':120,
	'Auto_user_S5-0198-10-NGS169NEW_35pM_Chef_SBT-colon-lung_v5_530_113':40,
	'Auto_user_S5-0198-10-NGS169NEW_35pM_Chef_SBT-colon-lung_v5_530_113_tn':41,
	'Auto_user_S5-0198-11-NGS170_171_35pM_Chef_SBT-colon-lung_v5_530_112':42,
	'Auto_user_S5-0198-11-NGS170_171_35pM_Chef_SBT-colon-lung_v5_530_112_tn':43,
	'Auto_user_S5-0198-12-NGS_168SDD_172_Chef_SBT-colon-lung_v5_530_114':44,
	'Auto_user_S5-0198-12-NGS_168SDD_172_Chef_SBT-colon-lung_v5_530_114_tn':45,
##	'Auto_user_S5-0198-13-Run54_Chef_LAMadj-ssadj-60-70pM-TP53_530_115':46,
##	'Auto_user_S5-0198-13-Run54_Chef_LAMadj-ssadj-60-70pM-TP53_530_115_tn':47,
	'Auto_user_S5-0198-14-NGS173_174_35pM_Chef_SBT-colon-lung_v5_530_117':48,
	'Auto_user_S5-0198-14-NGS173_174_35pM_Chef_SBT-colon-lung_v5_530_117_tn':49,
	'Auto_user_S5-0198-15-NGSValidPanel_175_35pM_Chef_SBT-colon-lung_v7_530_118':52,
	'Auto_user_S5-0198-15-NGSValidPanel_175_35pM_Chef_SBT-colon-lung_v7_530_118_tn':53,
	'Auto_user_S5-0198-16-NGS176_177_35pM_Chef_SBT-colon-lung_v7_530_119':54,
	'Auto_user_S5-0198-16-NGS176_177_35pM_Chef_SBT-colon-lung_v7_530_119_tn':55,
	'Auto_user_S5-0198-17-NGS178_179_35pM_Chef_SBT-colon-lung_v7_530_121':56,
	'Auto_user_S5-0198-17-NGS178_179_35pM_Chef_SBT-colon-lung_v7_530_121_tn':57,
##	'Auto_user_S5-0198-18-Run55-LAMv4-TP53-530-13-03-2017_120':58,
##	'Auto_user_S5-0198-18-Run55-LAMv4-TP53-530-13-03-2017_120_tn':59,
#	'Auto_user_S5-0198-19-NGS180_181_35pM_Chef_SBT-colon-lung_v7_530_122':60,
#	'Auto_user_S5-0198-19-NGS180_181_35pM_Chef_SBT-colon-lung_v7_530_122_tn':61,
#	'Auto_user_S5-0198-1-ctrl_S5_CEPH_formation_11_01_17_91':18,
#	'Auto_user_S5-0198-1-ctrl_S5_CEPH_formation_11_01_17_91_tn':19,
#	'Auto_user_S5-0198-20-NGStheseClemence_Chef_SBT-colon-lung_v7_530_123':62,
#	'Auto_user_S5-0198-20-NGStheseClemence_Chef_SBT-colon-lung_v7_530_123_tn':63,
	'Auto_user_S5-0198-21-NGS180_181bis_Chef_SBT-colon-lung_v7_530_124':64,
	'Auto_user_S5-0198-21-NGS180_181bis_Chef_SBT-colon-lung_v7_530_124_tn':65,
	'Auto_user_S5-0198-22-NGS_182_183_35pM_Chef_SBT-colon-lung_v7_530_126':66,
	'Auto_user_S5-0198-22-NGS_182_183_35pM_Chef_SBT-colon-lung_v7_530_126_tn':67,
##	'Auto_user_S5-0198-23-Run56-LAMv4-TP53-adjSRSF2_530_125':68,
##	'Auto_user_S5-0198-23-Run56-LAMv4-TP53-adjSRSF2_530_125_tn':69,
##	'Auto_user_S5-0198-24-Run56-LAMv4adjSRSF2-Tp53_530-30-03-2017_127':70,
##	'Auto_user_S5-0198-24-Run56-LAMv4adjSRSF2-Tp53_530-30-03-2017_127_tn':71,
	'Auto_user_S5-0198-25-NGS184_185_35pM_Chef_SBT-colon-lung_v7_530_128':72,
	'Auto_user_S5-0198-25-NGS184_185_35pM_Chef_SBT-colon-lung_v7_530_128_tn':73,
#	'Auto_user_S5-0198-26-NGS_theseClemencebis_35pM_Chef_SBT-colon-lung_v7_530_129':74,
#	'Auto_user_S5-0198-26-NGS_theseClemencebis_35pM_Chef_SBT-colon-lung_v7_530_129_tn':75,
#	'Auto_user_S5-0198-27-NGS_186_187_35pM_Chef_SBT-colon-lung_v7_530_132':76,
#	'Auto_user_S5-0198-27-NGS_186_187_35pM_Chef_SBT-colon-lung_v7_530_132_tn':77,
##	'Auto_user_S5-0198-28-run57_lamv4adjdnmt3-srsf2-tp53-20170410_133':78,
##	'Auto_user_S5-0198-28-run57_lamv4adjdnmt3-srsf2-tp53-20170410_133_tn':79,
#	'Auto_user_S5-0198-29-NGS186bis_NGS187bis_Chef_SBT-colon-lung_v7_530_134':81,
#	'Auto_user_S5-0198-29-NGS186bis_NGS187bis_Chef_SBT-colon-lung_v7_530_134_tn':82,
#	'Auto_user_S5-0198-2-Run_validpanel_PGM151_Chef_SBT_colon_lung_v5_530_90':20,
#	'Auto_user_S5-0198-2-Run_validpanel_PGM151_Chef_SBT_colon_lung_v5_530_90_tn':21,
##	'Auto_user_S5-0198-30-RUN57BIS-LAMv4adjDNMTSRSF2-TP53_530-12-04-2017_135':83,
##	'Auto_user_S5-0198-30-RUN57BIS-LAMv4adjDNMTSRSF2-TP53_530-12-04-2017_135_tn':84,
	'Auto_user_S5-0198-32-NGS_186ter_187ter_35pM_Chef_SBT-colon-lung_v7_530_138':87,
	'Auto_user_S5-0198-32-NGS_186ter_187ter_35pM_Chef_SBT-colon-lung_v7_530_138_tn':88,
	'Auto_user_S5-0198-33-NGS_188bis_189bis_35pM_Chef_SBT-colon-lung_v7_530_139':89,
	'Auto_user_S5-0198-33-NGS_188bis_189bis_35pM_Chef_SBT-colon-lung_v7_530_139_tn':90,
	'Auto_user_S5-0198-34-NGS190_191_35pM_Chef_SBT-colon-lung_v7_530_141':91,
	'Auto_user_S5-0198-34-NGS190_191_35pM_Chef_SBT-colon-lung_v7_530_141_tn':92,
#	'Auto_user_S5-0198-35-TS-1-Benjamin-Oncomine_Focus_DNA_142':93,
#	'Auto_user_S5-0198-35-TS-1-Benjamin-Oncomine_Focus_DNA_142_tn':94,
#	'Auto_user_S5-0198-36-Run1_30pM_Lymphome_B_Oceane_SBT_530_145':96,
#	'Auto_user_S5-0198-36-Run1_30pM_Lymphome_B_Oceane_SBT_530_145_tn':97,
#	'Auto_user_S5-0198-37-TP1-Benjamin-Oncomine_Focus_DNA_144':98,
#	'Auto_user_S5-0198-37-TP1-Benjamin-Oncomine_Focus_DNA_144_tn':99,
##	'Auto_user_S5-0198-38-Run57ter-Chef_LAMv4-27-04-2017_146':100,
##	'Auto_user_S5-0198-38-Run57ter-Chef_LAMv4-27-04-2017_146_tn':101,
##	'Auto_user_S5-0198-39-run58-lamv4-tp53-27-4-17_147':102,
##	'Auto_user_S5-0198-39-run58-lamv4-tp53-27-4-17_147_tn':103,
#	'Auto_user_S5-0198-3-Run_validation_PGM159CQ160_40pM_Chef_SBT-colon-lung_v5_530_99':22,
#	'Auto_user_S5-0198-3-Run_validation_PGM159CQ160_40pM_Chef_SBT-colon-lung_v5_530_99_tn':23,
	'Auto_user_S5-0198-40-NGS192_193_35pM_Chef_SBT-colon-lung_v7_530_148':104,
	'Auto_user_S5-0198-40-NGS192_193_35pM_Chef_SBT-colon-lung_v7_530_148_tn':105,
#	'Auto_user_S5-0198-41-TP2_Benjamin-Oncomine_Focus_DNA_149':106,
#	'Auto_user_S5-0198-41-TP2_Benjamin-Oncomine_Focus_DNA_149_tn':107,
#	'Auto_user_S5-0198-42-TS2_Benjamin-Oncomine_Focus_DNA_150':108,
#	'Auto_user_S5-0198-42-TS2_Benjamin-Oncomine_Focus_DNA_150_tn':109,
#	'Auto_user_S5-0198-43-Temoins_patients_Benjamin-Oncomine_Focus_DNA_151':110,
#	'Auto_user_S5-0198-43-Temoins_patients_Benjamin-Oncomine_Focus_DNA_151_tn':111,
#	'Auto_user_S5-0198-44-TS3_Benjamin-Oncomine_Focus_DNA_152':112,
#	'Auto_user_S5-0198-44-TS3_Benjamin-Oncomine_Focus_DNA_152_tn':113,
#	'Auto_user_S5-0198-45-TS4_Benjamin-Oncomine_Focus_DNA_153':114,
#	'Auto_user_S5-0198-45-TS4_Benjamin-Oncomine_Focus_DNA_153_tn':115,
	'Auto_user_S5-0198-46-NGS194_195_35pM_Chef_SBT-colon-lung_v7_530_155':116,
	'Auto_user_S5-0198-46-NGS194_195_35pM_Chef_SBT-colon-lung_v7_530_155_tn':117,
#	'Auto_user_S5-0198-47-Patients2_Benjamin-Oncomine_Focus_DNA_154':118,
#	'Auto_user_S5-0198-47-Patients2_Benjamin-Oncomine_Focus_DNA_154_tn':119,
	'Auto_user_S5-0198-48-NGS196_197_35pM_Chef_SBT-colon-lung_v7_530_156':122,
	'Auto_user_S5-0198-48-NGS196_197_35pM_Chef_SBT-colon-lung_v7_530_156_tn':123,
	'Auto_user_S5-0198-49-NGS198_199_35pM_Chef_SBT-colon-lung_v7_530_158':124,
	'Auto_user_S5-0198-49-NGS198_199_35pM_Chef_SBT-colon-lung_v7_530_158_tn':125,
##	'Auto_user_S5-0198-4-Run52-Chef-S5_copucageLAMv4-tp53_530-17-01-2017_100':24,
##	'Auto_user_S5-0198-4-Run52-Chef-S5_copucageLAMv4-tp53_530-17-01-2017_100_tn':25,
##	'Auto_user_S5-0198-5-NGS165_TP53serie7_35pM_Chef_SBT-colon-lung_v5_530_105':26,
##	'Auto_user_S5-0198-5-NGS165_TP53serie7_35pM_Chef_SBT-colon-lung_v5_530_105_tn':27,
#	'Auto_user_S5-0198-6-NGS168_169_35pM_Chef_SBT-colon-lung_v5_530_108':31,
#	'Auto_user_S5-0198-6-NGS168_169_35pM_Chef_SBT-colon-lung_v5_530_108_tn':32,
##	'Auto_user_S5-0198-7-Run53-Chef_LAMv4_530-70pM-7-2-2017_109':33,
##	'Auto_user_S5-0198-7-Run53-Chef_LAMv4_530-70pM-7-2-2017_109_tn':34,
#	'Auto_user_S5-0198-8-NGS168_169_bis_35pM_Chef_SBT-colon-lung_v5_530_110':35,
#	'Auto_user_S5-0198-8-NGS168_169_bis_35pM_Chef_SBT-colon-lung_v5_530_110_tn':36,
	'Auto_user_S5-0198-9-NGS168_169_complete_35pM_Chef_SBT-colon-lung_v5_530_111':37,
	'Auto_user_S5-0198-9-NGS168_169_complete_35pM_Chef_SBT-colon-lung_v5_530_111_tn':38,
#	'Name_corrected_NGS_186_187':80,
	'Run_validpanel_name_corrected':39,
#	'TEST_ALIGNMENT2':29,
#	'test_import_reanalyse':121,
#	'Test_new_Host_Name':13,
#	'Test-no-calibration-standard':51,
#	'TS-1-Benjamin-Oncomine_Focus_DNA_barcode_corrected':95,
}

# PGM 178, 181, 192, 211 : puces bleues
# PGM 200 : quel panel?

# PGM 207 et suivants : new panel 

HOSTNAME = 'http://chu06S501/rundb/api'
h = httplib2.Http()
h.add_credentials('ionadmin', 'ionadmin')
headers = {"Content-type": "application/json","Accept": "application/json"}
pluginName = "VariantBase"
pluginUpdate  = {"plugin": [pluginName]}

for run in Run2PK.keys():
	if '_tn' in run:
		print 'tn run, pass'
		continue
	PK = Run2PK[run]
	url = 'http://10.67.1.38' + '/rundb/api/v1/results/' + str(PK) + "/plugin/"
	print run
	print "-> Launching VariantBase..."
	resp, content = h.request(url, "POST", body=json.dumps(pluginUpdate),headers=headers )
	print resp
	print content
	
	time.sleep(120)
	PluginCompleted = False
	while not PluginCompleted:
		time.sleep(20)
		api_url = HOSTNAME + '/v1/pluginresult/?format=json&plugin__name=%s&result=%s' % (pluginName,str(PK))
		f = urllib2.urlopen(api_url)
		d = json.loads(f.read())
		PluginCompleted = d['objects'][0]['state']

