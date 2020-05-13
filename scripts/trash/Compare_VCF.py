#!/usr/bin/env python
import os
import json
import csv
from ion.plugin import *
from subprocess import *


class Compare_VCF(IonPlugin):
    """Plugin object which compare two variant caller output"""
    version = "0.1"
    envDict = dict(os.environ)
    json_dat = {}
    barcodeNames = []
    sampleNameLookup = {} # dictionary allowing us to get the sample name associated with a particular barcode name
    
    def getBarcodeNameFromFileName(self, fileName):
	for testBarcode in self.barcodeNames:
		testBarcode2 = testBarcode + '_'
		if testBarcode2 in fileName:
			barcodeName = testBarcode

	return barcodeName
	
    def compareConfig(self, config1, config2, html):
    	for key in config1.keys():
    		if key in config2.keys():
    			if (str(config1[key]) != str(config2[key])):
    				html.write('&#8594; %s : <b>VC A</b> = <span style=color:#B18904;>%s</span> &#8596; <b>VC B</b> = <span style=color:#B45F04;>%s</span> <br/>\n'% (key,str(config1[key]),str(config2[key])))
    			del config1[key],config2[key]
    		else:
    			html.write('&#8594; %s : <b>VC A</b> = <span style=color:#6A0888;>%s</span> &#8596; <b>VC B</b> = <span style=color:#6A0888;>%s</span> <br/>\n'% (key,str(config1[key]),"<i>Not found</i>"))#<span style=color:#6A0888;> %s</span>
    	for key in config2.keys():
    		html.write('&#8594; %s : <b>VC A</b> = <span style=color:#6A0888;>%s</span> &#8596; <b>VC B</b> = <span style=color:#6A0888;>%s</span> <br/>\n'% (key,"<i>Not found</i>",str(config2[key])))
    	return

    def launch(self):
        try:
		with open('startplugin.json', 'r') as fh:
			self.json_dat = json.load(fh)
	except:
		print 'Error reading plugin json.'

	try:
		htmlOut = open('%s/plugin_out/Compare_VCF_out/Compare_VCF_block.html'%self.envDict['ANALYSIS_DIR'], 'a+')
	except:
		htmlOut = open('Compare_VCF_block.html', 'w')
	htmlOut.write('<html><body>\n')
	
	if os.path.isfile(self.envDict['TSP_FILEPATH_BARCODE_TXT']):
		self.isBarcodedRun = True
			
	if not self.isBarcodedRun:
		print "This plugin only work for barcoded run for now. \n Exiting..."
		sys.exit(0)
		
	
	# set defaults for user options
	VCF_A = False
	VCF_B = False
	checkZygous = False
	checkNocall = False
	checkAbsent = False
	VCF1_config = False
	VCF2_config = False
	
	# Parse pluginconfig json.
	try:
		VCF_A = self.json_dat['pluginconfig']['select_VCF1']
		VCF_B = self.json_dat['pluginconfig']['select_VCF2']
		if (self.json_dat['pluginconfig']['zygousDiff'] == 'on'):
			checkZygous = True
		if (self.json_dat['pluginconfig']['nocallDiff'] == 'on'):
			checkNocall = True
		if (self.json_dat['pluginconfig']['absentDiff'] == 'on'):
			checkAbsent = True
		results_dir=self.json_dat['runinfo']['results_dir']
		VCF1_config = self.json_dat['pluginconfig']['VCF1_config']
		VCF2_config = self.json_dat['pluginconfig']['VCF2_config']
		
	except:
		print 'Warning: plugin does not appear to be configured'
		
	htmlOut.write('<b> VC A (base) :</b> variantCaller [%s]<br/>\n'%VCF_A)
	htmlOut.write('<b>   -</b> %s<br/>\n'%self.json_dat['pluginconfig']['LibraryLeft'])
	htmlOut.write('<b>   -</b> %s<br/>\n'%self.json_dat['pluginconfig']['TargetLeft'])
	htmlOut.write('<b>   -</b> %s<br/>\n'%self.json_dat['pluginconfig']['HotspotLeft'])
	htmlOut.write('<b>   -</b> %s<br/><br/>\n\n'%self.json_dat['pluginconfig']['ConfigLeft'])
	
	htmlOut.write('<b> VC B :</b> variantCaller [%s]<br/>\n'%VCF_B)
	htmlOut.write('<b>   -</b> %s<br/>\n'%self.json_dat['pluginconfig']['LibraryRight'])
	htmlOut.write('<b>   -</b> %s<br/>\n'%self.json_dat['pluginconfig']['TargetRight'])
	htmlOut.write('<b>   -</b> %s<br/>\n'%self.json_dat['pluginconfig']['HotspotRight'])
	htmlOut.write('<b>   -</b> %s<br/><br/>\n'%self.json_dat['pluginconfig']['ConfigRight'])
	
	#htmlOut.write('Check "Homozygous/Heterozygous" difference = %s<br/>\n'%checkZygous)
	#htmlOut.write('Check "No Call" difference = %s<br/>\n'%checkNocall)
	#htmlOut.write('Check "Absent" difference = %s<br/>\n'%checkAbsent)
	
	htmlOut.write('<br/><b> Differences in VariantCallers configurations :</b> <br/>\n')
	
	# Configuration comparison
	if (str(VCF1_config[0]) != str(VCF2_config[0])):
		htmlOut.write('&#8594; variantcaller version : <b>VC A</b> = %s &#8596; <b>VC B</b> = %s <br/>\n'% (str(VCF1_config[0]),str(VCF2_config[0])))
	
	self.compareConfig(VCF1_config[1],VCF2_config[1],htmlOut)
	self.compareConfig(VCF1_config[2],VCF2_config[2],htmlOut)
	## cleaning some data in meta parameters before compare
	for key in VCF1_config[3].keys():
		if type(VCF1_config[3][key]) is dict:
			VCF1_config[3].update(VCF1_config[3][key])
			del VCF1_config[3][key]
	for key in VCF2_config[3].keys():
		if type(VCF2_config[3][key]) is dict:
			VCF2_config[3].update(VCF2_config[3][key])
			del VCF2_config[3][key]
	## delete some pointless info
	#del VCF1_config[3]['targetloci_id'],VCF2_config[3]['targetloci_id'],VCF1_config[3]['targetloci_merge'],VCF2_config[3]['targetloci_merge'],VCF1_config[3]['targetregions_id'],VCF2_config[3]['targetregions_id'],VCF1_config[3]['targetregions_merge'],VCF2_config[3]['targetregions_merge']
	self.compareConfig(VCF1_config[3],VCF2_config[3],htmlOut)
	self.compareConfig(VCF1_config[4],VCF2_config[4],htmlOut)
	
	# Get bam filenames.
	with open(os.path.join(self.json_dat['runinfo']['basecaller_dir'], 'datasets_basecaller.json'), 'r') as f:
		json_basecaller = json.load(f)

	if isinstance(self.json_dat['plan']['barcodedSamples'],dict):
		samples = self.json_dat['plan']['barcodedSamples']
	else:
		samples = json.loads(self.json_dat['plan']['barcodedSamples'])
	
	bamPaths = []
	bams = []
	
	try:
		reference_path = self.envDict['TSP_FILEPATH_GENOME_FASTA']
	except:
		reference_path = ''
			
	for datum in json_basecaller['datasets']:
		if reference_path != '':
			tempPath = os.path.join(self.json_dat['runinfo']['alignment_dir'], datum['file_prefix']+'.bam')
		else:
			tempPath = os.path.join(self.json_dat['runinfo']['basecaller_dir'], datum['file_prefix']+'.basecaller.bam')

		if os.path.exists(tempPath):
			bamPaths.append(tempPath)
			if datum['dataset_name'][datum['dataset_name'].rfind('/')+1:] != 'No_barcode_match' and '/' in datum['dataset_name']:
				bams.append(datum['dataset_name'])

	# get the list of 'valid' barcodes or samples (could be either depending on whether user altered names with run planning
	# and sort of hacky, but extract this from the BAM file names we just got above
	for bamFileName in bamPaths:
		barcodeName = bamFileName.split('/')[-1] # get the last part, just the name with no path (probably can use os method here too)
		barcodeName = barcodeName.split('_rawlib')[0] # get just the barcode part of the name
		# find a possible matching sample name
		for sampleItemName in samples:
			sampleItem = samples[sampleItemName]
			if barcodeName in sampleItem['barcodes']:
				self.sampleNameLookup[barcodeName] = sampleItemName
		if barcodeName in self.sampleNameLookup.keys():
			sampleName = self.sampleNameLookup[barcodeName]
		else:
			sampleName = ''
			self.sampleNameLookup[barcodeName] = '' # makes it much easier later to do the lookup with no conditional tests
		# MGD note: I considered setting blank sample names to the barcode name instead, but might not be what customer intended
		print 'BARCODE FOUND: %s SAMPLE ID: %s' % (barcodeName, sampleName)
		self.barcodeNames.append(barcodeName)
		
	htmlOut.write('<br/><b> RESULTS :</b> <br/>\n')
	htmlOut.write('<ul type="none"> legend :<LI> <b>&#43;</b> &#8594; only found in VCF B</LI><LI> <b>&#8722;</b> &#8594; only found in VCF A</LI><LI> <b>&#8800;</b> &#8594; Allele Call different from VCF A to VCF B</LI></ul><br/>\n')
	
	self.OUTPUT_DIR = os.path.join(self.envDict['ANALYSIS_DIR'], 'plugin_out', 'downloads')
	if not os.path.isdir(self.OUTPUT_DIR):
		cmd = Popen(['mkdir', self.OUTPUT_DIR], stdout=PIPE, env=self.envDict)
		out, err = cmd.communicate()
		print 'OUT: %s\nERR: %s'%(out, err)
	
	filename = "Compare_VCF_%s_%s.csv" % (VCF_A,VCF_B)
	csvfile = open(self.OUTPUT_DIR+"/"+filename,"wb")
	results = csv.writer(csvfile, delimiter =";")
	results.writerow(["Barcode","Sample Name","Comparison","Chrom","Position", "Ref", "Variant", "Allele Call"])
	# parse VCFs and compare
	for barcode in self.barcodeNames:
		srcVcfA = '%s/plugin_out/variantCaller_out.%s/%s/alleles_%s.xls' % (self.envDict['ANALYSIS_DIR'], VCF_A, barcode, barcode)
		srcVcfB = '%s/plugin_out/variantCaller_out.%s/%s/alleles_%s.xls' % (self.envDict['ANALYSIS_DIR'], VCF_B, barcode, barcode)
		
		a = csv.reader(open(srcVcfA,"rb"),delimiter="\t")
		b = csv.reader(open(srcVcfB,"rb"),delimiter="\t")
		rows_a = [row_a[0:5] for row_a in a]
		rows_b = [row_b[0:5] for row_b in b]
		
		#########################################################################################################################################
		# eviter bug du variantcaller de duplication actuel, deux mutations identiques avec alleles call differentes. Souvent entre novel/hotspot
		# si on trouve duplication de mutation avec allele call diff, on garde seulement celle homo/hetero sinon cela fausse les resultats
		toRemove = []
		for row1 in rows_a:
			for row2 in rows_a:
				if (row1 != row2) and (row1[0:4] == row2[0:4]) and (row1[4] != row2[4]):
					print row1[0],row1[1],row1[2],row1[3],row1[4]+ "  === DUP BUG IN VCFA ===  " + row2[0],row2[1],row2[2],row2[3],row2[4]
					if (row1[4] == "Homozygous" or row1[4] == "Heterozygous"):
						if row2 not in toRemove:
							toRemove.append(row2)
					else:
						if row1 not in toRemove:
							toRemove.append(row1)
		for line in toRemove:
			rows_a.remove(line)
			
		toRemove = []				
		for row1 in rows_b:
			for row2 in rows_b:
				if (row1[0:4] == row2[0:4]) and (row1[4] != row2[4]):
					print row1[0],row1[1],row1[2],row1[3],row1[4]+ "  === DUP BUG IN VCFB ===  " + row2[0],row2[1],row2[2],row2[3],row2[4]
					if (row1[4] == "Homozygous" or row1[4] == "Heterozygous"):
						if row2 not in toRemove:
							toRemove.append(row2)
					else:
						if row1 not in toRemove:
							toRemove.append(row1)
		for line in toRemove:
			rows_b.remove(line)		
		#########################################################################################################################################
		
        	only_b = []
        	only_a = []
        	different_a_b = []
        	for row in rows_b:
        		only = True
        		for row_a in rows_a:
        			if row[0:4] == row_a[0:4]:
        				only = False
            				if row[4] != row_a[4]:
            					d = row_a
            					d.append(row[4])
                				different_a_b.append(d)
                	if only:
                		only_b.append(row[0:5])
                				
                		
               	for row in rows_a:
            		if row[0:5] not in rows_b:
                		only_a.append(row[0:5])
                
                htmlToWrite = []

		for row in only_a:
       			if ((row[4] == 'Absent') and checkAbsent):
       				htmlToWrite.append('<LI><b>&#45;</b><span style=color:#61210B;> %s</span>, <span style=color:#FFBF00;> %s</span> : <span style=color:#DBA901;>%s/%s</span> <span style=color:#B45F04;>%s</span>\n' % (row[0],row[1],row[2],row[3],row[4]))
       				results.writerow([barcode,self.sampleNameLookup[barcode],"Only found in VC "+VCF_A,row[0],row[1],row[2],row[3],row[4]])
       			elif ((row[4] == 'Homozygous' or row[4] == 'Heterozygous') and checkZygous):
       				htmlToWrite.append('<LI><b>&#45;</b><span style=color:#61210B;> %s</span>, <span style=color:#FFBF00;> %s</span> : <span style=color:#DBA901;>%s/%s</span> <span style=color:#74DF00;>%s</span>\n' % (row[0],row[1],row[2],row[3],row[4]))
       				results.writerow([barcode,self.sampleNameLookup[barcode],"Only found in VC "+VCF_A,row[0],row[1],row[2],row[3],row[4]])
       			elif ((row[4] == 'No Call') and checkNocall):
       				htmlToWrite.append('<LI><b>&#45;</b><span style=color:#61210B;> %s</span>, <span style=color:#FFBF00;> %s</span> : <span style=color:#DBA901;>%s/%s</span> <span style=color:#FF8000;>%s</span>\n' % (row[0],row[1],row[2],row[3],row[4]))
       				results.writerow([barcode,self.sampleNameLookup[barcode],"Only found in VC "+VCF_A,row[0],row[1],row[2],row[3],row[4]])
       		for row in only_b:
       			if ((row[4] == 'Absent') and checkAbsent):
       				htmlToWrite.append('<LI><b>&#43;</b><span style=color:#61210B;> %s</span>, <span style=color:#FFBF00;> %s</span> : <span style=color:#DBA901;>%s/%s</span> <span style=color:#B45F04;>%s</span>\n' % (row[0],row[1],row[2],row[3],row[4]))
       				results.writerow([barcode,self.sampleNameLookup[barcode],"Only found in VC "+VCF_B,row[0],row[1],row[2],row[3],row[4]])
       			elif ((row[4] == 'Homozygous' or row[4] == 'Heterozygous') and checkZygous):
       				htmlToWrite.append('<LI><b>&#43;</b><span style=color:#61210B;> %s</span>, <span style=color:#FFBF00;> %s</span> : <span style=color:#DBA901;>%s/%s</span> <span style=color:#74DF00;>%s</span>\n' % (row[0],row[1],row[2],row[3],row[4]))
       				results.writerow([barcode,self.sampleNameLookup[barcode],"Only found in VC "+VCF_B,row[0],row[1],row[2],row[3],row[4]])
       			elif ((row[4] == 'No Call') and checkNocall):
       				htmlToWrite.append('<LI><b>&#43;</b><span style=color:#61210B;> %s</span>, <span style=color:#FFBF00;> %s</span> : <span style=color:#DBA901;>%s/%s</span> <span style=color:#FF8000;>%s</span>\n' % (row[0],row[1],row[2],row[3],row[4]))
				results.writerow([barcode,self.sampleNameLookup[barcode],"Only found in VC "+VCF_B,row[0],row[1],row[2],row[3],row[4]])
		for row in different_a_b:
       			if ((row[4] == 'Absent' or row[5] == 'Absent') and checkAbsent):
       				htmlToWrite.append('<LI><b>&#8800;</b><span style=color:#61210B;> %s</span>, <span style=color:#FFBF00;> %s</span> : <span style=color:#DBA901;>%s/%s</span> <span style=color:#B45F04;>%s</span> &#10142; <span style=color:#B45F04;>%s</span>\n' % (row[0],row[1],row[2],row[3],row[4],row[5]))
       				results.writerow([barcode,self.sampleNameLookup[barcode],"Allele Call diff from VC "+VCF_A+" to "+VCF_B,row[0],row[1],row[2],row[3],row[4]+" --> "+row[5]])
       				continue
       			elif ((row[4] == 'Homozygous' or row[4] == 'Heterozygous' or row[5] == 'Homozygous' or row[5] == 'Heterozygous') and checkZygous):
       				htmlToWrite.append('<LI><b>&#8800;</b><span style=color:#61210B;> %s</span>, <span style=color:#FFBF00;> %s</span> : <span style=color:#DBA901;>%s/%s</span> <span style=color:#74DF00;>%s</span> &#10142; <span style=color:#74DF00;>%s</span>\n' % (row[0],row[1],row[2],row[3],row[4],row[5]))
       				results.writerow([barcode,self.sampleNameLookup[barcode],"Allele Call diff from VC "+VCF_A+" to "+VCF_B,row[0],row[1],row[2],row[3],row[4]+" --> "+row[5]])
       				continue
       			elif ((row[4] == 'No Call' or row[5] == 'No Call') and checkNocall):
       				htmlToWrite.append('<LI><b>&#8800;</b><span style=color:#61210B;> %s</span>, <span style=color:#FFBF00;> %s</span> : <span style=color:#DBA901;>%s/%s</span> <span style=color:#FF8000;>%s</span> &#10142; <span style=color:#FF8000;>%s</span>\n' % (row[0],row[1],row[2],row[3],row[4],row[5]))
				results.writerow([barcode,self.sampleNameLookup[barcode],"Allele Call diff from VC "+VCF_A+" to "+VCF_B,row[0],row[1],row[2],row[3],row[4]+" --> "+row[5]])
				continue

               	if htmlToWrite:
               		htmlOut.write('<UL TYPE="disc"><LI><b> %s  -  %s </b>\n<UL TYPE="none">' % (barcode, self.sampleNameLookup[barcode]))
               		if (len(htmlToWrite) <= 3):
               			for line in htmlToWrite:
               				htmlOut.write(line)
               			htmlOut.write('</UL></UL>')
               		else:
		       		for line in htmlToWrite[0:3]:
		       			htmlOut.write(line)
		       		htmlOut.write("</UL><div class='encart'><p class='titre'>...</p><div class='panneau' id='panneau-%s' style='display: none;'><UL TYPE='none'>" % barcode)
		       		for line in htmlToWrite[3:]:
		       			htmlOut.write(line)
		       		htmlOut.write('</UL></div></div></UL>')	
	htmlOut.write('<br/><a href="../%s">%s</a><br>'%("downloads/"+filename, filename))
	
	# JavaScript pour les panneaux
	htmlOut.write('''
	<script type="text/javascript" src="/site_media/jquery/js/jquery-1.6.1.min.js"></script>
	<script type='text/javascript'>
	$(document).ready(function()
{
	$panneaux = $('div.panneau').hide();
	$('p.titre').each(function(i)
	{
		$this = $(this);
		ancre = $this.next($panneaux)[0].id;
		lien = $('<a>',
		{
			'href':				'#' + ancre,
			'aria-expanded':	'false',
			'aria-controls':	ancre
		});
		$this.wrapInner(lien);
	});
	$('p.titre > a').click(function() 
	{
		if ($(this).attr('aria-expanded') == 'false') 
		{
            $(this).attr('aria-expanded', true).parent().next($panneaux).show();
		} 
		else 
		{
            $(this).attr('aria-expanded', false).parent().next($panneaux).hide();
		}
		return false;
	}); 
});
	</script>
	''')
	
	htmlOut.close()
	csvfile.close()
        return True


if __name__ == "__main__":
    PluginCLI(Compare_VCF())
