#!/usr/bin/python
import sqlite3
import pysam
import json
import csv
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.assemblymapper
import hgvs.variantmapper
import hgvs.validator
import hgvs.exceptions
import hgvs.normalizer
import hgvs.exceptions
import xlrd

def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d

def reverse(seq):
	seq = seq.replace('A','W').replace('T','X').replace('G','Y').replace('C','Z')
	seq = seq.replace('W','T').replace('X','A').replace('Y','C').replace('Z','G')
	seq = seq[::-1]
	return seq

with open('/DATA/work/global_parameters.json', 'r') as g:
	global_param = json.load(g)
db_path = global_param['VariantBase']
db_con = sqlite3.connect(db_path)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()

#hdp = hgvs.dataproviders.uta.connect()
#hp = hgvs.parser.Parser()
#hn = hgvs.normalizer.Normalizer(hdp)
#hv = hgvs.validator.Validator(hdp)
#am37 = easyvariantmapper = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh37')


db_cur.execute("SELECT * FROM Variant")
db_variants = db_cur.fetchall()

for db_variant in db_variants:
	db_cur.execute("UPDATE Variant SET hgvs=NULL, hgvsInfo=NULL WHERE variantID='%s'" % (db_variant['variantID']))

#custom_drug_sensitivity = global_param['run_type']['SBT']['sbt_sensitivity']
#sensi = []
#sensitivity_xls = xlrd.open_workbook(custom_drug_sensitivity)
#sensitivity_sheet = sensitivity_xls.sheet_by_index(0)
#i = 1
#for row in range(sensitivity_sheet.nrows-1):
	#sensi.append((str(sensitivity_sheet.cell(row+1,5).value),str(sensitivity_sheet.cell(row+1,7).value),str(sensitivity_sheet.cell(row+1,8).value),str(sensitivity_sheet.cell(row+1,9).value))) # gene,c,p,sensi
	
#print sensi

#db_cur.execute("SELECT * FROM Variant")
#db_variants = db_cur.fetchall()


#for db_variant in db_variants:
	
	#if (db_variant['gene'] == 'MET') and (db_variant['region'] == 'intronic') and ((116411873 <= db_variant['genomicStart'] <= 116411902) or (116412044 <= db_variant['genomicStart'] <= 116412087)):
		#print db_variant['variantID']
		#print db_variant['actionability']
	
	
	## DRUG SENSITIVITY
	#for s in sensi:
		#if s[2] not in ['','p.(=)','p.?','p.=']:
			#if ((db_variant['gene'],db_variant['proteinDescription']) == (s[0],s[2])) or ((db_variant['gene'],db_variant['annovarProteinDescription']) == (s[0],s[2])) :
				#print '%s-%s-%s-%s' % (db_variant['variantID'],db_variant['gene'],db_variant['proteinDescription'],s[3])
				##db_variant['actionability'] = s[3]
				#break
		#if s[1] != '':
			#if ((db_variant['gene'],db_variant['transcriptDescription']) == (s[0],s[1])) or ((db_variant['gene'],db_variant['annovarTranscriptDescription']) == (s[0],s[1])):
				#print '%s-%s-%s-%s' % (db_variant['variantID'],db_variant['gene'],db_variant['transcriptDescription'],s[3])
				##db_variant['actionability'] = s[3]
				#break





#db_cur.execute("SELECT * FROM Gene")
#db_genes = db_cur.fetchall()
#for db_gene in db_genes:
	#actual_version = db_gene[9]
	#if str(actual_version) == 'None':
		#actual_version = 0
	#geneName = db_gene[0]
	#chromosome = db_gene[1]
	#nc = NC_data[chromosome]
	#transcript = db_gene[8]
	#gene_tx = hdp.get_tx_for_gene(geneName)
	#old_v = 0
	#for item in gene_tx:
		#if (item[3].split('.')[0] == transcript) and (item[4] == nc):
			#v = max(old_v,int(item[3].split('.')[-1]))
	#print "- %s:%s\t\t%s->%s" % (geneName,transcript,actual_version,v)
	#if actual_version != v:
		#print "- %s:%s\t\t%s->%s" % (geneName,transcript,actual_version,v)
		#db_cur.execute("UPDATE Gene SET transcriptVersion=%s WHERE GeneName='%s'" % (v,geneName))
	
##db_con.commit()
#db_con.close()
#exit()

# GET VARIANTS HGVS-ED WITHOUT GDESC
#db_cur.execute("SELECT * FROM Variant")
#db_variants = db_cur.fetchall()
#db_variant_corrected = []
#for db_variant in db_variants:
	#variantID = db_variant[0]
	#hgvs = db_variant[33]
	#if (hgvs != 'checked') and (hgvs != 'warning') and (str(hgvs) != 'None'):
		#db_cur.execute("SELECT * FROM Variant WHERE variantID='%s'"%hgvs)
		#db_variant_to_correct = db_cur.fetchone()
		#if str(db_variant_to_correct[35]) == 'None':
			#db_variant_corrected.append(variantID)
	
#db_cur.execute("SELECT * FROM Variant")
#db_variants = db_cur.fetchall()

##go = False
##z = 0
#for db_variant in db_variants:
	##z += 1
	#variantID = db_variant[0]
	##if variantID == 'chr4:106157606-106157607:->A':
		##print z
		##go = True
	##if not go:
		##continue	
	#genomeBuild = db_variant[1]
	#chromosome = db_variant[2]
	#genomicStart = db_variant[3]
	#genomicStop = db_variant[4]
	#ref = db_variant[5]
	#alt = db_variant[6]
	#variantType = db_variant[7]
	#genomicDescription = db_variant[35].split(':')[-1]
	#nc = NC_data[chromosome]
	#gene = db_variant[8]
	#db_cur.execute("SELECT * FROM Gene WHERE geneName='%s'" % gene)
	#db_gene = db_cur.fetchone()
	#transcript = "%s.%s" % (db_gene[8],db_gene[9])

	#var_g = '%s:%s' % (nc,genomicDescription)
	## STEP 1 : PARSE
	#var_g1 = hp.parse_hgvs_variant(var_g)
	#print "- %s (%s:%s-%s:%s>%s)" % (var_g1,chromosome,genomicStart,genomicStop,ref,alt)
	## STEP 2 : VALIDATE
	#try:
		#hv.validate(var_g1)
	#except hgvs.exceptions.HGVSError as e:
		#print (e)
		#db_cur.execute("UPDATE Variant SET hgvs='warning', hgvsWarning='%s' WHERE variantID='%s'" % (str(e),variantID))
		#continue
	## STEP 3 : NORMALIZE
	#var_g2 = hn.normalize(var_g1)
	#if var_g2 != var_g1:
		#var_g2_start = str(var_g2.posedit.pos.start)
		#var_g2_end = str(var_g2.posedit.pos.end)
		#if 'dup' in str(var_g2):
			#var_g2_start = str(var_g2.posedit.pos.start)
			#var_g2_end = str(int(str(var_g2.posedit.pos.start))+1)
			#var_g2_ref = '-'
			#var_g2_alt = str(var_g2.posedit.edit.ref)
			#variantType = 'DUP'
		#elif 'inv' in str(var_g2):
			#var_g2_ref = str(var_g2.posedit.edit.ref)
			#var_g2_alt = reverse(str(var_g2.posedit.edit.ref))
			#variantType = 'INV'
		#else:
			#var_g2_ref = str(var_g2.posedit.edit.ref)
			#var_g2_alt = str(var_g2.posedit.edit.alt)
			#if var_g2_ref == '-':
				#variantType = 'INS'
			#elif var_g2_alt == '-':
				#variantType = 'DEL'
			#elif len(var_g2_ref) > 1 or len(var_g2_alt) > 1: 
				#variant_type = 'DELINS'
			#else:
				#variant_type = 'SNV'
		#if var_g2_ref == 'None':
			#var_g2_ref = '-'
		#if var_g2_alt == 'None':
			#var_g2_alt = '-'
		#print "\t- normalyzed : %s (%s:%s-%s:%s>%s)" % (var_g2,chromosome,var_g2_start,var_g2_end,var_g2_ref,var_g2_alt)
		#if (genomicStart != var_g2_start) or (genomicStop != var_g2_end):
			#bad_variantID = variantID
			#variantID = '%s:%s-%s:%s>%s' % (chromosome,var_g2_start,var_g2_end,var_g2_ref,var_g2_alt)
			#genomicDescription = "%s:%s" % (chromosome,str(var_g2).split(':')[-1])
			#print "\t- HGVS shift, variantID is now %s" % variantID 
			#db_cur.execute("SELECT * FROM Variant WHERE variantID='%s'"%variantID)
			#if db_cur.fetchone() is None:
				#print "\t\t- creating new entry in db"
				#db_cur.execute("INSERT INTO Variant (variantID, genomeBuild, chromosome, genomicStart, genomicStop, referenceAllele, alternativeAllele, variantType, gene, genomicDescription) VALUES ('%s','%s','%s',%s, %s,'%s','%s','%s','%s','%s')" % (variantID, genomeBuild, chromosome, var_g2_start, var_g2_end, var_g2_ref, var_g2_alt, variantType,gene,genomicDescription))
			#db_cur.execute("UPDATE Variant SET hgvs='%s' WHERE variantID = '%s'" % (variantID,bad_variantID))
			#db_cur.execute("SELECT * FROM VariantMetrics WHERE variant = '%s'" % bad_variantID)
			#varmetrics = db_cur.fetchall()
			#for varmetric in varmetrics:
				#VariantMetricsID = varmetric[0]
				#db_cur.execute("UPDATE VariantMetrics SET variant = '%s' WHERE VariantMetricsID = '%s'" % (variantID,VariantMetricsID))
			#db_con.commit()
		#elif (ref != var_g2_ref) or (alt != var_g2_alt):
			#print "\t- WARNING : no change in position but REF ALT different : %s>%s => %s>%s" % (ref,alt,var_g2_ref,var_g2_alt)
	## STEP 4 : g_to_c
	#try:
		#c = am37.g_to_c(var_g2,transcript)
		#c_pos = str(c).split(':')[-1]
		#print "\t- %s" % c_pos
	#except hgvs.exceptions.HGVSError as e:
		#print (e)
		#db_cur.execute("UPDATE Variant SET hgvs='warning', hgvsWarning='%s' WHERE variantID='%s'" % (str(e),variantID))
		#db_con.commit()
		#continue
	#db_cur.execute("UPDATE Variant SET transcriptDescription='%s', hgvs='checked' WHERE variantID='%s'" % (c_pos,variantID))
	#db_con.commit()
	## STEP 5 : c_to_p
	#try:
		#p = am37.c_to_p(c)
		#p_pos = str(p).split(':')[-1]
		#print "\t- %s" % p_pos
	#except hgvs.exceptions.HGVSError as e:
		#print (e)
		#db_cur.execute("UPDATE Variant SET hgvs='warning', hgvsWarning='%s' WHERE variantID='%s'" % (str(e),variantID))
		#db_con.commit()
		#continue
	#db_cur.execute("UPDATE Variant SET proteinDescription='%s' WHERE variantID='%s'" % (p_pos,variantID))
	#db_con.commit()
#db_con.close()
#exit()


#db_cur.execute("SELECT * FROM Gene")
#gene_data = db_cur.fetchall()

#for i in range(len(gene_data)):
	#geneName = gene_data[i][0]
	#nm = gene_data[i][8]
	#refGene_file.seek(0)
	#for rgline in refGene_reader:
		#if rgline[1] == nm:
			#strand = rgline[3]
			#if strand == '+':
				#strand = 'forward'
			#else:
				#strand = 'reverse'
			#transcriptionStart = rgline[4]
			#transcriptionStop = rgline[5]			# rgline[6] et rgline[7] sont CodingRegionStart et codingRegionStop
			#exons = int(rgline[8])
			#exonsStart = rgline[9]
			#exonsStop =  rgline[10]
			#if exonsStart.endswith(','):
				#exonsStart = exonsStart[:-1]
			#if exonsStop.endswith(','):
				#exonsStop = exonsStop[:-1]			
			#if strand == 'reverse':
				#exonsStart = exonsStart.split(',')[::-1] 	# pour mettre dans l'odre des exons : 1, 2, 3..
				#exonsStop = exonsStop.split(',')[::-1]
				#exonsStart = ','.join(exonsStart)
				#exonsStop = ','.join(exonsStop)
			#break 
				
	#db_cur.execute("UPDATE Gene SET exons=%s, exonsStart='%s', exonsStop='%s' WHERE geneName = '%s'" % (exons,exonsStart,exonsStop,geneName))

#db_cur.execute("SELECT * FROM Gene")
#gene_data = db_cur.fetchall()



# bad variant to correct
#shit_variants = [
#'chr7:55242469-55242486:TAAGAGAAGCAACATCTC>-',
#'chr9:133729624-133729627:GGTG>AGGTGA',
#'chr9:133738339-133738419:AAGCTGGGCGGGGGCCAGTACGGGGAGGTGTACGAGGGCGTGTGGAAGAAATACAGCCTGACGGTGGCCGTGAAGACCTTG>-',
#'chr9:133748422-133748493:CAGAGATCTTGCTGCCCGAAACTGCCTGGTAGGGGAGAACCACTTGGTGAAGGTAGCTGATTTTGGCCTGAG>-',
#'chr11:534331-534440:AGGGGCCTGCGGCCCGGGGTCCTCCTACAGGGTCTCCTGCCCCACCTGCCAAGGAGGGCCCTGCTCAGCCAGGCCCAGGCCCAGCCCCAGGCCCCACAGGGCAGCTGCTG>-',
#'chr7:116411842-116411914:AAGTCTCCTGGGGCCCATGATAGCCGTCTTTAACAAGCTCTTTCTTTCTCTCTGTTTTAAGATCTGGGCAGTG>-',
#'chr1:2489255-2489358:TGCCCCAAGTGCAGTCCAGGTAGGTGCAGCCCTTTGGCGGGCCAGCTCTGTGGGCCGAGGGCAGACACTCTTGCCCCCTTCTGCCCCAGACACCCCTGTGTTCT>-',
#'chr10:31816312-31816382:GTAAAAACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACACAA>AGTAAAAACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACAAC',
#'chr10:31816312-31816382:GTAAAAACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACACAA>AGTAAAAACTAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACAAC',
#'chr10:31816312-31816381:GTAAAAACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACACA>AGTAAAAACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACAAC',
#'chr17:7753023-7753135:TGCTGACGTCGTGCGCGCCAGCAGGTGAGTCGGCTGCCTGCTTGCTTGTCCCGGAGACAGGCTCCCTTCCCCCATCACCCTGATGTTTCTGTCTTCATTCCACTGTCCTCCCG>CGCCGATGTCATGCGTGCCAGCAAGTGAGTTCGGAGGCTGCCGCACTTGTCCTGGGGCTAGGCTTGTGCTCGTTACCCCATGTTCCTGACCTCATTCTATGGTCTTCTCA',
#'chr17:7751226-7751329:TCACACCCCTCCCACTCCCCCAACCCCAACCACCAGCAGTAGCAACAGCAACAGTGGCAGCCACAGCAGCAGCCCTGCTGGGCCTGTGTCCTTTCCCCCACCAC>GCATAACCCTCCCATTCCCCCAACCACCACCAGCAGCAGCAGCAGCAGCAACAGCCACAGCAGTAGTCCTACTGGGCCGGTGCCCTTTCCACCACCCT',
#'chr10:31816312-31816381:GTAAAAACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACACA>AGTAAAAACTAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACAAC',
#'chr10:31816312-31816382:GTAAAAACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACACAA>AGTAAAAACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACACAAC',
#'chr10:31816312-31816381:GTAAAAACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACACA>AGTAAAAAACTAAAAAAATACAAAAATACAAAACACACACACACACACACACACACACACACACACACACAAC',
#'chr10:31816312-31816382:GTAAAAACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACACAA>AGTAAAAACTAAAAAATACAAAAATACAAAACACACACACACACACACACACACACACACACACAAC',
#'chr10:31816319-31816382:CTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACACAA>ACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACAAC',
#'chr20:39788406-39788500:TTAAGGGGGTAGAGGAGGTAGAGGATAGTTAGGGGAATGCCTGCTGGCTCCTGCCCAGTGGGAGGTATGTGCCCTCGGGGCAGCTATTGATACCT>ACATGGGCTGAGAGGAGGTGGGGTGTAGGTAAGGGGAGCCCAAGCACTTGCTCTTGTGGGAGCTGTGTGCCCTTGGTTCAGCCACTGCTC',
#'chr10:31816312-31816381:GTAAAAACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACACA>AGTAAAAACTAAAAAAATACAAAAATACAAAACACACACACACACACACACACACACACACACACACACAAC',
#'chr10:31816312-31816382:GTAAAAACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACACAA>AGTAAAAACTAAAAAAATACAAAAATACAAAACACACACACACACACACACACACACACACACACAAC',
#]
#for sv in shit_variants:
	#db_cur.execute("SELECT * FROM VariantMetrics WHERE variant='%s'" % sv)
	#varmetrics = db_cur.fetchall()
	#for varmetric in varmetrics:
		#db_cur.execute("DELETE FROM VariantMetrics WHERE VariantMetricsID='%s'" % varmetric[0])
	#db_cur.execute("DELETE FROM Variant WHERE variantID='%s'" % sv)
#db_con.commit()

## TEST UNIQUE
#db_cur.execute("SELECT * FROM Variant")
#all_variants = db_cur.fetchall()
#unique = {}
#for i in range(len(all_variants)):
	#genomicDescription = genomicStart = all_variants[i][35]
	#if genomicDescription not in unique :
		#unique[genomicDescription] = all_variants[i][0]
	#else:
		#print "- %s & %s" % (all_variants[i][0],unique[genomicDescription])

## TEST VARIANTTYPE
#for i in range(len(all_variants)):
	#variantID = all_variants[i][0]
	#genomicStart = all_variants[i][3]
	#genomicStop = all_variants[i][4]
	#ref = all_variants[i][5]
	#alt = all_variants[i][6]
	#variantType = all_variants[i][7]
	
	##if (ref != '-' and len(ref) == 1) and (alt != '-' and len(alt) == 1):
		##if variantType != 'SNV':
			##print variantID
			##print variantType
			##db_cur.execute("UPDATE Variant SET variantType='SNV' WHERE variantID='%s'" % variantID)
	#if variantType == 'DUP' and (genomicStart != genomicStop):
		#print variantID
		#print all_variants[i][13]

#db_con.commit()
#db_con.close()

#with open('/media/stuff/NMbatch.txt','w') as batch:
	#db_cur.execute("SELECT * FROM Variant")
	#db_variants = db_cur.fetchall()
	#for db_variant in db_variants:
		#gene = db_variant[8]
		#genomicDescription = db_variant[35]
		#gg = genomicDescription.split(':')[-1]
		#db_cur.execute("SELECT nc FROM Gene WHERE geneName='%s'" % gene)
		#nc = db_cur.fetchone()[0]
		#batch.write("%s:%s\n"%(nc,gg))
#exit()
	

## COMPARAISON MUTALYZER RESULTS
#with open('/DATA/work/variantBase/reannotate/mutalyzer_all_results.txt','r') as mutalyzer:
	#mutalyzer.next()
	#for line in mutalyzer:
		#line = line.split('\t')
		#genomicDescription = line[0]
		#db_cur.execute("SELECT * FROM Variant WHERE genomicDescription='%s'"%genomicDescription)
		#db_variant = db_cur.fetchone()
		#if db_variant is None:
			#print "- (%s not found in DB)" % genomicDescription
			#continue
		#variantID = db_variant[0]
		#chromosome = db_variant[2]
		#genomicStart = db_variant[3]
		#genomicStop = db_variant[4]
		#ref = db_variant[5]
		#alt = db_variant[6]
		#gene = db_variant[8]
		#db_cur.execute("SELECT * FROM Gene WHERE geneName='%s'"%gene)
		#db_gene = db_cur.fetchone()
		#strand = db_gene[4]
		##print "- %s" % variantID
		#message = line[3]
		#transcriptDescription = line[4]
		#posconv_transcriptDescription = line[2]
		#proteinDescription = line[5]
		#correct_genomicDescription = line[7]
		#nm = line[6]
		#if genomicDescription != correct_genomicDescription:
			#if 'delins' in correct_genomicDescription:
				#if '_' in correct_genomicDescription:
					#correct_start = int(correct_genomicDescription.split('g.')[-1].split('_')[0])
					#correct_stop = int(correct_genomicDescription.split('_')[-1].split('delins')[0])
				#else:
					#correct_start = int(correct_genomicDescription.split('g.')[-1].split('delins')[0])
					#correct_stop = int(correct_genomicDescription.split('g.')[-1].split('delins')[0])
				#correct_ref = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chromosome,correct_start,correct_stop)).split('\n')[1]
				#correct_alt = correct_genomicDescription.split('delins')[-1]
			#elif 'dup' in correct_genomicDescription:
				#if strand == 'reverse':
					#if '_' in correct_genomicDescription:
						#correct_start = int(correct_genomicDescription.split('g.')[-1].split('_')[0])-1
						#correct_stop = int(correct_genomicDescription.split('g.')[-1].split('_')[0])
						#correct_ref = '-'
						#correct_alt = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chromosome,int(correct_genomicDescription.split('g.')[-1].split('_')[0]),int(correct_genomicDescription.split('_')[-1].split('dup')[0]))).split('\n')[1]
					#else:
						#correct_start = int(correct_genomicDescription.split('g.')[-1].split('dup')[0])-1
						#correct_stop = int(correct_genomicDescription.split('g.')[-1].split('dup')[0])
						#correct_ref = '-'
						#correct_alt = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chromosome,int(correct_genomicDescription.split('g.')[-1].split('dup')[0]),int(correct_genomicDescription.split('g.')[-1].split('dup')[0]))).split('\n')[1]
				#else:
					#if '_' in correct_genomicDescription:
						#correct_start = int(correct_genomicDescription.split('_')[-1].split('dup')[0])
						#correct_stop = int(correct_genomicDescription.split('_')[-1].split('dup')[0])+1
						#correct_ref = '-'
						#correct_alt = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chromosome,int(correct_genomicDescription.split('g.')[-1].split('_')[0]),int(correct_genomicDescription.split('_')[-1].split('dup')[0]))).split('\n')[1]
					#else:
						#correct_start = int(correct_genomicDescription.split('g.')[-1].split('dup')[0])
						#correct_stop = int(correct_genomicDescription.split('g.')[-1].split('dup')[0])+1
						#correct_ref = '-'
						#correct_alt = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chromosome,int(correct_genomicDescription.split('g.')[-1].split('dup')[0]),int(correct_genomicDescription.split('g.')[-1].split('dup')[0]))).split('\n')[1]
			#elif 'ins' in correct_genomicDescription:
				#correct_start = int(correct_genomicDescription.split('g.')[-1].split('_')[0])
				#correct_stop = int(correct_genomicDescription.split('_')[-1].split('ins')[0])
				#correct_ref = '-'
				#correct_alt = correct_genomicDescription.split('ins')[-1]
			#elif 'inv' in correct_genomicDescription:
				#correct_start = int(correct_genomicDescription.split('g.')[-1].split('_')[0])
				#correct_stop = int(correct_genomicDescription.split('_')[-1].split('inv')[0])
				#correct_ref = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chromosome,correct_start,correct_stop)).split('\n')[1]
				#correct_alt = reverse(correct_ref)
			#elif 'del' in correct_genomicDescription:
				#if '_' in correct_genomicDescription:
					#correct_start = int(correct_genomicDescription.split('g.')[-1].split('_')[0])
					#correct_stop = int(correct_genomicDescription.split('_')[-1].split('del')[0])
				#else:
					#correct_start = int(correct_genomicDescription.split('g.')[-1].split('del')[0])
					#correct_stop = int(correct_genomicDescription.split('g.')[-1].split('del')[0])
				#correct_ref = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chromosome,correct_start,correct_stop)).split('\n')[1]
				#correct_alt = '-'
			#else: #SNV
				#correct_start = int(correct_genomicDescription.split('g.')[-1].split('A')[0].split('T')[0].split('G')[0].split('C')[0])
				#correct_stop = int(correct_genomicDescription.split('g.')[-1].split('A')[0].split('T')[0].split('G')[0].split('C')[0])
				#correct_ref = correct_genomicDescription[-1:]
				#correct_alt = correct_genomicDescription[-3:-2]
		
			#if (correct_ref != ref) or (correct_alt!= alt):
				#if (len(correct_ref) != len(ref)) or (len(correct_alt) != len(alt)):
					#print "- %s (%s>%s) & %s (%s>%s)" % (genomicDescription, ref, alt, correct_genomicDescription, correct_ref, correct_alt)
			#if (correct_start != genomicStart) or (correct_stop!= genomicStop):
				#print "- %s (%s-%s) & %s (%s-%s)" % (genomicDescription, genomicStart, genomicStop, correct_genomicDescription, correct_start, correct_stop)
				
			 #TEST : ref alt different ET start stop different? 	
			#if (correct_ref != ref) or (correct_alt!= alt):
				#if (len(correct_ref) != len(ref)) or (len(correct_alt) != len(alt)):
					#continue
				#if (correct_start == genomicStart) and (correct_stop == genomicStop): 
					#print "- %s : %s-%s:%s>%s" % (genomicDescription, genomicStart, genomicStop, ref, alt)
					#print "& %s : %s-%s:%s>%s" % (correct_genomicDescription, correct_start, correct_stop, correct_ref, correct_alt)
					

# VERIFICATION REF ET ALT PAR RAPPORT A HG19
#db_cur.execute("SELECT * FROM Variant")
#all_variants = db_cur.fetchall()
#for i in range(len(all_variants)):
	#variantID = all_variants[i][0]
	#gene = all_variants[i][8]
	#chromosome = all_variants[i][2]
	#genomicStart = all_variants[i][3]
	#genomicStop = all_variants[i][4]
	#ref = all_variants[i][5]
	#alt = all_variants[i][6]
	#if ref != '-':
		#refcheck = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chromosome,genomicStart,genomicStop)).split('\n')[1]
		#if ref != refcheck:
			#print "- %s : ref is %s, not %s" % (variantID,refcheck,ref)
			#print gene

# VERIFICATION REF ET ALT PAR RAPPORT A HG19
#db_cur.execute("SELECT * FROM Variant")
#all_variants = db_cur.fetchall()
#for i in range(len(all_variants)):
	#variantID = all_variants[i][0]
	#gene = all_variants[i][8]
	#chromosome = all_variants[i][2]
	#genomicStart = all_variants[i][3]
	#genomicStop = all_variants[i][4]
	#ref = all_variants[i][5]
	#alt = all_variants[i][6]
	#if ref != '-':
		#refcheck = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chromosome,genomicStart,genomicStop)).split('\n')[1]
		#if ref != refcheck:
			#print "- %s : ref is %s, not %s" % (variantID,refcheck,ref)
			#print gene
	#if '' in variantID:
		#refcheck = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chromosome,genomicStart,genomicStop)).split('\n')[1]
		#altcheck = variantID.split('>')[-1]
		#if alt != altcheck:
			#print "- %s : alt is %s, not %s" % (variantID,altcheck,alt)
			#db_cur.execute("UPDATE Variant SET alternativeAllele = '%s' WHERE variantID = '%s'" % (altcheck,variantID))
		#if ref != refcheck:
			#print "- %s : ref is %s, not %s" % (variantID,refcheck,ref)
			#db_cur.execute("UPDATE Variant SET referenceAllele = '%s' WHERE variantID = '%s'" % (refcheck,variantID))	
	#elif 'dup' in variantID:
		#altcheck = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chromosome,genomicStart,genomicStop)).split('\n')[1]
		#if alt != altcheck:
			#print "- %s : alt is %s, not %s" % (variantID,altcheck,alt)
			#db_cur.execute("UPDATE Variant SET alternativeAllele = '%s' WHERE variantID = '%s'" % (altcheck,variantID))
		#if ref != '-':
			#print "- %s : ref is %s, not %s" % (variantID,'-',ref)
			#db_cur.execute("UPDATE Variant SET referenceAllele = '%s' WHERE variantID = '%s'" % (refcheck,variantID))	
	#elif 'delins' in variantID:
		#refcheck = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chromosome,genomicStart,genomicStop)).split('\n')[1]
		#altcheck = variantID.split('delins')[-1]
		#if alt != altcheck:
			#print "- %s : alt is %s, not %s" % (variantID,altcheck,alt)
			#db_cur.execute("UPDATE Variant SET alternativeAllele = '%s' WHERE variantID = '%s'" % (altcheck,variantID))
		#if ref != refcheck:
			#print "- %s : ref is %s, not %s" % (variantID,refcheck,ref)
			#db_cur.execute("UPDATE Variant SET referenceAllele = '%s' WHERE variantID = '%s'" % (refcheck,variantID))	
	#elif 'ins' in variantID:
		#altcheck = variantID.split('ins')[-1]
		#if alt != altcheck:
			#print "- %s : alt is %s, not %s" % (variantID,altcheck,alt)
			#db_cur.execute("UPDATE Variant SET alternativeAllele = '%s' WHERE variantID = '%s'" % (altcheck,variantID))
		#if ref != '-':
			#print "- %s : ref is %s, not %s" % (variantID,'-',ref)
			#db_cur.execute("UPDATE Variant SET referenceAllele = '%s' WHERE variantID = '%s'" % (refcheck,variantID))	
	#elif 'del' in variantID:
		#refcheck = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chromosome,genomicStart,genomicStop)).split('\n')[1]
		#if alt != '-':
			#print "- %s : alt is %s, not %s" % (variantID,'-',alt)
			#db_cur.execute("UPDATE Variant SET alternativeAllele = '%s' WHERE variantID = '%s'" % (altcheck,variantID))
		#if ref != refcheck:
			#print "- %s : ref is %s, not %s" % (variantID,refcheck,ref)
			#db_cur.execute("UPDATE Variant SET referenceAllele = '%s' WHERE variantID = '%s'" % (refcheck,variantID))	
	#else:
		#print "- %s : ?????" % (variantID)

db_con.close()
