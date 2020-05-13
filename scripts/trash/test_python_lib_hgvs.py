#!/usr/bin/env python
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.assemblymapper
import hgvs.variantmapper
import hgvs.validator
import hgvs.exceptions
import hgvs.normalizer
import hgvs.exceptions

hdp = hgvs.dataproviders.uta.connect()
hp = hgvs.parser.Parser()
hn = hgvs.normalizer.Normalizer(hdp)
hv = hgvs.validator.Validator(hdp)

#vm = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh37', alt_aln_method='splign')
am37 = easyvariantmapper = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh37')
am37no = easyvariantmapper = hgvs.assemblymapper.AssemblyMapper(hdp, normalize=False, assembly_name='GRCh37')

#gene_info = hdp.get_gene_info('PDGFRA')
#print "- gene_info : %s \n" % gene_info

#gene_tx = hdp.get_tx_for_gene('PDGFRA')
#print "- gene_tx : %s \n" % gene_tx

#tx_exons = hdp.get_tx_exons(tx_ac, alt_ac, alt_aln_method) # argument ??

#tx_identity = hdp.get_tx_identity_info('NM_001127500.1')
#print "- tx_identity : %s \n" % tx_identity

#tx_info = hdp.get_tx_info(tx_ac, alt_ac, alt_aln_method)

#tx_mapping_opt = hdp.get_tx_mapping_options('NM_006206.5')
#print "- tx_mapping_opt : %s \n" % tx_mapping_opt

#var_g = 'NC_000004.11:g.55141054_55141055delinsTG'
##var_g = 'NC_000004.11:g.55141054G>A' # FAUX

## STEP 1 : PARSE
#g = hp.parse_hgvs_variant(var_g)
#print "- parsed_variant_g : %s \n" % g

## STEP 2 : VALIDATE
#unvalid = False
#try:
	#hv.validate(g)
#except hgvs.exceptions.HGVSError as e:
	#print (e)
	#unvalid = True
	
## STEP 3 : NORMALIZE
#if not unvalid:
	#g = hn.normalize(g)
	#print "- parsed_variant_g_normalyzed : %s \n" % g

## STEP 4 : ???
## STEP 5 :  PROFIT

	##try:
		##transcripts = am37.relevant_transcripts(g)
		##print "- transcripts : %s \n" % transcripts
	##except Exception as e:
		##print (e)

	#try:
		#c = am37.g_to_c(g,'NM_006206.5')
		#print "- c : %s \n" % c
	#except hgvs.exceptions.HGVSError as e:
		#print (e)

#vc = hp.parse_hgvs_variant('NM_002524.4:c.186_188dup')
#vg1 = am37.c_to_g(vc)
#vc1 = am37.g_to_c(vg1,'NM_002524.4')

#vg2 = am37no.c_to_g(vc)
#vc2 = am37no.g_to_c('NC_000001.10:g.115256527_115256529dup','NM_002524.4')

#vg3 = hp.parse_hgvs_variant('NC_000001.10:g.115256525_115256526insTCT')
#vc3 = am37.g_to_c(vg3,'NM_002524.4')
#vg4 = am37no.c_to_g(vc3)

#vg3 = hp.parse_hgvs_variant('NC_000001.10:g.115256560G>T')
#vc3 = am37.g_to_c(vg3,'NM_002524.4')
#vg4 = am37no.c_to_g(vc3)

#g0 = hp.parse_hgvs_variant('NC_000004.11:g.55592185_55592186insGCCTAT')
#c = am37.g_to_c(vg3,'NM_000222.2')
#g = am37no.c_to_g(vc3)

g0 = hp.parse_hgvs_variant('NC_000005.9:g.170837547_170837548insTCTG')
c = am37.g_to_c(g0,'NM_002520.6')
g = am37no.c_to_g(c)
gbis = am37.c_to_g(c)

print '%s : start=%s, end=%s' % (g,g.posedit.pos.start.base,g.posedit.pos.end.base)
print '%s : start=%s, end=%s' % (gbis,gbis.posedit.pos.start.base,gbis.posedit.pos.end.base)
