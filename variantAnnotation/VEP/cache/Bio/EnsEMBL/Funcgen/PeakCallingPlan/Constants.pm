package Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants;

=head1 Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants

Here is where all the constants are declared. To use them, add:

  use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );
  
to your module, then use:
  
  if ($peak_calling_strategy eq CALL_BROAD_PEAKS ) {
    convert_bam_to_bed();
  }

=cut

use base qw( Exporter );
use vars qw( @EXPORT_OK );
use strict;

our @EXPORT_OK = qw(

  CALL_BROAD_PEAKS
  CALL_NARROW_PEAKS
  CALL_TIGHT_PEAKS
  
  SKIP_IDR
  RUN_IDR_ON_BIOLOGICAL_REPLICATES 
  RUN_IDR_ON_TECHNICAL_REPLICATES  
  
  BAM_FORMAT
  BIGWIG_FORMAT
  BED_FORMAT
  
  SINGLE_END
  PAIRED_END
  
  ALIGNMENT_TYPE
  PEAK_CALLING_TYPE
  
  ALIGNMENT_ANALYSIS
  REMOVE_DUPLICATES_ANALYSIS
  IDR_ANALYSIS
  CALL_PEAKS_ANALYSIS
  CONVERT_BAM_TO_BED_ANALYSIS
  CONVERT_BAM_TO_BIGWIG_ANALYSIS
  
  SIGNAL_EXPERIMENT
  CONTROL_EXPERIMENT
  NO_CONTROL_FLAG
  
  TRUE
  FALSE
  
  ENSEMBL_SINGLE_END_ALIGNMENT_ANALYSIS
  ENSEMBL_PAIRED_END_ALIGNMENT_ANALYSIS
  ENSEMBL_HODGEPODGE_ALIGNMENT_ANALYSIS
  ENSEMBL_REMOVE_DUPLICATES_ANALYSIS
  
  ENSEMBL_BROAD_PEAK_CALLING_ANALYSIS
  ENSEMBL_NARROW_PEAK_CALLING_ANALYSIS_PERMISSIVE
  ENSEMBL_NARROW_PEAK_CALLING_ANALYSIS_DEFAULT
  ENSEMBL_TIGHT_PEAK_CALLING_ANALYSIS_DEFAULT
  
  NA

  SEGMENTATION_CLASS_NONE
  SEGMENTATION_CLASS_CTCF
  SEGMENTATION_CLASS_NO_CTCF
  SEGMENTATION_CLASS_CTCF_NO_H3K27AC
  SEGMENTATION_CLASS_NO_CTCF_NO_H3K27AC

);

our %EXPORT_TAGS = (all => \@EXPORT_OK);

use constant {

  SEGMENTATION_CLASS_NONE               => 'none',
  SEGMENTATION_CLASS_CTCF               => 'ctcf',
  SEGMENTATION_CLASS_NO_CTCF            => 'no_ctcf',
  SEGMENTATION_CLASS_CTCF_NO_H3K27AC    => 'ctcf_no_H3K27ac',
  SEGMENTATION_CLASS_NO_CTCF_NO_H3K27AC => 'no_ctcf_no_H3K27ac',

  CALL_BROAD_PEAKS  => 'CALL_BROAD_PEAKS',
  CALL_NARROW_PEAKS => 'CALL_NARROW_PEAKS',
  CALL_TIGHT_PEAKS  => 'CALL_TIGHT_PEAKS',

  SKIP_IDR                         => 'SKIP_IDR',
  RUN_IDR_ON_BIOLOGICAL_REPLICATES => 'RUN_IDR_ON_BIOLOGICAL_REPLICATES',
  RUN_IDR_ON_TECHNICAL_REPLICATES  => 'RUN_IDR_ON_TECHNICAL_REPLICATES',

  BAM_FORMAT    => 'bam',
  BIGWIG_FORMAT => 'bigwig',
  BED_FORMAT    => 'bed',
  
  SINGLE_END => 'single_end',
  PAIRED_END => 'paired_end',
  
  ENSEMBL_SINGLE_END_ALIGNMENT_ANALYSIS => 'bwa_samse',
  ENSEMBL_PAIRED_END_ALIGNMENT_ANALYSIS => 'bwa_sampe',
  ENSEMBL_HODGEPODGE_ALIGNMENT_ANALYSIS => 'bwa_hodgepodge',
  ENSEMBL_REMOVE_DUPLICATES_ANALYSIS    => 'remove_duplicates',
  
  ENSEMBL_BROAD_PEAK_CALLING_ANALYSIS             => 'ccat_histone',
  ENSEMBL_NARROW_PEAK_CALLING_ANALYSIS_PERMISSIVE => 'SWEmbl_R0005',
  ENSEMBL_NARROW_PEAK_CALLING_ANALYSIS_DEFAULT    => 'SWEmbl_default',
  ENSEMBL_TIGHT_PEAK_CALLING_ANALYSIS_DEFAULT     => 'SWEmbl_R0025',
  
  # Result type
  #
  ALIGNMENT_TYPE    => 'alignment',
  PEAK_CALLING_TYPE => 'peak',
  
  # Anlaysis type
  #
  ALIGNMENT_ANALYSIS             => 'alignment',
  REMOVE_DUPLICATES_ANALYSIS     => 'remove_duplicates',
  IDR_ANALYSIS                   => 'idr',
  CALL_PEAKS_ANALYSIS            => 'call_peaks',
  CONVERT_BAM_TO_BED_ANALYSIS    => 'convert bam to bed',
  CONVERT_BAM_TO_BIGWIG_ANALYSIS => 'create_bigwig',
  
  SIGNAL_EXPERIMENT  => 'signal',
  CONTROL_EXPERIMENT => 'control',
  
  NO_CONTROL_FLAG => 'no control',
  
  TRUE  => 'true',
  FALSE => 'false',
  
  NA => 'NA',
};

1;