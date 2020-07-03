Trimmer
===============
Prior to alignment, Trimmer processes the read sequences to identify and remove the adaptor sequences and also extracts molecular barcodes (for SureSelect XT HS2).

*Note: This jar was compiled using Java version 8. Please make sure your Java Runtime Environment is at least at version 8 by running the command "java -version".*
<br>
<br>

**Command-line syntax:**
```
java -jar trimmer.jar [mandatory options] [options] -fq1 <read1_filename> -fq2 <read2_filename>
```
<br>

**Required parameters:**


Parameter              | Description                        
---------              | -----------                         
`-fq1 <filename>`      | Read1 FASTQ file (Multiple files can be provided separated by a comma).
`-fq2 <filename>`      | Read2 FASTQ file (Multiple files can be provided separated by a comma).

***Note:***  *Even though -fq1 and -fq2 accept multiple files separated by a comma, the program will output results in a single file for each read.*

<br>
At least one of the following available library prep types is also
mandatory to set the correct adaptor sequences for trimming.

Mandatory Option| Library Prep Type           
------| ----------------------------
`-halo` | HaloPlex
`-hs`   | HaloPlexHS
`-xt`   | SureSelect XT, XT2, XT HS
`-v2`   | SureSelect XT HS2
`-qxt`  | SureSelect QXT

<br>

**Optional Parameters:**

| Option | Description |
| -------|-------------|
| `-minFractionRead <n>`   | Sets the minimum read length as a fraction of the original read length after trimming.<br>Value range permitted is 0 to 99. Default value is 30. |
| `-idee_fixe` | Indicates that the fastq files are in the older Illumina fastq format (v1.5 or earlier). In addition to handling the older style read names, this option also assumes that the base qualities are encoded using the Illumina v1.5+ Phred+64 format and will attempt to convert bases to Phred+33. |
| `-out_loc` | Directory path for output files.  |

<br>

**Usage Examples:**<br>
```
java -jar trimmer.jar \
     -fq1 ./ICCG-repl1_S1_L001_R1_001.fastq.gz,./ICCG-repl1_S1_L001_R1_002.fastq.gz
     -fq2 ./ICCG-repl1_S1_L001_R2_001.fastq.gz,./ICCG-repl1_S1_L001_R2_002.fastq.gz
     -halo -minFractionRead 50 -idee_fixe \
     -out_loc result/outputFastqs/
```
<br>

**Tags for SureSelect XT HS2:**

For SureSelect XT HS2 option, trimmed molecular barcodes (MBCs) will be annotated in the readname. These annotation tags are:
* BC:Z:*sample barcode*
* ZA:Z:*3 bases of MBC (first half of dual MBC) followed by 1 or 2 dark base(s)*
* ZB:Z:*3 bases of MBC (second half of dual MBC) followed by 1 or 2 dark base(s)*
* RX:Z:*first half of MBC + second half of MBC concatenated with a "-")*
* QX:Z:*base quality of sequence in RX:Z (concatenated with a space)* 

e.g.
`@D00266:1113:HTWK5BCX2:1:1102:9976:2206 BC:Z:CTACCGAA+AAGTGTCT ZA:Z:TTAGT ZB:Z:TCCT RX:Z:TTA-TCC QX:Z:DDD DDA`

note:
The MBC bases are masked as **N** and corresponding base qualities marked as **$** in some annotations if they are not recognized as a valid XT HS2 MBC.

e.g.
`@K00336:80:HW7GLBBXX:7:1115:1184:3688 BC:Z:CTACCGAA+AGACACTT ZA:Z:NNNNN ZB:Z:AAAGT RX:Z:NNN-AAA QX:Z:$$$ <AA`

**Output for SureSelect XT HS2:**

In SureSelect XT HS2 mode (-v2), for every two FASTQ files (read 1 FASTQ file and read 2 FASTQ file) the program outputs three compressed files: 
* trimmed read 1 FASTQ file (.fastq.gz)
* trimmed read 2 FASTQ file (.fastq.gz)
* MBC sequence file (.txt.gz).
