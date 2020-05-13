    ###############################################################
    #                 This is HotCount verion 0                   #
    #                                                             #
    # released under GNU General Public Licence (see LICENCE.md)  #
    #                                                             #
    #                                                             #
    #                                                             #
    ###############################################################

 # This program is free software: you can redistribute it and/or modify
 # it under the terms of the GNU General Public License as published by
 # the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.
 #
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 #
 # You should have received a copy of the GNU General Public License
 # along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 # WARNING: doMreps use a slightly modified version of mreps. 
 # The original version should be downloaded from 
 #  http://mreps.univ-mlv.fr/.
 # The aim of the modification was to have a less verbose output to 
 # accelerate the processing (a more specific output for fastq data)
 # 
 # mreps is released by the author under the GNU GPL licence 
 #  (see LICENCE.md)
version=0
## usage : Rscript do_it.R design.txt CTL MUT {MUT}
##
## usage : Rscript do_it.R count.txt ALL delC dupC

## 
## Compute a few stats on all count data
## Warning : will substract MUT from CTL (CTL must be 'ALL')

## at least one count for each sample muste be over this threshold
MIN_COUNTS=100

cat ("## script do_it.R running in version", version,"\n");
cat ("## date is", date (), "\n");
args=commandArgs(TRUE)
cat ("openning file", args[1], "...");
csv=read.csv (args[1], sep=" ", header=TRUE, as.is=TRUE, comment.char="#");
cat ("done\n");

if (args[2] %in% colnames(csv)) {
  ctl=which (args[2] == colnames(csv));
} else {
  cat ("Erreur argument for control : ", args[2], " not found in count file", args[1], "\n");
  cat ("colnames = ", colnames(csv),"\n");
  exit (1)
}

if (length(args) < 3) {
  cat ("Erreur usage, need at least 3 arguments : Rscript do_it.R count.txt CTL MUT1 {MUT2}\n");
  exit (1);
}

for (mutationsArg in 3:length(args)) {
  if (args[mutationsArg] %in% colnames(csv)) {
    mut=which (args[mutationsArg] == colnames(csv));
  } else {
    cat ("Erreur argument for mutation : ", args[mutationsArg], " not found in count file", args[1], "\n");
    cat ("colnames = ", colnames(csv),"\n");
    exit (1)
  }
  cat ("\n##\n");
  cat ("## Testing mutation", args[mutationsArg], "against control", args[2],"\n");
  
  rows=c()
  for (sample in 1:nrow(csv)) {
    if (any (csv[sample,c(mut, ctl)] >= MIN_COUNTS)) {
      rows=c(rows, sample);
    } else {
      cat ("WARNING : Sample ", csv[sample,1], " doesn't have enough counts (", MIN_COUNTS, ") and has been excluded\n", sep="");
    }
  }
  
  cat ("\n");
  cat ("Fisher Matrice :\n");
  pMatrix = matrix(0,nrow=length(rows), ncol=length(rows));
  colnames(pMatrix)=rownames(pMatrix)=substr(csv[rows, 1], 1, 200);
  i=j=1;
  
  cat ("Samples", substr(csv[rows,1], 1,200), sep="\t"); cat ("\n");
  for (sample1 in rows) {
    cat (substr(csv[sample1, 1], 1,200), "\t", sep="");
    for (sample2 in rows) {
      m=matrix (c(csv[sample2, ctl]-csv[sample2, mut], csv[sample1, ctl]-csv[sample1, mut], csv[sample2, mut],csv[sample1, mut]), nrow=2);
      p=fisher.test (m, alternative="gr")$p.value;
                                        #    p=fisher.test (m, alternative="t")$p.value;
                                        #    p=chisq.test (m)$p.value;
      pMatrix[i,j]=p;
      cat (p, "\t", sep="");
      j=j+1;
    }
    cat ("\n");
    i=i+1;
    j=1;
  }
  
  cat ("\n");
  
  for (sample in 1:nrow(pMatrix)) {
    pMatrix[sample,sample] = NA;
  }
  
  for (pos in 1:min(6,nrow(pMatrix)-1)) {
    if (pos==1) {
      cat ("Who is significant with only", pos, "positive  sample  in this library (p<=0.001):");
    } else {
      cat ("Who is significant with only", pos, "positives samples in this library (p<=0.001):");
    }
    pass=0;
    for (i in 1:nrow(pMatrix)) {
      id=order (pMatrix[i,], decreasing=TRUE);
      if (pMatrix[i,id[pos]] <= 0.001) {
        cat ("\n", rownames(pMatrix)[i], "with a p-value of", pMatrix[i, id[pos]]);
        pass=pass+1;
      }
    }
    if (pass==0) {
      cat ("none !\n");
    } else {
      cat ("\n");
    }
    cat ("\n");
  }
}
