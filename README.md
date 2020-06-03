# ngs-pipeline-bx

(work in progress)

ngs-pipeline-bx is NGS analysis pipeline dealing with sequencing data.

## Installation

### Download

```bash
git clone https://github.com/tbandres/ngs-pipeline-bx.git
cd ngs-pipeline-bx
```

### Dependencies

```bash
sudo apt-get install sqlite3 numpy samtools bedtools r-base \
 fonts-freefont-ttf python-backports.functools-lru-cache
sudo pip install openpyxl==2.4.8 docker matplotlib pandas Pilow
sudo R
> install.packages("stringr")

```

### Environnement variable

Define these variables in /etc/bashrc or /etc/bash.bashrc

```bash
export NGS_PIPELINE_BX_DIR=/home/t8illumy/ngs-pipeline-bx
export HGVS_SEQREPO_DIR=$NGS_PIPELINE_BX_DIR/seqrepo/2019-06-20
export UTA_DB_URL=postgresql://anonymous@localhost:15032/uta/uta_20180821
```

### Docker

```bash
docker pull broadinstitute/gatk
 docker run -dit --name gatk --restart always -v /media/n06lbth:/media/n06lbth -v /home/t8illumy:/home/t8illumy broadinstitute/gatk
```

### Python HGVS module

```bash
sudo apt-get install libpq-dev python-pip zlib1g-dev libbz2-dev liblzma-dev
sudo pip install hgvs biocommons.seqrepo
```
If firewall not blocking :
```bash
seqrepo -r /usr/local/share/seqrepo pull 
```
Else :
```bash
mkdir ~/ngs-pipeline-bx/hgvs
cd ~/ngs-pipeline-bx/hgvs 
wget -r --no-parent --reject="index.html*" http://54.201.113.125/seqrepo/2019-06-20/
# nettoyage
mv ~/ngs-pipeline-bx/54.201.113.125/2019-06-20 ~/ngs-pipeline-bx/seqrepo/
rm -r ~/ngs-pipeline-bx/seqrepo/54.201.113.125
export HGVS_SEQREPO_DIR=~/ngs-pipeline-bx/seqrepo/2019-06-20
# OU plutot, dans /etc/bash.bashrc avant "#if not running interactively..."
export HGVS_SEQREPO_DIR=~/ngs-pipeline-bx/seqrepo/2019-06-20
# en ensuite recharger
source ~/.bashrc
```

Install UTA db :
```bash
https://github.com/biocommons/uta/tree/master/misc/docker
#Port 50827:5432 is for python HGVS module
```

### Locales

Locales should be set in english UTF8 (if not already) :

```bash
sudo apt-get update
sudo dpkg-reconfigure locales 
# check with : locales -a
```


## Usage

```bash
python run_analysis.py --run /path/to/run
```


```bash
```
