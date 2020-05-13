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

#!/bin/bash -f
version=0.0
## usage : ./do_it.sh design.txt *.fastq

d=`date`
echo "## script do_it.sh running in version $version"
echo "## date is $d"
design=$1;
shift;

echo -n "Sample"
while read line
do
    mut=`echo $line | sed -e s/=\.\*//g`
    echo -n " $mut";
done < $design
echo ""


for file in $*
do

    echo -n $(basename $file)
    total=0

    while read line  
    do  
	Forward=`echo $line | sed -e s/\.\*=//g`
	Reversed=`echo $Forward | rev`
	Reversed=`echo $Reversed | tr A W | tr C X | tr G Y | tr T Z`
	Reversed=`echo $Reversed | sed -e s/\(/V+/g | sed -e s/\+\)/S/g`
	Reversed=`echo $Reversed | tr W T | tr X G | tr Y C | tr Z A | tr V ")" | tr S "("`
	nb=$(grep -c -E "${Forward}"\|"${Reversed}" "$file")
	echo -n " $nb"
    done < $design

    echo ""

done
