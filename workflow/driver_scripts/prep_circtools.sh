#!/bin/bash
# @Author: Tobias Jakobi <tjakobi@arizona.edu>

#$ -N prep_circtools
#$ -cwd
#$ -pe smp 2
#$ -l h_vmem=10G

# Check if we have 2 arguments
if [ ! $# == 2 ]; then
  echo "Usage: $0 [STAR source dir] [circtools destination dir]"
  exit 1
fi

SRC=$(realpath "$1")
DEST=$(realpath "$2")

# Check that directories exist
if [ ! -d "$SRC" ]; then
  echo "Source directory $SRC does not exist!"
  exit 1
fi

if [ ! -d "$DEST" ]; then
  echo "Destination directory $DEST does not exist!"
  exit 1
fi

cd "$SRC" || exit 1

# Link files using GNU parallel
parallel ln -s "$SRC"/{1}/mate{2}/Chimeric.out.junction "$DEST"/{1}.mate{2}.Chimeric.out.junction ::: * ::: 1 2
parallel ln -s "$SRC"/{1}/mate{2}/Aligned.noS.bam "$DEST"/{1}.mate{2}.bam ::: * ::: 1 2 
parallel ln -s "$SRC"/{1}/mate{2}/Aligned.noS.bam.bai "$DEST"/{1}.mate{2}.bam.bai ::: * ::: 1 2

parallel ln -s "$SRC"/{1}/Chimeric.out.junction "$DEST"/{1}.Chimeric.out.junction ::: * 
parallel ln -s "$SRC"/{1}/Aligned.noS.bam "$DEST"/{1}.bam ::: * 
parallel ln -s "$SRC"/{1}/Aligned.noS.bam.bai "$DEST"/{1}.bam.bai ::: *
parallel ln -s "$SRC"/{1}/SJ.out.tab "$DEST"/{1}.SJ.out.tab ::: *

# Go back to the root directory
cd "$DEST" || exit 1

# Create metadata files
ls *.bam | grep -v bai | grep -v mate | grep -v bam_files > bam_files.txt
ls *.junction | grep Chimeric.out.junction | grep mate1 | grep -v fixed > mate1
ls *.junction | grep Chimeric.out.junction | grep mate2 | grep -v fixed > mate2
ls *.junction | grep Chimeric.out.junction | grep -v mate > samplesheet
