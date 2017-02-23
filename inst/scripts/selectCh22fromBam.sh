#!/bin/bash
#
# Copy the u12Regions.bed file to $MAPPINGPATH and run the following scripts
# The folder used here was ~/analysis/zrsr2/mapping 
# samtools 0.1.18 was used

# Copying the u12Regions.bed file to MAPPINGPATH
cp $SCRIPTSPATH/u12Regions.bed $MAPPINGPATH/

cd $MAPPINGPATH
# Choosing Chr22 
samtools view -b -o SRR1691633_ZRSR2Mut.bam SRR1691633/accepted_hits.sorted.bam chr22
samtools view -b -o SRR1691634_ZRSR2Mut.bam SRR1691634/accepted_hits.sorted.bam chr22
samtools view -b -o SRR1691635_ZRSR2Mut.bam SRR1691635/accepted_hits.sorted.bam chr22
samtools view -b -o SRR1691636_ZRSR2Mut.bam SRR1691636/accepted_hits.sorted.bam chr22
samtools view -b -o SRR1691637_ZRSR2Mut.bam SRR1691637/accepted_hits.sorted.bam chr22
samtools view -b -o SRR1691638_ZRSR2Mut.bam SRR1691638/accepted_hits.sorted.bam chr22
samtools view -b -o SRR1691639_ZRSR2Mut.bam SRR1691639/accepted_hits.sorted.bam chr22
samtools view -b -o SRR1691640_ZRSR2Mut.bam SRR1691640/accepted_hits.sorted.bam chr22
samtools view -b -o SRR1691641_WT.bam SRR1691641/accepted_hits.sorted.bam chr22
samtools view -b -o SRR1691642_WT.bam SRR1691642/accepted_hits.sorted.bam chr22
samtools view -b -o SRR1691643_WT.bam SRR1691643/accepted_hits.sorted.bam chr22
samtools view -b -o SRR1691644_WT.bam SRR1691644/accepted_hits.sorted.bam chr22
samtools view -b -o SRR1691645_Normal.bam SRR1691645/accepted_hits.sorted.bam chr22
samtools view -b -o SRR1691646_Normal.bam SRR1691646/accepted_hits.sorted.bam chr22
samtools view -b -o SRR1691647_Normal.bam SRR1691647/accepted_hits.sorted.bam chr22
samtools view -b -o SRR1691648_Normal.bam SRR1691648/accepted_hits.sorted.bam chr22
# Sorting the bam files
samtools sort SRR1691633_ZRSR2Mut.bam SRR1691633_ZRSR2Mut
samtools sort SRR1691634_ZRSR2Mut.bam SRR1691634_ZRSR2Mut
samtools sort SRR1691635_ZRSR2Mut.bam SRR1691635_ZRSR2Mut
samtools sort SRR1691636_ZRSR2Mut.bam SRR1691636_ZRSR2Mut
samtools sort SRR1691637_ZRSR2Mut.bam SRR1691637_ZRSR2Mut
samtools sort SRR1691638_ZRSR2Mut.bam SRR1691638_ZRSR2Mut
samtools sort SRR1691639_ZRSR2Mut.bam SRR1691639_ZRSR2Mut
samtools sort SRR1691640_ZRSR2Mut.bam SRR1691640_ZRSR2Mut
samtools sort SRR1691641_WT.bam SRR1691641_WT
samtools sort SRR1691642_WT.bam SRR1691642_WT
samtools sort SRR1691643_WT.bam  SRR1691643_WT
samtools sort SRR1691644_WT.bam  SRR1691644_WT
samtools sort SRR1691645_Normal.bam SRR1691645_Normal
samtools sort SRR1691646_Normal.bam SRR1691646_Normal
samtools sort SRR1691647_Normal.bam SRR1691647_Normal
samtools sort SRR1691648_Normal.bam SRR1691648_Normal

#Indexing bam files
samtools index SRR1691633_ZRSR2Mut.bam
samtools index SRR1691634_ZRSR2Mut.bam
samtools index SRR1691635_ZRSR2Mut.bam
samtools index SRR1691636_ZRSR2Mut.bam
samtools index SRR1691637_ZRSR2Mut.bam
samtools index SRR1691638_ZRSR2Mut.bam
samtools index SRR1691639_ZRSR2Mut.bam
samtools index SRR1691640_ZRSR2Mut.bam
samtools index SRR1691641_WT.bam
samtools index SRR1691642_WT.bam
samtools index SRR1691643_WT.bam
samtools index SRR1691644_WT.bam
samtools index SRR1691645_Normal.bam
samtools index SRR1691646_Normal.bam
samtools index SRR1691647_Normal.bam
samtools index SRR1691648_Normal.bam

#Converting bam to sam files
samtools view -h -o SRR1691633_ZRSR2Mut.sam SRR1691633_ZRSR2Mut.bam
samtools view -h -o SRR1691634_ZRSR2Mut.sam SRR1691634_ZRSR2Mut.bam
samtools view -h -o SRR1691635_ZRSR2Mut.sam SRR1691635_ZRSR2Mut.bam
samtools view -h -o SRR1691636_ZRSR2Mut.sam SRR1691636_ZRSR2Mut.bam
samtools view -h -o SRR1691637_ZRSR2Mut.sam SRR1691637_ZRSR2Mut.bam
samtools view -h -o SRR1691638_ZRSR2Mut.sam SRR1691638_ZRSR2Mut.bam
samtools view -h -o SRR1691639_ZRSR2Mut.sam SRR1691639_ZRSR2Mut.bam
samtools view -h -o SRR1691640_ZRSR2Mut.sam SRR1691640_ZRSR2Mut.bam
samtools view -h -o SRR1691641_WT.sam SRR1691641_WT.bam
samtools view -h -o SRR1691642_WT.sam SRR1691642_WT.bam
samtools view -h -o SRR1691643_WT.sam SRR1691643_WT.bam
samtools view -h -o SRR1691644_WT.sam SRR1691644_WT.bam
samtools view -h -o SRR1691645_Normal.sam SRR1691645_Normal.bam
samtools view -h -o SRR1691646_Normal.sam SRR1691646_Normal.bam
samtools view -h -o SRR1691647_Normal.sam SRR1691647_Normal.bam
samtools view -h -o SRR1691648_Normal.sam SRR1691648_Normal.bam

# Correct the read names to identify the paired reads correctly
sed -ri -- 's/\.2\t|\.1\t/\t/' SRR1691633_ZRSR2Mut.sam
sed -ri -- 's/\.2\t|\.1\t/\t/' SRR1691634_ZRSR2Mut.sam
sed -ri -- 's/\.2\t|\.1\t/\t/' SRR1691635_ZRSR2Mut.sam
sed -ri -- 's/\.2\t|\.1\t/\t/' SRR1691636_ZRSR2Mut.sam
sed -ri -- 's/\.2\t|\.1\t/\t/' SRR1691637_ZRSR2Mut.sam
sed -ri -- 's/\.2\t|\.1\t/\t/' SRR1691638_ZRSR2Mut.sam
sed -ri -- 's/\.2\t|\.1\t/\t/' SRR1691639_ZRSR2Mut.sam
sed -ri -- 's/\.2\t|\.1\t/\t/' SRR1691640_ZRSR2Mut.sam
sed -ri -- 's/\.2\t|\.1\t/\t/' SRR1691641_WT.sam
sed -ri -- 's/\.2\t|\.1\t/\t/' SRR1691642_WT.sam
sed -ri -- 's/\.2\t|\.1\t/\t/' SRR1691643_WT.sam
sed -ri -- 's/\.2\t|\.1\t/\t/' SRR1691644_WT.sam
sed -ri -- 's/\.2\t|\.1\t/\t/' SRR1691645_Normal.sam
sed -ri -- 's/\.2\t|\.1\t/\t/' SRR1691646_Normal.sam
sed -ri -- 's/\.2\t|\.1\t/\t/' SRR1691647_Normal.sam
sed -ri -- 's/\.2\t|\.1\t/\t/' SRR1691648_Normal.sam

# Convert corrected sam to bam
samtools view -bS -o SRR1691633_ZRSR2Mut.bam SRR1691633_ZRSR2Mut.sam
samtools view -bS -o SRR1691634_ZRSR2Mut.bam SRR1691634_ZRSR2Mut.sam
samtools view -bS -o SRR1691635_ZRSR2Mut.bam SRR1691635_ZRSR2Mut.sam
samtools view -bS -o SRR1691636_ZRSR2Mut.bam SRR1691636_ZRSR2Mut.sam
samtools view -bS -o SRR1691637_ZRSR2Mut.bam SRR1691637_ZRSR2Mut.sam
samtools view -bS -o SRR1691638_ZRSR2Mut.bam SRR1691638_ZRSR2Mut.sam
samtools view -bS -o SRR1691639_ZRSR2Mut.bam SRR1691639_ZRSR2Mut.sam
samtools view -bS -o SRR1691640_ZRSR2Mut.bam SRR1691640_ZRSR2Mut.sam
samtools view -bS -o SRR1691641_WT.bam SRR1691641_WT.sam
samtools view -bS -o SRR1691642_WT.bam SRR1691642_WT.sam
samtools view -bS -o SRR1691643_WT.bam SRR1691643_WT.sam
samtools view -bS -o SRR1691644_WT.bam SRR1691644_WT.sam
samtools view -bS -o SRR1691645_Normal.bam SRR1691645_Normal.sam
samtools view -bS -o SRR1691646_Normal.bam SRR1691646_Normal.sam
samtools view -bS -o SRR1691647_Normal.bam SRR1691647_Normal.sam
samtools view -bS -o SRR1691648_Normal.bam SRR1691648_Normal.sam

# Convert Sort the corrected bam
samtools sort SRR1691633_ZRSR2Mut.bam SRR1691633_ZRSR2Mut
samtools sort SRR1691634_ZRSR2Mut.bam SRR1691634_ZRSR2Mut
samtools sort SRR1691635_ZRSR2Mut.bam SRR1691635_ZRSR2Mut
samtools sort SRR1691636_ZRSR2Mut.bam SRR1691636_ZRSR2Mut
samtools sort SRR1691637_ZRSR2Mut.bam SRR1691637_ZRSR2Mut
samtools sort SRR1691638_ZRSR2Mut.bam SRR1691638_ZRSR2Mut
samtools sort SRR1691639_ZRSR2Mut.bam SRR1691639_ZRSR2Mut
samtools sort SRR1691640_ZRSR2Mut.bam SRR1691640_ZRSR2Mut
samtools sort SRR1691641_WT.bam SRR1691641_WT
samtools sort SRR1691642_WT.bam SRR1691642_WT
samtools sort SRR1691643_WT.bam  SRR1691643_WT
samtools sort SRR1691644_WT.bam  SRR1691644_WT
samtools sort SRR1691645_Normal.bam SRR1691645_Normal
samtools sort SRR1691646_Normal.bam SRR1691646_Normal
samtools sort SRR1691647_Normal.bam SRR1691647_Normal
samtools sort SRR1691648_Normal.bam SRR1691648_Normal
# Index the corrected bam
samtools index SRR1691633_ZRSR2Mut.bam
samtools index SRR1691634_ZRSR2Mut.bam
samtools index SRR1691635_ZRSR2Mut.bam
samtools index SRR1691636_ZRSR2Mut.bam
samtools index SRR1691637_ZRSR2Mut.bam
samtools index SRR1691638_ZRSR2Mut.bam
samtools index SRR1691639_ZRSR2Mut.bam
samtools index SRR1691640_ZRSR2Mut.bam
samtools index SRR1691641_WT.bam
samtools index SRR1691642_WT.bam
samtools index SRR1691643_WT.bam
samtools index SRR1691644_WT.bam
samtools index SRR1691645_Normal.bam
samtools index SRR1691646_Normal.bam
samtools index SRR1691647_Normal.bam
samtools index SRR1691648_Normal.bam

#Remove the sam files
rm *.sam

#make a tmp directory (folder)
mkdir tmp
# Limit the data even further to the U12 genes regions
samtools view -b -L u12Regions.bed -o tmp/SRR1691633_ZRSR2Mut.bam SRR1691633_ZRSR2Mut.bam
samtools view -b -L u12Regions.bed -o tmp/SRR1691634_ZRSR2Mut.bam SRR1691634_ZRSR2Mut.bam
samtools view -b -L u12Regions.bed -o tmp/SRR1691635_ZRSR2Mut.bam SRR1691635_ZRSR2Mut.bam
samtools view -b -L u12Regions.bed -o tmp/SRR1691636_ZRSR2Mut.bam SRR1691636_ZRSR2Mut.bam
samtools view -b -L u12Regions.bed -o tmp/SRR1691637_ZRSR2Mut.bam SRR1691637_ZRSR2Mut.bam
samtools view -b -L u12Regions.bed -o tmp/SRR1691638_ZRSR2Mut.bam SRR1691638_ZRSR2Mut.bam
samtools view -b -L u12Regions.bed -o tmp/SRR1691639_ZRSR2Mut.bam SRR1691639_ZRSR2Mut.bam
samtools view -b -L u12Regions.bed -o tmp/SRR1691640_ZRSR2Mut.bam SRR1691640_ZRSR2Mut.bam
samtools view -b -L u12Regions.bed -o tmp/SRR1691641_WT.bam SRR1691641_WT.bam
samtools view -b -L u12Regions.bed -o tmp/SRR1691642_WT.bam SRR1691642_WT.bam
samtools view -b -L u12Regions.bed -o tmp/SRR1691643_WT.bam SRR1691643_WT.bam
samtools view -b -L u12Regions.bed -o tmp/SRR1691644_WT.bam SRR1691644_WT.bam
samtools view -b -L u12Regions.bed -o tmp/SRR1691645_Normal.bam SRR1691645_Normal.bam
samtools view -b -L u12Regions.bed -o tmp/SRR1691646_Normal.bam SRR1691646_Normal.bam
samtools view -b -L u12Regions.bed -o tmp/SRR1691647_Normal.bam SRR1691647_Normal.bam
samtools view -b -L u12Regions.bed -o tmp/SRR1691648_Normal.bam SRR1691648_Normal.bam

samtools index tmp/SRR1691633_ZRSR2Mut.bam
samtools index tmp/SRR1691634_ZRSR2Mut.bam
samtools index tmp/SRR1691635_ZRSR2Mut.bam
samtools index tmp/SRR1691636_ZRSR2Mut.bam
samtools index tmp/SRR1691637_ZRSR2Mut.bam
samtools index tmp/SRR1691638_ZRSR2Mut.bam
samtools index tmp/SRR1691639_ZRSR2Mut.bam
samtools index tmp/SRR1691640_ZRSR2Mut.bam
samtools index tmp/SRR1691641_WT.bam
samtools index tmp/SRR1691642_WT.bam
samtools index tmp/SRR1691643_WT.bam
samtools index tmp/SRR1691644_WT.bam
samtools index tmp/SRR1691645_Normal.bam
samtools index tmp/SRR1691646_Normal.bam
samtools index tmp/SRR1691647_Normal.bam
samtools index tmp/SRR1691648_Normal.bam

# Remove all previously generatyed files
rm *.bam*
# move everything from tmp directory to current path
mv tmp/*.* ./
#Remove tmp directory
rm -r tmp/

