WHITELIST=/home/data/737K-cratac-v1.txt
WHITELIST=/media/data/AtteR/scifi-analysis/737K-cratac-v1.txt

    
solo_exec () {
DIR=${EXP}_${ID}             
STAR \
--genomeDir $MAP \
--runThreadN 24 \
--soloType CB_UMI_Simple \
--soloCBwhitelist $WHITELIST \
--soloCBstart 1 \
--soloCBlen 16 \
--soloUMIstart 17 \
--soloUMIlen 8 \
--soloBarcodeReadLength 0 \
--soloFeatures GeneFull \
--readFilesCommand zcat \
--soloMultiMappers Uniform EM PropUnique Unique \
--soloUMIdedup 1MM_All \
--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
--genomeLoad LoadAndKeep \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--readFilesIn $RV $FW 
mv Solo.out $DIR
mv Log.final.out Log.out Log.progress.out SJ.out.tab Aligned.sortedByCoord.out.bam Aligned.out.sam $DIR/
chmod -R 755 $DIR
rm -rf _STARtmp
# --soloCellFilter TopCells 1 0.99 10 \
# --soloOutFileNames $DIR/ features.tsv barcodes.tsv matrix.mtx \
# --outSAMtype BAM SortedByCoordinate \
# --outSAMattributes UR CB UB GX GN \
# --limitBAMsortRAM 220000000000 \
}

DATAPTH=/home/data/output_cutadapt_novaseq_r2
DATAPTH=/media/data/AtteR/scifi-analysis/pilot1_output_trim
DATAOUT=/media/data/AtteR/scifi-analysis/outputs_starsolo

DATAPTH=/media/data/AtteR/scifi-analysis/scifi6/trimmed
DATAOUT=/media/data/AtteR/scifi-analysis/scifi6/outputs

mkdir $DATAOUT
cd $DATAOUT


MAP=/media/data/AtteR/scifi-analysis/ref_genome/rat_index

# ####

EXP=Scifi_library
EXP=Scifi_library_2_S2
EXP=SciFi_
ID=SciFi5_oDT
#we take the files with the above expressions, then run the solo_exec function
FW=$DATAPTH/${EXP}_trimmed-SciFi5_oDT_R21.fastq.gz
RV=$DATAPTH/${EXP}_trimmed-SciFi5_oDT_R3.fastq.gz


#solo_exec 

ID=SciFi6_oDT

FW=$DATAPTH/${EXP}_trimmed-SciFi6_oDT_R21.fastq.gz
RV=$DATAPTH/${EXP}_trimmed-SciFi6_oDT_R3.fastq.gz


#solo_exec

ID=SciFi7_oDT

FW=$DATAPTH/${EXP}_trimmed-SciFi7_oDT_R21.fastq.gz
RV=$DATAPTH/${EXP}_trimmed-SciFi7_oDT_R3.fastq.gz

#solo_exec

ID=SciFi8_oDT

FW=$DATAPTH/${EXP}_trimmed-SciFi8_oDT_R21.fastq.gz
RV=$DATAPTH/${EXP}_trimmed-SciFi8_oDT_R3.fastq.gz

#solo_exec



echo "STARTING..."

# #####
ID=SciFi5_WP

FW=$DATAPTH/${EXP}_trimmed-SciFi5_WP_R21.fastq.gz
RV=$DATAPTH/${EXP}_trimmed-SciFi5_WP_R3.fastq.gz
solo_exec 

ID=SciFi6_WP

FW=$DATAPTH/${EXP}_trimmed-SciFi6_WP_R21.fastq.gz
RV=$DATAPTH/${EXP}_trimmed-SciFi6_WP_R3.fastq.gz
solo_exec

ID=SciFi7_WP

FW=$DATAPTH/${EXP}_trimmed-SciFi7_WP_R21.fastq.gz
RV=$DATAPTH/${EXP}_trimmed-SciFi7_WP_R3.fastq.gz
solo_exec

ID=SciFi8_WP

FW=$DATAPTH/${EXP}_trimmed-SciFi8_WP_R21.fastq.gz
RV=$DATAPTH/${EXP}_trimmed-SciFi8_WP_R3.fastq.gz
solo_exec

# ####


STAR --genomeDir $MAP --genomeLoad Remove
