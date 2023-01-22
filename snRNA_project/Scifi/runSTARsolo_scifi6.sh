WHITELIST=/home/data/737K-cratac-v1.txt
WHITELIST=/media/data/AtteR/scifi-analysis/737K-cratac-v1_RC.txt

    
solo_exec () {
DIR=${EXP}_${ID}             
/media/data/AtteR/Attes_bin/STAR-2.7.10a/bin/Linux_x86_64/STAR \
--genomeDir $MAP \
--runThreadN 24 \
--soloType CB_UMI_Simple \
--soloCBwhitelist $WHITELIST \
--soloCBstart 1 \
--soloCBlen 16 \
--soloUMIstart 17 \
--soloUMIlen 12 \  #earlier was 8
--soloBarcodeMate 0 --soloBarcodeReadLength 0 \
--soloFeatures GeneFull \
--readFilesCommand zcat \
--soloMultiMappers Uniform EM PropUnique Unique \
--soloUMIdedup 1MM_All \
--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
--genomeLoad LoadAndKeep \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--readFilesIn $RV, $FW 
mkdir $DIR
mv Solo.out $DIR
mv Log.final.out Log.out Log.progress.out SJ.out.tab Aligned.sortedByCoord.out.bam Aligned.out.sam $DIR/
chmod -R 755 $DIR
rm -rf _STARtmp
--soloCellFilter TopCells 1 0.99 10 \
# --soloOutFileNames $DIR/ features.tsv barcodes.tsv matrix.mtx \
# --outSAMtype BAM SortedByCoordinate \
# --outSAMattributes UR CB UB GX GN \
# --limitBAMsortRAM 220000000000 \
}

#DATAPTH=/media/data/AtteR/scifi-analysis
DATAPTH=/media/data/AtteR/scifi-analysis/scifi6/2nd_try/output
DATAOUT=/media/data/AtteR/scifi-analysis/scifi6/2nd_try/starsolo_outputs2

#DATAPTH=/media/data/AtteR/scifi-analysis/scifi6/output
#DATAOUT=/media/data/AtteR/scifi-analysis/scifi6/scifi_output

mkdir $DATAOUT
cd $DATAOUT


MAP=/media/data/AtteR/scifi-analysis/ref_genome/rat_index_aav

# #### A9, G8, H8 F8_Asyn, 5, 6

#POOL1
EXP="SciFi_new_cDNA_Pool1__trimmed-Pool1"
ID="_oDT_A9_"
#we take the files with the above expressions, then run the solo_exec function

#A9
FW=$DATAPTH/${EXP}${ID}R21.fastq.gz

echo "FW strand: $FW"
RV=$DATAPTH/${EXP}${ID}R3.fastq.gz
solo_exec 

ID="_WP_A9_"
FW=$DATAPTH/${EXP}${ID}R21.fastq.gz
echo $FW
RV=$DATAPTH/${EXP}${ID}R3.fastq.gz
solo_exec 

echo "Done A9, pool1"
#F8_Asyn
ID="_oDT_F8_Asyn_"
FW=$DATAPTH/${EXP}${ID}R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}R3.fastq.gz
solo_exec 

ID="_WP_F8_Asyn_"
FW=$DATAPTH/${EXP}${ID}R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}R3.fastq.gz
solo_exec 

echo "Done F8_Asyn, pool1"

#G8
ID="_oDT_G8_"
FW=$DATAPTH/${EXP}${ID}R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}R3.fastq.gz
solo_exec 
ID="_WP_G8_"
FW=$DATAPTH/${EXP}${ID}R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}R3.fastq.gz
solo_exec 

#H8
ID="_oDT_H8_"
FW=$DATAPTH/${EXP}${ID}R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}R3.fastq.gz
solo_exec 

ID="_WP_H8_"
FW=$DATAPTH/${EXP}${ID}R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}R3.fastq.gz
solo_exec 

# #5
# ID="_oDT_5_"
# FW=$DATAPTH/${EXP}${ID}R21.fastq.gz
# RV=$DATAPTH/${EXP}${ID}R3.fastq.gz
# solo_exec 

# ID="_WP_5_"
# FW=$DATAPTH/${EXP}${ID}R21.fastq.gz
# RV=$DATAPTH/${EXP}${ID}R3.fastq.gz
# solo_exec 

# #6
# ID="_oDT_6_"
# FW=$DATAPTH/${EXP}${ID}R21.fastq.gz
# RV=$DATAPTH/${EXP}${ID}R3.fastq.gz
# solo_exec 

# ID="_WP_6_"
# FW=$DATAPTH/${EXP}${ID}R21.fastq.gz
# RV=$DATAPTH/${EXP}${ID}R3.fastq.gz
# solo_exec 


echo "Pool1 done!"

#Pool2
EXP="SciFi_old_cDNA_Pool2__trimmed-Pool2"
ID="_oDT_A9_"
#we take the files with the above expressions, then run the solo_exec function

# #A9
# FW=$DATAPTH/${EXP}${ID}R21.fastq.gz
# echo $FW
# RV=$DATAPTH/${EXP}${ID}R3.fastq.gz
# solo_exec 

# ID="_WP_A9_"
# FW=$DATAPTH/${EXP}${ID}R21.fastq.gz
# echo $FW
# RV=$DATAPTH/${EXP}${ID}R3.fastq.gz
# solo_exec 

# echo "Done A9, pool1"
# #F8_Asyn
# ID="_oDT_F8_Asyn_"
# FW=$DATAPTH/${EXP}${ID}R21.fastq.gz
# RV=$DATAPTH/${EXP}${ID}R3.fastq.gz
# solo_exec 

# ID="_WP_F8_Asyn_"
# FW=$DATAPTH/${EXP}${ID}R21.fastq.gz
# RV=$DATAPTH/${EXP}${ID}R3.fastq.gz
# solo_exec 

# echo "Done F8_Asyn, pool1"

# #G8
# ID="_oDT_G8_"
# FW=$DATAPTH/${EXP}${ID}R21.fastq.gz
# RV=$DATAPTH/${EXP}${ID}R3.fastq.gz
# solo_exec 
# ID="_WP_G8_"
# FW=$DATAPTH/${EXP}${ID}R21.fastq.gz
# RV=$DATAPTH/${EXP}${ID}R3.fastq.gz
# solo_exec 

# #H8
# ID="_oDT_H8_"
# FW=$DATAPTH/${EXP}${ID}R21.fastq.gz
# RV=$DATAPTH/${EXP}${ID}R3.fastq.gz
# solo_exec 

# ID="_WP_H8_"
# FW=$DATAPTH/${EXP}${ID}R21.fastq.gz
# RV=$DATAPTH/${EXP}${ID}R3.fastq.gz
# solo_exec 

#5
ID="_oDT_5_"
FW=$DATAPTH/${EXP}${ID}R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}R3.fastq.gz
solo_exec 

ID="_WP_5_"
FW=$DATAPTH/${EXP}${ID}R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}R3.fastq.gz
solo_exec 

#6
ID="_oDT_6_"
FW=$DATAPTH/${EXP}${ID}R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}R3.fastq.gz
solo_exec 

ID="_WP_6_"
FW=$DATAPTH/${EXP}${ID}R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}R3.fastq.gz
solo_exec 


echo "Pool2 done!"




STAR --genomeDir $MAP --genomeLoad Remove
