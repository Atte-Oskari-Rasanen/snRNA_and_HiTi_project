#!/bin/bash
#Rat
#split-pipe --mode all --kit WT_mini --genome_dir /Parse/genomes/rat72_AAV_mix/ \
#    --fq1 /Parse/expdata/ParseBio_WT-mini_Sample1_S1_R1_001.fastq.gz \
#    --fq2 /Parse/expdata/ParseBio_WT-mini_Sample1_S1_R2_001.fastq.gz \
#    --output_dir /Parse/analysis_humrat_s1/ \
#    --sample rat_bfp_sn_gfp A1 \
#    --sample rat_asyn1_sn_nogfp A2 \
#    --sample rat_asyn3_sn_gfp A3 \
#    --sample rat_asyn4_sn_nogfp A4 \
#    --sample mus_asyn1_sn_gfp A5 \
#    --sample mus_asyn2_sn_nogfp A6 \
#    --sample mus_asyn3_sn_nogfp A7 \
#    --sample humrat_asyn1_crtx_gfp A9 \
#    --sample humrat_asyn1_sn_nogfp A10 \
#    --sample humrat_asyn2_sn_nogfp A11 \
#    --sample humrat_asyn3_sn_nogfp A12

split-pipe --mode all --kit WT_mini --genome_dir /Parse/genomes/rat72_AAV_mix/ \
    --fq1 /Parse/expdata/ParseBio_WT-mini_Sample2_S2_R1_001.fastq.gz \
    --fq2 /Parse/expdata/ParseBio_WT-mini_Sample2_S2_R2_001.fastq.gz \
    --output_dir /Parse/analysis_humrat_2/ \
    --sample rat_bfp_sn_gfp A1 \
    --sample rat_asyn1_sn_nogfp A2 \
    --sample rat_asyn3_sn_gfp A3 \
    --sample rat_asyn4_sn_nogfp A4 \
    --sample mus_asyn1_sn_gfp A5 \
    --sample mus_asyn2_sn_nogfp A6 \
    --sample mus_asyn3_sn_nogfp A7 \
    --sample humrat_asyn1_crtx_gfp A9 \
    --sample humrat_asyn1_sn_nogfp A10 \
    --sample humrat_asyn2_sn_nogfp A11 \
    --sample humrat_asyn3_sn_nogfp A12
