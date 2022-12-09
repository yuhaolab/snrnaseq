
export PATH=/home/yh/2-nano/backup/scrnaseq/cellranger/yard/apps/cellranger-6.1.2:$PATH

## make reference
cd /home/yh/2-nano/backup/scrnaseq/ref

cellranger mkgtf \
Oryza_sativa.IRGSP-1.0.52.gtf \
Oryza_sativa.IRGSP-1.0.52.filtered.gtf \
--attribute=gene_biotype:protein_coding \
--attribute=gene_biotype:antisense_RNA \
--attribute=gene_biotype:nontranslating_CDS \
--attribute=gene_biotype:pre_miRNA \
--attribute=gene_biotype:rRNA \
--attribute=gene_biotype:sense_intronic \
--attribute=gene_biotype:snoRNA \
--attribute=gene_biotype:snRNA

cellranger mkref \
--genome=Oryza_sativa.IRGSP_1.0 \
--fasta=Oryza_sativa.IRGSP-1.0.dna.toplevel.fa \
--genes=Oryza_sativa.IRGSP-1.0.52.filtered.gtf \
--nthreads=30 --memgb=128


## reads count
cd /home/yh/2-nano/backup/scrnaseq/result
cellranger count --id=rice_test1 \
                 --transcriptome=/home/yh/2-nano/backup/scrnaseq/ref/Oryza_sativa.IRGSP_1.0 \
                 --fastqs=/home/yh/2-nano/backup/scrnaseq/X401SC21115819-Z01-F001/raw_data/Rice \
                 --sample=Rice-SCI7T034-AK33418_H3C7CDSX3 \
                 --expect-cells=10000 \
                 --include-introns \
                 --localcores=16 \
                 --localmem=128