featureCounts \
-T 8 \
-p \
-a /group/sbms006/rrubimsilvadocouto/reference_genome/mouse_GRCm39/gencode.vM32.annotation.gtf \
-o /group/sbms006/rrubimsilvadocouto/MiSeq_Data/gene_counts.txt \
$(find /group/sbms006/rrubimsilvadocouto/MiSeq_Data/STAR_output \
-name "Aligned.sortedByCoord.out.bam")
