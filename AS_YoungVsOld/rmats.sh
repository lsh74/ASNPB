cd /data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/AS/file/
/data3/lishuhan/copy/lishuhan/mini/bin/python /data3/lishuhan/copy/lishuhan/mini/rMATS/rmats1.py \
--b1 ./b1.txt \
--b2 ./b2.txt \
--gtf /nas01/Genome/fasta_gtf/Homo_sapiens.GRCh38.84.gtf \
--od /data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/AS/file/output_rmats \
-t paired \
--nthread 25 \
--readLength 150 \
--variable-read-length \
--tmp /data3/lishuhan/30GSE_health_SRR/v3AS/v2_old_you/repair/AS/file/out_rmats \
--drop-zero-read-replicates-for-stat


