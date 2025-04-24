#!/bin/bash

work_dir="/data3/lishuhan/30GSE_health_SRR"
output_fasta="/data3/lishuhan/30GSE_health_SRR/mid_peptide/combined_peptides.fa"
temp_fasta="/data3/lishuhan/30GSE_health_SRR/mid_peptide/temp_peptides.fa"

# ��ջ򴴽�������ļ�
> "$output_fasta"
> "$temp_fasta"

# ��ȡ CSV �ļ�
tail -n +2 /data3/lishuhan/30GSE_health_SRR/v2AS/v2SRR/v2_df_17_59.csv | while IFS=',' read -r run sample_name status age sex dataset
do
    run=$(echo $run | tr -d '"')
    dataset=$(echo $dataset | tr -d '"')
    input_pep="$work_dir/$dataset/tabfile/pep/${run}_pep.fa"

    if [ -f "$input_pep" ]; then
        echo "Processing $input_pep..."

        # ��� output_fasta Ϊ�գ�ֱ�Ӹ��� input_pep �� output_fasta
        if [ ! -s "$output_fasta" ]; then
            cp "$input_pep" "$output_fasta"
        else
            # �ϲ� input_pep �� output_fasta �� temp_fasta
            awk '
            BEGIN {
                RS = ">"
                ORS = ""
                FS = "\n"
                OFS = "\n"
            }

            # ��ȡ��һ���ļ�������������ӳ��
            FNR == NR && NR > 1 {
                id = $1
                peptides = $2
                split(peptides, pep_array, ",")
                for (i in pep_array) {
                    if (pep_array[i] != "") {
                        peptides_map[id][pep_array[i]] = 1
                    }
                }
                next
            }

            # ��ȡ�ڶ����ļ�������������ӳ��
            NR > FNR {
                id = $1
                peptides = $2
                split(peptides, pep_array, ",")
                for (i in pep_array) {
                    if (pep_array[i] != "") {
                        peptides_map[id][pep_array[i]] = 1
                    }
                }
            }

            END {
                for (id in peptides_map) {
                    output_pep = ""
                    for (pep in peptides_map[id]) {
                        output_pep = (length(output_pep) > 0 ? output_pep "," : "") pep
                    }
                    print ">" id "\n" output_pep "\n"
                }
            }' "$output_fasta" "$input_pep" > "$temp_fasta"
            
            # ���� output_fasta
            mv "$temp_fasta" "$output_fasta"
            > "$temp_fasta"
        fi
    else
        echo "Peptide file not found for $run at $input_pep"
    fi
done

echo "All peptide sequences have been processed and combined into $output_fasta."
