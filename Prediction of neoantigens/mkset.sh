#!/bin/bash

work_dir="/data3/lishuhan/30GSE_health_SRR"
output_fasta="/data3/lishuhan/30GSE_health_SRR/mid_peptide/combined_peptides.fa"
temp_fasta="/data3/lishuhan/30GSE_health_SRR/mid_peptide/temp_peptides.fa"

# 清空或创建总输出文件
> "$output_fasta"
> "$temp_fasta"

# 读取 CSV 文件
tail -n +2 /data3/lishuhan/30GSE_health_SRR/v2AS/v2SRR/v2_df_17_59.csv | while IFS=',' read -r run sample_name status age sex dataset
do
    run=$(echo $run | tr -d '"')
    dataset=$(echo $dataset | tr -d '"')
    input_pep="$work_dir/$dataset/tabfile/pep/${run}_pep.fa"

    if [ -f "$input_pep" ]; then
        echo "Processing $input_pep..."

        # 如果 output_fasta 为空，直接复制 input_pep 到 output_fasta
        if [ ! -s "$output_fasta" ]; then
            cp "$input_pep" "$output_fasta"
        else
            # 合并 input_pep 和 output_fasta 到 temp_fasta
            awk '
            BEGIN {
                RS = ">"
                ORS = ""
                FS = "\n"
                OFS = "\n"
            }

            # 读取第一个文件并构建肽序列映射
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

            # 读取第二个文件并更新肽序列映射
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
            
            # 更新 output_fasta
            mv "$temp_fasta" "$output_fasta"
            > "$temp_fasta"
        fi
    else
        echo "Peptide file not found for $run at $input_pep"
    fi
done

echo "All peptide sequences have been processed and combined into $output_fasta."
