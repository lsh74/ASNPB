from Bio.Seq import Seq
from Bio import SeqIO
from functools import reduce
import os
from multiprocessing.dummy import Pool as MyPool
from time import sleep, ctime

# 1. 提取序列信息hg38.fa   ----------------------------------------------------
hg38_dict = {}
for seq_record in SeqIO.parse("/home/Jupyter/Python/Alternative_splicing/data/hg38.fa", "fasta"):
    hg38_dict[seq_record.id.strip('chr')] = seq_record
# 判断染色体5是否为纯N
chromosome_5_sequence = hg38_dict['5'].seq
is_all_N = all(base == 'N' for base in chromosome_5_sequence)
print(is_all_N)
# 查看其中一个序列
hg38_dict['5'].seq
chromosome_5_sequence = hg38_dict['5'].seq
# 统计非 'N' 的数量
non_N_count = sum(base != 'N' for base in chromosome_5_sequence)
# 计算非 'N' 的比例
non_N_ratio = non_N_count / len(chromosome_5_sequence)
print("非 'N' 的数量:", non_N_count)
print("非 'N' 的比例:", non_N_ratio)


# 2. 提取基因组注释信息hg38.gtf ---------------------------------------------------------------------
h38gtf_file_path = '/home/Jupyter/Python/Alternative_splicing/data/Homo_sapiens.GRCh38.84.gtf'
# 查看前几行
with open(h38gtf_file_path, 'r') as file:
    for i, line in enumerate(file):
        if i < 2:
            columns = line.strip().split('\t')
            print(f"Line {i + 1}: {columns}, Number of elements: {len(columns)}")
        else:
            break

from collections import Counter
# 统计gtf文件中序列（feature）类型的频数
type_counter = Counter()
# 遍历GTF文件中的每一行
for line in open(h38gtf_file_path, 'r'):
    # 提取第三列的元素
    type_element = line.strip().split('\t')[2]
    type_counter[type_element] += 1
# 打印各种元素出现的频数
for element, count in type_counter.items():
    print(f"{element}: {count}")

# 1）Gene_dict（每个基因的所有转录本及对应CDS位置） ------------------------------------
import re

chr_list = [str(each) for each in range(1, 23)] + ['X', 'Y']
type_list = ['CDS', 'stop_codon', 'start_codon']
Gene_dict = {}

mark = ''
# 遍历GTF文件中的每一行
for line in open(h38gtf_file_path, 'r'):
    # 提取第九列的属性信息
    attributes = line.strip().split('\t')[8].strip(';')
    # 使用正则表达式通过空格或分号拆分属性，忽略引号内的空格
    fields = dict(re.findall(r'(\S+)\s+"([^"]+)"', attributes))

    # 从行中提取相关信息
    gene_id = fields['gene_id']
    transcript_id = fields.get('transcript_id', '')
    exon_number = fields.get('exon_number', '')
    strand = line.strip().split('\t')[6]
    chr = line.strip().split('\t')[0].strip('chr')
    type = line.strip().split('\t')[2]
    start = int(line.strip().split('\t')[3])
    end = int(line.strip().split('\t')[4])

    # 跳过非基因条目
    if 'gene_name' not in fields:
        continue

    if (chr not in chr_list) or (type not in type_list):
        continue

    if (type == "start_codon"):
        type = "startcodon"
    if (type == "stop_codon"):
        type = "stopcodon"

    if transcript_id != mark:
        mark = transcript_id
        last_n = 0

    gene_name = fields['gene_name']

    #  Gene_dict.setdefault(gene_name, {})
    # 创建一个新的键值对，键是 gene_name，对应的值是一个空字典 {}。如果前一行的 gene_name 和当前行相同，它不会删除前一行记录的信息，
    # 而是在原有的基因条目上继续添加新的信息。这样可以确保一个基因名下包含了所有相关的转录本和外显子信息。
    Gene_dict.setdefault(gene_name, {})
    if transcript_id:
        Gene_dict[gene_name].setdefault(transcript_id, {'CDS': [strand]})  # 同上 Gene_dict.setdefault(gene_name, {})
        if exon_number:
            CDSID = '_'.join(
                [str(chr), str(type), str(start), str(end), str(strand), str(exon_number)])  # +exon_number，
            Gene_dict[gene_name][transcript_id]['CDS'].append(CDSID)
            Gene_dict[gene_name][transcript_id][CDSID] = [chr, max(int(start), last_n + 1), end]
            last_n = int(end)

# 输出 Gene_dict 的前几个键值对
for gene_name, transcripts in list(Gene_dict.items())[:2]:
    print(f"Gene: {gene_name}")
    for transcript_id, info in transcripts.items():
        print(f"  Transcript: {transcript_id}")
        print(f"    CDS: {info['CDS']}")
        for cds_id, cds_info in info.items():
            if cds_id != 'CDS':
                print(f"    {cds_id}: {cds_info}")
    print("\n")

# 2）Max_Gene_ditc 遍历 Gene_dict 中的每个基因及其转录本，找到每个基因编码蛋白质最长的转录本，并记录其相关信息到 Max_Gene_ditc 字典中
Max_Gene_ditc = {}
# 遍历 Gene_dict 中的每个基因
for eachgene in Gene_dict:
    # 初始化变量 exon_len 记录最长的编码序列长度，start 记录最小的起始位置
    exon_len = 0
    start = 9999999999999999999999999999999999999999999999
    # 遍历当前基因的每个转录本
    for eachtran in Gene_dict[eachgene]:
        # 获取转录本的编码方向和当前转录本的所有CDS信息
        stand = Gene_dict[eachgene][eachtran]['CDS'][0]
        AllCDS = Gene_dict[eachgene][eachtran]['CDS'][1:]
        # 计算当前转录本的编码序列长度
        LenCDS = sum([int(each2.split('_')[3]) - int(each2.split('_')[2]) for each2 in AllCDS])

        # 获取编码序列的起始位置
        if stand == '+':
            startexon = Gene_dict[eachgene][eachtran]['CDS'][1]
            startpos = Gene_dict[eachgene][eachtran][startexon][1]
        else:
            startexon = Gene_dict[eachgene][eachtran]['CDS'][-1]
            startpos = Gene_dict[eachgene][eachtran][startexon][2]

        # 如果当前转录本的起始位置比已记录的最小位置还小，则更新 Max_Gene_ditc 中当前基因的记录为当前转录本
        if int(stand + str(startpos)) < start:
            start = int(stand + str(startpos))
            Max_Gene_ditc[eachgene] = Gene_dict[eachgene][eachtran]
            Max_Gene_ditc[eachgene]['Rid'] = eachtran
            exon_len = LenCDS
        # 如果当前转录本的起始位置与已记录的最小位置相同，但编码序列长度比已记录的最长长度还长，则更新 Max_Gene_ditc 中当前基因的记录为当前转录本
        elif int(stand + str(startpos)) == start and LenCDS > exon_len:
            exon_len = LenCDS
            Max_Gene_ditc[eachgene] = Gene_dict[eachgene][eachtran]
            Max_Gene_ditc[eachgene]['Rid'] = eachtran

# 输出 Max_Gene_ditc 的前几个键值对
for gene_name, transcripts in list(Gene_dict.items())[:2]:
    print(f"Gene: {gene_name}")
    for transcript_id, info in transcripts.items():
        print(f"  Transcript: {transcript_id}")
        for cds_id, cds_info in info.items():
            if cds_id != 'CDS':
                print(f"    {cds_id}: {cds_info}")
    print("\n")

#Gene_dict.get("SAMD11")
#Max_Gene_ditc.get("SAMD11")
# 代码测试
stand = Gene_dict["SAMD11"]['ENST00000622503']['CDS'][0]
AllCDS = Gene_dict["SAMD11"]['ENST00000622503']['CDS'][1:]
# CDSID = '_'.join([chr, exon_number, start, end, strand])
LenCDS = sum([int(each2.split('_')[3]) - int(each2.split('_')[2]) for each2 in AllCDS])
stand = Gene_dict["SAMD11"]['ENST00000341065']['CDS'][0]
AllCDS = Gene_dict["SAMD11"]['ENST00000341065']['CDS'][1:]
# CDSID = '_'.join([chr, exon_number, start, end, strand])
LenCDS2 = sum([int(each2.split('_')[3]) - int(each2.split('_')[2]) for each2 in AllCDS])
print(LenCDS)
LenCDS2


# 3）Tran_pos_dict（每个基因的所有转录本的位置信息）----------------------------------------------------------------------------
##Get thr pos of each tran 提取每个转录本的位置
type_list = ['exon']
Tran_pos_dict = {}
# Tran_pos_dict
for line in open(h38gtf_file_path, 'r'):
    # 提取第九列的属性信息
    attributes = line.strip().split('\t')[8].strip(';')
    # 使用正则表达式通过空格或分号拆分属性，忽略引号内的空格
    fields = dict(re.findall(r'(\S+)\s+"([^"]+)"', attributes))

    # 从行中提取相关信息
    gene_id = fields['gene_id']
    transcript_id = fields.get('transcript_id', '')
    exon_number = fields.get('exon_number', '')
    strand = line.strip().split('\t')[6]
    chr = line.strip().split('\t')[0].strip('chr')
    type = line.strip().split('\t')[2]
    start = int(line.strip().split('\t')[3])
    end = int(line.strip().split('\t')[4])

    # 跳过非基因条目
    if 'gene_name' not in fields:
        continue

    if (chr not in chr_list) or (type not in type_list) or (transcript_id.startswith('NR')):
        continue
    gene_name = fields['gene_name']
    Tran_pos_dict.setdefault(gene_name, {})
    Tran_pos_dict[gene_name].setdefault(transcript_id, {})
    Tran_pos_dict[gene_name][transcript_id].setdefault('start', []).append(start)
    Tran_pos_dict[gene_name][transcript_id].setdefault('end', []).append(end)

# 输出 Tran_pos_dict 的前几个键值对
for gene_name, transcripts in list(Tran_pos_dict.items())[:2]:
    print(f"Gene: {gene_name}")
    for transcript_id, info in transcripts.items():
        print(f"  Transcript: {transcript_id}")
        for cds_id, cds_info in info.items():
            print(f"    {cds_id}: {cds_info}")
    print("\n")

# 3）Get the exon pos 提取外显子的位置信息  ----------------------------------------------------------------------------
Gene_pos = {}

# 遍历 GTF 文件的每一行
for line in open(h38gtf_file_path, 'r'):
    # 解析每一行的字段
    fields = line.strip().split('\t')
    chr = fields[0].strip('chr')
    type = fields[2]
    start = int(fields[3])
    end = int(fields[4])
    strand = fields[6]

    # 仅处理基因和外显子类型的行
    if type not in ['gene', 'exon']:
        continue

    # 提取基因名和外显子编号
    attributes = fields[8]
    gene_id = re.search(r'gene_id "(.*?)";', attributes).group(1)
    exon_number_match = re.search(r'exon_number "(.*?)";', attributes)
    exon_number = exon_number_match.group(1) if exon_number_match else None

    # 如果是基因行，则创建基因的字典条目
    if type == 'gene':
        Gene_pos.setdefault(gene_id, {})
    # 如果是外显子行，则将外显子的位置信息存储到基因的字典中
    elif type == 'exon':
        if gene_id in Gene_pos:
            Gene_pos[gene_id].setdefault(exon_number, [])
            if strand == '+':
                Gene_pos[gene_id][exon_number] = [start, end]
            else:
                Gene_pos[gene_id][exon_number] = [end, start]

# 输出 Gene_pos 的前几个键值对
for gene_name, exons in list(Gene_pos.items())[:2]:
    print(f"Gene: {gene_name}")
    for exon_number, pos in exons.items():
        print(f"  Exon {exon_number}: {pos}")
    print("\n")

