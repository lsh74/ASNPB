import pandas as pd

# 读取txt文件
SE_data = pd.read_csv("/data3/lishuhan/30GSE_health_SRR/v3AS/old_young/output_rmats/SE.MATS.JC.txt", sep='\t')
MXE_data = pd.read_csv("/data3/lishuhan/30GSE_health_SRR/v3AS/old_young/output_rmats/MXE.MATS.JC.txt", sep='\t')
RI_data = pd.read_csv("/data3/lishuhan/30GSE_health_SRR/v3AS/old_young/output_rmats/RI.MATS.JC.txt", sep='\t')
A5SS_data = pd.read_csv("/data3/lishuhan/30GSE_health_SRR/v3AS/old_young/output_rmats/A5SS.MATS.JC.txt", sep='\t')
A3SS_data = pd.read_csv("/data3/lishuhan/30GSE_health_SRR/v3AS/old_young/output_rmats/A3SS.MATS.JC.txt", sep='\t')
SE_data['splicing_type'] = 'SE'
MXE_data['splicing_type'] = 'MXE'
A5SS_data['splicing_type'] = 'A5SS'
A3SS_data['splicing_type'] = 'A3SS'
RI_data['splicing_type'] = 'RI'

# 查看数据的前几行，默认为前5行
print(SE_data.head())
print("SE列名：", SE_data.columns)
print("MXE列名：", MXE_data.columns)
print("RI列名：", RI_data.columns)
print("A5SS列名：", A5SS_data.columns)
print("A3SS列名：", A3SS_data.columns)

# 最终合并的数据集的所有列名
col_names = ['ID', 'splicing_type', 'GeneID', 'geneSymbol', 'chr', 'strand', 'PValue', 'FDR', 'IncLevelDifference',
             'exonStart_0base','exonEnd', 'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE',   # SE
             'longExonStart_0base', 'longExonEnd', 'shortES', 'shortEE', 'flankingES', 'flankingEE',    # A5SS和A3SS
             '1stExonStart_0base', '1stExonEnd', '2ndExonStart_0base', '2ndExonEnd',                    # MXE
             'riExonStart_0base', 'riExonEnd'                                                           # RI
             ]

type(SE_data)

# 合并数据集
rMAts_results = pd.concat([SE_data, MXE_data, RI_data, A5SS_data, A3SS_data], ignore_index=True)

# 处理缺失的列
for col_name in col_names:
    if col_name not in rMAts_results.columns:
        rMAts_results[col_name] = '.'


rMAts_results = rMAts_results[col_names]

# 保存结果到文件
rMAts_results.to_csv('/data3/lishuhan/30GSE_health_SRR/v3AS/old_young/files/rMAts_results.csv', index=False)
