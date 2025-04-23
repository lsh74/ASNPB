# -*- coding: utf-8 -*-
import os
import wget
import subprocess

# 读取包含GSM号码的文件
with open('gds_result.txt', 'r') as f:
    lines = f.readlines()

# 提取所有的GSM号码
gsm_numbers = [line.split(': ')[1].split()[0] for line in lines if line.startswith('Sample')]

# 创建一个新的文件夹来存储下载的文件
os.makedirs('downloaded_files', exist_ok=True)

# 获取GSM号码的总数，以便计算进度
total_gsm = len(gsm_numbers)

# 对每个GSM号码，下载对应的txt文件
for i, gsm_number in enumerate(gsm_numbers, start=1):
    print(f'Downloading {gsm_number} ({i}/{total_gsm})...')
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm_number}&targ=self&form=text&view=quick"
    output_file = f'downloaded_files/{gsm_number}.txt'
    
    # 使用wget下载文件，并使用-c选项实现断点续传
    subprocess.run(['wget', '-c', url, '-O', output_file])
    
    print(f'Done with {gsm_number} ({i}/{total_gsm})')