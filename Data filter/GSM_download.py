# -*- coding: utf-8 -*-
import os
import wget
import subprocess

# ��ȡ����GSM������ļ�
with open('gds_result.txt', 'r') as f:
    lines = f.readlines()

# ��ȡ���е�GSM����
gsm_numbers = [line.split(': ')[1].split()[0] for line in lines if line.startswith('Sample')]

# ����һ���µ��ļ������洢���ص��ļ�
os.makedirs('downloaded_files', exist_ok=True)

# ��ȡGSM������������Ա�������
total_gsm = len(gsm_numbers)

# ��ÿ��GSM���룬���ض�Ӧ��txt�ļ�
for i, gsm_number in enumerate(gsm_numbers, start=1):
    print(f'Downloading {gsm_number} ({i}/{total_gsm})...')
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm_number}&targ=self&form=text&view=quick"
    output_file = f'downloaded_files/{gsm_number}.txt'
    
    # ʹ��wget�����ļ�����ʹ��-cѡ��ʵ�ֶϵ�����
    subprocess.run(['wget', '-c', url, '-O', output_file])
    
    print(f'Done with {gsm_number} ({i}/{total_gsm})')