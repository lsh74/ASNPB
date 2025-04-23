# -*- coding: utf-8 -*-
import os

num_sh_files = 10
num_lines_per_sh = 20000

# ��ȡ����GSM������ļ�
with open('gds_result.txt', 'r') as f:
    lines = f.readlines()

# ��ȡ���е�GSM����
gsm_numbers = [line.split(': ')[1].split()[0] for line in lines if line.startswith('Sample')]

num_sh_files_needed = len(gsm_numbers) // num_lines_per_sh
if len(gsm_numbers) % num_lines_per_sh != 0:
    num_sh_files_needed += 1

for i in range(num_sh_files_needed):
    # ���Ŀ¼�Ƿ���ڣ��������򴴽�
    if not os.path.exists(f'wget{i+1}'):
        os.makedirs(f'wget{i+1}')

    with open(f'wget{i+1}.sh', 'w') as f:
        start = i * num_lines_per_sh
        end = min((i + 1) * num_lines_per_sh, len(gsm_numbers))
        for idx, gsm_number in enumerate(gsm_numbers[start:end]):
            cmd = f'wget -c -t 0 "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm_number}&targ=self&form=text&view=brief" -O "./wget{i+1}/{gsm_number}.txt" &\n'
            f.write(cmd)

            # ÿ100����������һ��wait����
            if (idx+1) % 100 == 0:
                f.write('wait\n')

        # ������һ��wait�����Ӳ���ִ�е�����
        f.write("wait\n")
        f.write("cat wget_commands.txt | xargs -n 1 -P 200 sh\n")
