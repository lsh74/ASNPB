# -*- coding: utf-8 -*-
import csv
import re

# �������ļ�
with open('/data2/lishuhan/GSEdir/age_reslut/X.Sample_characteristics_ch1..age..csv', mode='r') as infile:
    reader = csv.DictReader(infile)
    
    # ����һ���µ��б���������ĺ����
    new_rows = []

    # �����ļ���ÿһ��
    for row in reader:
        content = row["Content"]

        # �������ֺ��Y, y, years��
        match = re.search(r'(\d+(\.\d+)?)\s*(y|Y|years?)', content)
        if match:
            row["Content"] = match.group(1)
        
        # �������ֺ��M, m, months��
        match = re.search(r'(\d+(\.\d+)?)\s*(m|M|months?)', content)
        if match:
            number = float(match.group(1))
            row["Content"] = "{:.2f}".format(number / 12)
        
        # �������ֺ��D, d, days��
        match = re.search(r'(\d+(\.\d+)?)\s*(d|D|days?)', content)
        if match:
            number = float(match.group(1))
            row["Content"] = "{:.2f}".format(number / 365)

        # ������������ӵ����б���
        new_rows.append(row)

# д�뵽�µ�CSV�ļ�
with open('/data2/lishuhan/GSEdir/age_reslut/output.csv', mode='w', newline='') as outfile:
    writer = csv.DictWriter(outfile, fieldnames=["Content", "Row_Number"])
    writer.writeheader()
    writer.writerows(new_rows)
