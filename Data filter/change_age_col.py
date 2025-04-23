# -*- coding: utf-8 -*-
import csv
import re

# 打开输入文件
with open('/data2/lishuhan/GSEdir/age_reslut/X.Sample_characteristics_ch1..age..csv', mode='r') as infile:
    reader = csv.DictReader(infile)
    
    # 创建一个新的列表来保存更改后的行
    new_rows = []

    # 遍历文件的每一行
    for row in reader:
        content = row["Content"]

        # 搜索数字后跟Y, y, years等
        match = re.search(r'(\d+(\.\d+)?)\s*(y|Y|years?)', content)
        if match:
            row["Content"] = match.group(1)
        
        # 搜索数字后跟M, m, months等
        match = re.search(r'(\d+(\.\d+)?)\s*(m|M|months?)', content)
        if match:
            number = float(match.group(1))
            row["Content"] = "{:.2f}".format(number / 12)
        
        # 搜索数字后跟D, d, days等
        match = re.search(r'(\d+(\.\d+)?)\s*(d|D|days?)', content)
        if match:
            number = float(match.group(1))
            row["Content"] = "{:.2f}".format(number / 365)

        # 将处理后的行添加到新列表中
        new_rows.append(row)

# 写入到新的CSV文件
with open('/data2/lishuhan/GSEdir/age_reslut/output.csv', mode='w', newline='') as outfile:
    writer = csv.DictWriter(outfile, fieldnames=["Content", "Row_Number"])
    writer.writeheader()
    writer.writerows(new_rows)
