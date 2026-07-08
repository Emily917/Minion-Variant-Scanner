import matplotlib.pyplot as plt
import csv
import os
import glob

# ============================================================
# 1. 自动找到最新的 CSV 文件
# ============================================================
csv_files = glob.glob("*_k*_top*.csv")
if not csv_files:
    print("❌ 错误: 找不到任何 *_k*_top*.csv 文件")
    print("   请先运行 ./scanner 生成数据文件")
    exit()

# 按修改时间排序，取最新的
csv_files.sort(key=os.path.getmtime, reverse=True)
csv_filename = csv_files[0]
print(f"📂 读取数据文件: {csv_filename}")

# ============================================================
# 2. 从 CSV 读取数据
# ============================================================
kmer_data = []
filename = "unknown"
k_value = "unknown"
topN = "unknown"

with open(csv_filename, 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        if not row:
            continue
        if row[0].startswith('# filename:'):
            filename = row[0].replace('# filename:', '').strip()
        elif row[0].startswith('# k:'):
            k_value = row[0].replace('# k:', '').strip()
        elif row[0].startswith('# topN:'):
            topN = row[0].replace('# topN:', '').strip()
        elif row[0] != 'kmer' and len(row) == 2:
            kmer_data.append((row[0], int(row[1])))

if not kmer_data:
    print("❌ 错误: CSV 文件为空或格式不正确")
    exit()

# ============================================================
# 3. 生成动态文件名
# ============================================================
# 从 CSV 文件名提取信息：SRR390728_1_k11_top10.csv
base_name = csv_filename.replace('.csv', '')  # 去掉 .csv

# 组装图片文件名
png_filename = base_name + ".png"
print(f"📊 生成图片: {png_filename}")

# ============================================================
# 4. 画图
# ============================================================
labels, values = zip(*kmer_data)

plt.figure(figsize=(12, 6))
bars = plt.bar(labels, values, color='steelblue')

for bar, val in zip(bars, values):
    plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
             str(val), ha='center', va='bottom', fontsize=10)

plt.xlabel('K-mer Sequence', fontsize=12)
plt.ylabel('Frequency', fontsize=12)
plt.title(f'Top {topN} Most Frequent {k_value}-mers ({filename})', fontsize=14)
plt.xticks(rotation=45, ha='right')
plt.tight_layout()

# ============================================================
# 5. 保存到桌面 Minion_Results 文件夹
# ============================================================
desktop_path = os.path.expanduser("~/Desktop")
target_folder = os.path.join(desktop_path, "Minion_Results")

if not os.path.exists(target_folder):
    os.makedirs(target_folder)
    print(f"📁 创建文件夹: {target_folder}")

save_path = os.path.join(target_folder, png_filename)
plt.savefig(save_path, dpi=150)
print(f"✅ 图表已保存到: {save_path}")

# 同时在项目文件夹保存一份
plt.savefig(png_filename, dpi=150)
print(f"✅ 图表已保存到本地: {png_filename}")