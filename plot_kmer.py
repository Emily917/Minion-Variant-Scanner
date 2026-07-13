import matplotlib.pyplot as plt
import csv
import os
import glob

# 1. 自动找到最新的 CSV 文件 （Automatically locate the latest CSV file.）
csv_files = glob.glob("*_k*_top*.csv")
if not csv_files:
    print("❌ error: can't find any *_k*_top*.csv file")
    print("   please run python plot_kmer.py to generate files")
    exit()

# 按修改时间排序，取最新的 （Sort by modification time and select the most recent one.）
csv_files.sort(key=os.path.getmtime, reverse=True)
csv_filename = csv_files[0]
print(f"📂 read documents file: {csv_filename}")

# 2. 从 CSV 读取数据 （Read data from CSV）
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
    print("❌ error: CSV file is blank/Incorrect format")
    exit()

# 3. 生成动态文件名 （Generate dynamic file names）
# 从 CSV 文件名提取信息 （Extract information from CSV filenames）
base_name = csv_filename.replace('.csv', '')  # 去掉 .csv

png_filename = base_name + ".png"
print(f"📊 generate picture : {png_filename}") #设置生成的图片文件名 （Set the filename for the generated image.）

labels, values = zip(*kmer_data) # draw the graph

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

# 5. 保存到桌面 Minion_Results 文件夹
desktop_path = os.path.expanduser("~/Desktop")
target_folder = os.path.join(desktop_path, "Minion_Results")

if not os.path.exists(target_folder):
    os.makedirs(target_folder)
    print(f"📁 创建文件夹: {target_folder}")

save_path = os.path.join(target_folder, png_filename)
plt.savefig(save_path, dpi=150)
print(f"✅ graph saved at: {save_path}")

# 同时在项目文件夹保存一份，显示在 VS code 屏幕上 (At the same time, save a copy in the project folder and display it on the VS Code screen.)
plt.savefig(png_filename, dpi=150)
print(f"✅ graph saved nameis : {png_filename}")