#!/usr/bin/env python3
import matplotlib.pyplot as plt
import csv
import os
import glob
import sys

def find_latest_csv():
    csv_files = glob.glob("*_k*_top*.csv")
    if not csv_files:
        print("❌ Error: No CSV files found")
        print("   Please run scanner first")
        sys.exit(1)
    csv_files.sort(key=os.path.getmtime, reverse=True)
    return csv_files[0]

def main():
    csv_file = find_latest_csv()
    print(f"📂 Reading: {csv_file}")
    
    # 读取CSV文件
    high_data = []
    low_data = []
    metadata = {}
    is_new_format = False
    
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            # 读取metadata
            if row[0].startswith('# '):
                if ': ' in row[0]:
                    key_val = row[0][2:].strip().split(': ', 1)
                    if len(key_val) == 2:
                        metadata[key_val[0]] = key_val[1]
            # 检查是否有 category 列（新格式）
            elif row[0] == 'category':
                is_new_format = True
            # 读取数据 - 新格式 (category,kmer,frequency)
            elif is_new_format and len(row) == 3 and row[0] != 'category':
                category = row[0]
                kmer = row[1]
                freq = int(row[2])
                if category == 'high':
                    high_data.append((kmer, freq))
                elif category == 'low':
                    low_data.append((kmer, freq))
            # 读取数据 - 旧格式 (kmer,frequency)
            elif not is_new_format and len(row) == 2 and row[0] != 'kmer':
                try:
                    kmer = row[0]
                    freq = int(row[1])
                    high_data.append((kmer, freq))  # 旧格式全部当作高频
                except ValueError:
                    pass
    
    # 如果没有数据，尝试另一种读取方式
    if not high_data and not low_data:
        print("⚠️  Trying alternative read method...")
        import pandas as pd
        try:
            df = pd.read_csv(csv_file, comment='#')
            if 'category' in df.columns:
                high_data = df[df['category'] == 'high'][['kmer', 'frequency']].values.tolist()
                low_data = df[df['category'] == 'low'][['kmer', 'frequency']].values.tolist()
            elif 'kmer' in df.columns and 'frequency' in df.columns:
                high_data = df[['kmer', 'frequency']].values.tolist()
        except Exception as e:
            print(f"Error reading with pandas: {e}")
    
    if not high_data and not low_data:
        print("❌ No data found in CSV file")
        print("   Please ensure the CSV file has the correct format")
        return
    
    # 创建图形
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # 绘制高频
    if high_data:
        # 按频率排序
        high_data.sort(key=lambda x: x[1], reverse=True)
        kmers = [item[0] for item in high_data]
        freqs = [item[1] for item in high_data]
        colors = plt.cm.Reds_r([i/len(high_data) for i in range(len(high_data))])
        bars = axes[0].bar(kmers, freqs, color=colors)
        axes[0].set_title(f'High-Frequency K-mers (Top {len(high_data)})', 
                         fontsize=14, fontweight='bold')
        axes[0].set_xlabel('K-mer Sequence', fontsize=12)
        axes[0].set_ylabel('Frequency', fontsize=12)
        axes[0].tick_params(axis='x', rotation=45)
        
        for bar, val in zip(bars, freqs):
            axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                        str(val), ha='center', va='bottom', fontsize=9)
    else:
        axes[0].text(0.5, 0.5, 'No high-frequency data', 
                    ha='center', va='center', transform=axes[0].transAxes, fontsize=14)
        axes[0].set_title('High-Frequency K-mers', fontsize=14)
    
    # 绘制低频
    if low_data:
        low_data.sort(key=lambda x: x[1])
        kmers = [item[0] for item in low_data]
        freqs = [item[1] for item in low_data]
        colors = plt.cm.Blues_r([i/len(low_data) for i in range(len(low_data))])
        bars = axes[1].bar(kmers, freqs, color=colors)
        axes[1].set_title(f'Low-Frequency K-mers (Top {len(low_data)})', 
                         fontsize=14, fontweight='bold')
        axes[1].set_xlabel('K-mer Sequence', fontsize=12)
        axes[1].set_ylabel('Frequency', fontsize=12)
        axes[1].tick_params(axis='x', rotation=45)
        
        for bar, val in zip(bars, freqs):
            axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                        str(val), ha='center', va='bottom', fontsize=9)
    else:
        axes[1].text(0.5, 0.5, 'No low-frequency data', 
                    ha='center', va='center', transform=axes[1].transAxes, fontsize=14)
        axes[1].set_title('Low-Frequency K-mers', fontsize=14)
    
    # 主标题
    filename = metadata.get('filename', 'unknown')
    k_val = metadata.get('k', 'unknown')
    fig.suptitle(f'K-mer Distribution Analysis\n{filename} (k={k_val})', 
                 fontsize=16, fontweight='bold')
    
    plt.tight_layout()
    
    # 保存
    base_name = csv_file.replace('.csv', '')
    png_filename = base_name + '.png'
    
    desktop = os.path.expanduser("~/Desktop")
    target = os.path.join(desktop, "Minion_Results")
    os.makedirs(target, exist_ok=True)
    
    save_path = os.path.join(target, png_filename)
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"✅ Graph saved: {save_path}")
    
    plt.savefig(png_filename, dpi=150, bbox_inches='tight')
    print(f"✅ Local copy: {png_filename}")
    
    # 显示图形
    print("📊 Displaying plot...")
    plt.show()

if __name__ == "__main__":
    main()