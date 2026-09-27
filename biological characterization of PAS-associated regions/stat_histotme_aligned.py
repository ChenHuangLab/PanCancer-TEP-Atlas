import os
import argparse
import pandas as pd
import numpy as np
from tqdm import tqdm

# --- 1. 路径配置 ---
hover_base_path = None
histotme_base_dir = None
mining_csv_path = None
output_path = None

def analyze_aligned_histotme():
    df_mining = pd.read_csv(mining_csv_path)
    all_matched_data = []

    groups = ['High', 'Medium', 'Low']

    print(">>> 正在启动样本对齐提取程序...")

    for group in groups:
        group_path = os.path.join(hover_base_path, group)
        if not os.path.exists(group_path):
            continue

        # 获取同步基准：HoVer-Net 处理过的 WSI
        processed_wsis = [d for d in os.listdir(group_path) if os.path.isdir(os.path.join(group_path, d))]

        for wsi_name in tqdm(processed_wsis, desc=f"Processing {group}"):
            tme_file = os.path.join(histotme_base_dir, wsi_name, 'uni', f"{wsi_name}_5fold.csv")
            if not os.path.exists(tme_file):
                continue

            mining_row = df_mining[df_mining['slide_id'] == wsi_name]
            if mining_row.empty:
                continue

            # 读取 HistoTME 数据
            try:
                df_tme = pd.read_csv(tme_file)
                df_tme['x'] = df_tme['x'].astype(int)
                df_tme['y'] = df_tme['y'].astype(int)
            except:
                continue

            # 解析坐标
            coords_list = str(mining_row.iloc[0]['top_k_coords']).split(';')[:5]

            for i, coord_str in enumerate(coords_list):
                try:
                    target_x, target_y = map(int, coord_str.split(','))
                    matched_patch = df_tme[(df_tme['x'] == target_x) & (df_tme['y'] == target_y)]

                    if not matched_patch.empty:
                        patch_stats = matched_patch.iloc[0].to_dict()
                        patch_stats['Group'] = group
                        patch_stats['Slide_ID'] = wsi_name
                        all_matched_data.append(patch_stats)
                except:
                    continue

    return pd.DataFrame(all_matched_data)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Summarize HistoTME signatures on slides also processed by HoVer-Net')
    parser.add_argument('--hovernet-dir', required=True)
    parser.add_argument('--histotme-dir', required=True)
    parser.add_argument('--mining-csv', required=True)
    parser.add_argument('--matched-csv', required=True)
    parser.add_argument('--summary-csv', required=True)
    args = parser.parse_args()
    hover_base_path = args.hovernet_dir
    histotme_base_dir = args.histotme_dir
    mining_csv_path = args.mining_csv
    output_path = args.summary_csv
    df_aligned = analyze_aligned_histotme()
    if df_aligned.empty:
        raise SystemExit('No HistoTME rows matched the selected HoVer-Net slides')

    # 统计信息计算
    exclude_cols = ['Unnamed: 0', 'x', 'y', 'Slide_ID', 'Group', 'Rank']
    score_cols = [c for c in df_aligned.select_dtypes(include='number').columns if c not in exclude_cols]

    group_counts = df_aligned.groupby('Group')['Slide_ID'].nunique()
    group_summary = df_aligned.groupby('Group')[score_cols].mean().T

    # --- 终端打印部分 ---
    print("\n" + "="*60)
    print("📊 HistoTME 样本对齐统计分析报告")
    print("="*60)
    print(f"参与统计的 WSI 样本分布:")
    for g in ['High', 'Medium', 'Low']:
        count = group_counts.get(g, 0)
        print(f"  - {g:8}: {count} WSIs")

    print("-" * 60)
    print("各组 Signature 评分对比 (Top 15 特征):")
    # 按照分值差异较大的特征进行排序显示（或者直接显示前15个）
    print(group_summary.head(30).round(4))
    print("-" * 60)

    # 保存文件
    df_aligned.to_csv(args.matched_csv, index=False)
    group_summary.to_csv(output_path)

    print(f"✅ 结果已保存至: {output_path}")
    print("="*60)
