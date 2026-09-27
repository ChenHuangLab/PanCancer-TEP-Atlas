import os
import json
import argparse
import cv2
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.stats import mannwhitneyu

# --- 1. 配置区 ---
parser = argparse.ArgumentParser(description='Summarize morphology from PanNuke HoVer-Net JSON outputs')
parser.add_argument('--hovernet-dir', required=True)
parser.add_argument('--output-csv', required=True)
parser.add_argument('--summary-csv', required=True)
parser.add_argument('--mpp', type=float, default=0.5, help='Micrometers per pixel of the analyzed patches')
parser.add_argument('--patch-size', type=int, default=256)
args = parser.parse_args()
RESULT_ROOT = args.hovernet_dir
OUTPUT_CSV = args.output_csv
SUMMARY_CSV = args.summary_csv

# 物理参数 (20x, MPP=0.5)
MPP = args.mpp
PATCH_SIZE = args.patch_size
PATCH_AREA_MM2 = ((PATCH_SIZE * MPP) / 1000) ** 2

# PanNuke 映射表
type_map = {
    1: 'Neoplastic',
    2: 'Inflammatory',
    3: 'Connective',
    4: 'Necrosis',
    5: 'Non-Neoplastic'
}

metrics = ['Density_ln', 'Ratio', 'Mean_Area_um2', 'Median_Perimeter_um']

def process_json(json_path):
    """处理单个 JSON 文件提取细胞几何特征"""
    with open(json_path, 'r') as f:
        data = json.load(f)
    nucs = data.get('nuc', {})
    patch_cells = []
    for nuc_id, info in nucs.items():
        c_type = info['type']
        if c_type not in type_map: continue

        contour = np.array(info['contour']).astype(np.int32)
        area_um2 = cv2.contourArea(contour) * (MPP**2)
        perimeter_um = cv2.arcLength(contour, True) * MPP

        patch_cells.append({
            'type': type_map[c_type],
            'area': area_um2,
            'perimeter': perimeter_um
        })
    return patch_cells

# --- 2. 核心逻辑：形态学计算 ---
print(">>> 步骤 1: 正在从 JSON 提取形态学特征并计算 Patch 指标...")
all_results = []
groups = ['High', 'Medium', 'Low']

for group in groups:
    group_path = os.path.join(RESULT_ROOT, group)
    if not os.path.exists(group_path): continue

    slides = [d for d in os.listdir(group_path) if os.path.isdir(os.path.join(group_path, d))]
    for slide_id in tqdm(slides, desc=f"Processing {group}"):
        json_dir = os.path.join(group_path, slide_id, 'json')
        if not os.path.exists(json_dir): continue

        wsi_cells = []
        json_files = [f for f in os.listdir(json_dir) if f.endswith('.json')]

        for j_file in json_files:
            cells = process_json(os.path.join(json_dir, j_file))
            wsi_cells.extend(cells)

        if not wsi_cells: continue
        df_wsi = pd.DataFrame(wsi_cells)
        total_nuc = len(df_wsi)

        for t_idx, t_name in type_map.items():
            type_df = df_wsi[df_wsi['type'] == t_name]
            count = len(type_df)

            # 计算密度与比例
            total_area = len(json_files) * PATCH_AREA_MM2
            density_ln = np.log((count / total_area) + 1)
            ratio = count / total_nuc if total_nuc > 0 else 0

            # 计算形态学统计量
            m_area = type_df['area'].mean() if count > 0 else 0
            m_peri = type_df['perimeter'].median() if count > 0 else 0

            all_results.append({
                'Group': group,
                'Slide_ID': slide_id,
                'Cell_Type': t_name,
                'Count': count,
                'Density_ln': density_ln,
                'Ratio': ratio,
                'Mean_Area_um2': m_area,
                'Median_Perimeter_um': m_peri
            })

df_final = pd.DataFrame(all_results)
if df_final.empty:
    raise SystemExit('No HoVer-Net nuclei JSON files found')
df_final.to_csv(OUTPUT_CSV, index=False)
print(f"✅ 基础数据已保存至: {OUTPUT_CSV}")

# --- 3. 核心逻辑：组别汇总与统计检验 ---
print("\n>>> 步骤 2: 正在进行组别汇总统计与显著性检验 (High vs Low)...")
print("-" * 60)

summary_data = []

for c_type in df_final['Cell_Type'].unique():
    print(f"\n[细胞类型: {c_type}]")
    type_df = df_final[df_final['Cell_Type'] == c_type]

    # 打印每组的均值和标准差
    group_stats = type_df.groupby('Group')[metrics].agg(['mean', 'std']).round(4)
    print(group_stats)

    # 执行 Mann-Whitney U 检验 (High vs Low)
    high_vals = type_df[type_df['Group'] == 'High']['Density_ln']
    low_vals = type_df[type_df['Group'] == 'Low']['Density_ln']

    if len(high_vals) > 0 and len(low_vals) > 0:
        u_stat, p_val = mannwhitneyu(high_vals, low_vals, alternative='two-sided')
        sig_star = "*" if p_val < 0.05 else ""
        print(f"  >> Density_ln 显著性检验 (High vs Low): P = {p_val:.4f} {sig_star}")

# 保存最终汇总均值表
final_summary = df_final.groupby(['Group', 'Cell_Type'])[metrics].mean().reset_index()
final_summary.to_csv(SUMMARY_CSV, index=False)
print(f"\n✅ 汇总统计完成！最终结果查看: {SUMMARY_CSV}")
