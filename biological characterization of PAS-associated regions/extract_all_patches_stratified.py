import os
import argparse
import pandas as pd
import openslide
from tqdm import tqdm

# --- 配置区 ---
CSV_PATH = None
SVS_SEARCH_DIRS = []

# 更改路径，避免覆盖之前的全量结果，方便对比
PATCH_SAVE_DIR = None

PATCH_SIZE = 256
PATCHES_PER_WSI = 5
STRATIFY_PERCENT = 0.20

# 组别与概率列名映射
# 注意：请确认 CSV 中的列名是否为 prob_Low, prob_Mid, prob_High
group_info = {
    'Low': {'folder': 'Low', 'prob_col': 'prob_Low'},
    'Medium': {'folder': 'Medium', 'prob_col': 'prob_Mid'},
    'High': {'folder': 'High', 'prob_col': 'prob_High'}
}

def find_svs_path(slide_id):
    for root_dir in SVS_SEARCH_DIRS:
        potential_path = os.path.join(root_dir, f"{slide_id}.svs")
        if os.path.exists(potential_path):
            return potential_path
    return None

def extract_stratified_wsis():
    if not os.path.exists(CSV_PATH):
        print(f"[Error] 未找到 CSV 文件: {CSV_PATH}")
        return

    df = pd.read_csv(CSV_PATH)
    print(f"原始数据共包含 {len(df)} 条记录。")

    for csv_label, info in group_info.items():
        folder_name = info['folder']
        prob_col = info['prob_col']

        # 在同一 PAS 组内按训练样本上的概率排序，仅用于选择待表征切片。
        group_df = df[df['group'] == csv_label].copy()

        # 保留原分析中的 20% 比例及至少 10 张切片规则。
        n_to_extract = max(10, int(len(group_df) * STRATIFY_PERCENT))
        stratified_df = group_df.sort_values(by=prob_col, ascending=False).head(n_to_extract)

        print(f"\n--- {folder_name} 组: 原始 {len(group_df)} 个, 筛选后 {len(stratified_df)} 个 (Top {STRATIFY_PERCENT*100}%) ---")

        os.makedirs(os.path.join(PATCH_SAVE_DIR, folder_name), exist_ok=True)

        for _, row in tqdm(stratified_df.iterrows(), total=len(stratified_df), desc=f"提取 {folder_name}"):
            slide_id = row['slide_id']
            svs_path = find_svs_path(slide_id)
            if svs_path is None: continue

            sample_patch_dir = os.path.join(PATCH_SAVE_DIR, folder_name, slide_id)
            os.makedirs(sample_patch_dir, exist_ok=True)

            try:
                slide = openslide.OpenSlide(svs_path)
                coords_list = str(row['top_k_coords']).split(';')[:PATCHES_PER_WSI]

                for i, coord_str in enumerate(coords_list):
                    if not coord_str or ',' not in coord_str: continue
                    x, y = map(int, coord_str.split(','))
                    patch = slide.read_region((x, y), 0, (PATCH_SIZE, PATCH_SIZE)).convert('RGB')
                    patch_filename = f"top{i+1}_{slide_id}_{x}_{y}.png"
                    patch.save(os.path.join(sample_patch_dir, patch_filename))

                slide.close()
            except Exception as e:
                with open("stratified_errors.log", "a") as log_f:
                    log_f.write(f"{slide_id}: {str(e)}\n")

    print(f"\n✅ 分层 Patch 提取完成！存储路径: {PATCH_SAVE_DIR}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract top-scoring patches for exploratory characterization')
    parser.add_argument('--mining-csv', required=True)
    parser.add_argument('--luad-slides', required=True)
    parser.add_argument('--lusc-slides', required=True)
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    CSV_PATH = args.mining_csv
    SVS_SEARCH_DIRS = [args.luad_slides, args.lusc_slides]
    PATCH_SAVE_DIR = args.output_dir
    extract_stratified_wsis()
