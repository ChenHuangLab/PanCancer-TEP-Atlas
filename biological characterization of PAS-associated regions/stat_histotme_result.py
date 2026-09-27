import os
import argparse
import pandas as pd
import numpy as np

# --- 1. 路径配置（请根据实际情况修改） ---

def extract_tme_features(mining_csv, base_dir, top_k=5):
    # 读取 AttriMIL 挖掘出的 top patch 结果
    df_mining = pd.read_csv(mining_csv)
    all_matched_data = []

    print(f"开始处理 {len(df_mining)} 条 WSI 记录...")

    for _, row in df_mining.iterrows():
        slide_id = row['slide_id']
        group = row['group']

        # 构建 HistoTME 结果文件的路径
        # 路径结构: HistoTME_result/slide_id/uni/slide_id_5fold.csv
        tme_file = os.path.join(base_dir, slide_id, 'uni', f"{slide_id}_5fold.csv")

        if not os.path.exists(tme_file):
            # print(f"跳过: 找不到 HistoTME 文件 {slide_id}")
            continue

        # 读取该 WSI 的 HistoTME 评分表
        try:
            df_tme = pd.read_csv(tme_file)
        except Exception as e:
            print(f"读取文件失败 {slide_id}: {e}")
            continue

        # 确保坐标列是整数类型
        df_tme['x'] = df_tme['x'].astype(int)
        df_tme['y'] = df_tme['y'].astype(int)

        # 提取 top_k 的坐标字符串并解析
        # 原始格式: "16896,59392;19456,58880;..."
        coords_list = row['top_k_coords'].split(';')[:top_k]

        for i, coord_str in enumerate(coords_list):
            try:
                # 解析出当前 patch 的 x, y
                target_x, target_y = map(int, coord_str.split(','))

                # 在 HistoTME 表中匹配坐标
                matched_patch = df_tme[(df_tme['x'] == target_x) & (df_tme['y'] == target_y)]

                if not matched_patch.empty:
                    # 获取该 patch 的所有特征评分
                    patch_stats = matched_patch.iloc[0].to_dict()
                    # 补充元数据
                    patch_stats['group_label'] = group
                    patch_stats['slide_id'] = slide_id
                    patch_stats['patch_rank'] = i + 1
                    all_matched_data.append(patch_stats)
                # else:
                #    print(f"未匹配到坐标: {slide_id} at ({target_x}, {target_y})")

            except ValueError:
                continue

    return pd.DataFrame(all_matched_data)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Match top PAS-associated coordinates to HistoTME signatures')
    parser.add_argument('--mining-csv', required=True)
    parser.add_argument('--histotme-dir', required=True, help='Run once for each cohort-specific HistoTME root')
    parser.add_argument('--matched-csv', required=True)
    parser.add_argument('--summary-csv', required=True)
    parser.add_argument('--top-k', type=int, default=5)
    args = parser.parse_args()
    full_stats_df = extract_tme_features(args.mining_csv, args.histotme_dir, top_k=args.top_k)
    if full_stats_df.empty:
        raise SystemExit('No coordinates matched HistoTME results')
    exclude_cols = ['Unnamed: 0', 'x', 'y', 'slide_id', 'group_label', 'patch_rank']
    score_cols = [c for c in full_stats_df.select_dtypes(include='number').columns if c not in exclude_cols]
    summary = full_stats_df.groupby('group_label')[score_cols].mean().T
    full_stats_df.to_csv(args.matched_csv, index=False)
    summary.to_csv(args.summary_csv)
    print(summary.head(10))
