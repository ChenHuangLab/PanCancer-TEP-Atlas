import os
import json
import argparse
import pandas as pd
from collections import Counter

# 基础路径
groups = ['High', 'Medium', 'Low']

# 细胞类型映射表 (根据你的 HoverNet 训练集设定，以下为常见映射示例)
# 0: 背景/其他, 1: 肿瘤细胞, 2: 间质细胞, 3: 淋巴细胞, 4: 巨噬细胞, 5: 中性粒细胞
# --- PanNuke 模型细胞类型映射表 ---
# 0: nolabe (背景/无标签)
# 1: neopla (Neoplastic - 肿瘤细胞)
# 2: inflam (Inflammatory - 炎症细胞)
# 3: connec (Connective - 结缔组织/间质细胞)
# 4: necros (Necrosis - 坏死区域)
# 5: no-neo (Non-Neoplastic Epithelial - 非肿瘤上皮细胞)

type_map = {
    0: 'Other',
    1: 'Neoplastic',
    2: 'Inflammatory',
    3: 'Connective',
    4: 'Necrosis',
    5: 'Non-Neoplastic'
}

def analyze_hovernet_results(base_dir):
    all_stats = []

    for group in groups:
        group_path = os.path.join(base_dir, group)
        if not os.path.exists(group_path):
            continue

        print(f"正在分析 {group} 组...")

        # 遍历 WSI 文件夹
        for wsi_name in os.listdir(group_path):
            wsi_path = os.path.join(group_path, wsi_name, 'json')
            if not os.path.exists(wsi_path):
                continue

            # 遍历 Top5 JSON 文件
            for json_file in os.listdir(wsi_path):
                if not json_file.endswith('.json'):
                    continue

                with open(os.path.join(wsi_path, json_file), 'r') as f:
                    data = json.load(f)

                # 提取所有细胞核的类型
                # HoverNet json 结构通常为 {'nuc': {'id': {'type': int, ...}}}
                nucs = data.get('nuc', {})
                types = [nuc_info.get('type') for nuc_info in nucs.values()]
                counts = Counter(types)

                # 记录该 patch 的统计结果
                patch_stat = {
                    'Group': group,
                    'WSI': wsi_name,
                    'Patch': json_file,
                    'Total_Nuclei': len(types)
                }

                # 填充各个细胞类型的数量
                for type_idx, type_name in type_map.items():
                    patch_stat[type_name] = counts.get(type_idx, 0)

                all_stats.append(patch_stat)

    return pd.DataFrame(all_stats)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Summarize PanNuke HoVer-Net nuclei calls for selected patches')
    parser.add_argument('--hovernet-dir', required=True)
    parser.add_argument('--output-csv', required=True)
    args = parser.parse_args()
    df_results = analyze_hovernet_results(args.hovernet_dir)
    if df_results.empty:
        raise SystemExit('No HoVer-Net JSON files found')
    group_summary = df_results.groupby('Group')[list(type_map.values())].mean()
    print(group_summary)
    print(group_summary.div(group_summary.sum(axis=1), axis=0) * 100)
    df_results.to_csv(args.output_csv, index=False)
