# PAS-associated region mining adapted to the MedCAI/AttriMIL framework.
# Exploratory in-sample analysis; this script does not train or validate a classifier.
import os
import glob
import gc
import argparse
from pathlib import Path
import h5py
import pandas as pd
import numpy as np
import cv2
import openslide
import torch
import torch.nn.functional as F
from tqdm import tqdm

# --- 修正 Matplotlib 导入顺序 ---
import matplotlib
matplotlib.use('Agg')  # 必须在 plt 导入之前
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# --- 科学计算与图像处理 ---
from scipy.stats import rankdata

# --- 自定义模块 ---
from models.AttriMIL import AttriMIL
from dataloader import Generic_MIL_Dataset
# 基础设置
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

def calculate_full_weights(dataset):
    labels = dataset.slide_data['label'].values.astype(int)
    counts = np.bincount(labels, minlength=3)
    print(f"\n[Data Info] Full Dataset Counts: Low: {counts[0]}, Mid: {counts[1]}, High: {counts[2]}")
    weights = torch.FloatTensor([2.5, 2.5, 0.8]).to(device)
    return weights

def mine_regions_from_checkpoint(dataset, save_path, checkpoint_path, feature_dim=1024, n_classes=3):
    if not os.path.exists(save_path): os.makedirs(save_path)

    # 定义输出文件路径
    output_csv = os.path.join(save_path, 'mining_results_top100_final.csv')

    # 1. 模型初始化
    model = AttriMIL(dim=feature_dim, n_classes=n_classes)
    model.to(device)

    # 2. 加载权重
    print(f"\n--- Loading Trained Model from {checkpoint_path} ---")
    model.load_state_dict(torch.load(checkpoint_path, map_location=device, weights_only=True))

    # 3. 特征挖掘
    print("\n--- Starting Feature Mining & Top-K Extraction ---")
    model.eval()

    # --- 逻辑增强：断点续传检查 ---
    processed_ids = []
    extraction_results = []
    if os.path.exists(output_csv):
        existing_df = pd.read_csv(output_csv)
        processed_ids = existing_df['slide_id'].tolist()
        extraction_results = existing_df.to_dict('records')
        print(f"Detected existing CSV. Skipping {len(processed_ids)} processed slides.")

    with torch.no_grad():
        for i in tqdm(range(len(dataset))):
            slide_id = str(dataset.slide_data['slide_id'].iloc[i]).strip()

            # 跳过已处理的 ID
            if slide_id in processed_ids:
                continue

            try:
                batch_data = dataset[i]
                # 显式使用 .detach() 确保不追踪计算图
                data = batch_data[0].to(device).detach()
                label = int(batch_data[1])
                group_name = dataset.slide_data['group'].iloc[i]

                # --- 路径查找逻辑 ---
                patho_key, h5_path = dataset.get_feature_path(i)

                # --- 模型推理 ---
                outputs = model(data)
                logits, attribute_score = outputs[0], outputs[3]

                # 立即将结果转为 CPU 标量/数组，释放显存张量
                probs = F.softmax(logits, dim=1).cpu().numpy()[0]

                # 修复原代码中的缩进和重复逻辑
                if len(attribute_score.shape) == 3:
                    current_scores = attribute_score[0, label].detach().cpu().numpy()
                else:
                    current_scores = attribute_score[label].detach().cpu().numpy()

                with h5py.File(h5_path, 'r') as f:
                    coords_h5 = f['coords'][:]

                if len(coords_h5) != len(current_scores):
                    raise ValueError(f'{slide_id}: {len(coords_h5)} coordinates but {len(current_scores)} scores')

                top_k_indices = np.argsort(current_scores)[-100:][::-1]

                # 提取数据并转化为字符串存储
                extraction_results.append({
                    'slide_id': slide_id,
                    'group': group_name,
                    'patho': patho_key,
                    'true_label': label,
                    'prob_Low': probs[0],
                    'prob_Mid': probs[1],
                    'prob_High': probs[2],
                    'top_k_coords': ";".join([f"{int(x)},{int(y)}" for x, y in coords_h5[top_k_indices]]),
                    'top_k_scores': ";".join([f"{s:.4f}" for s in current_scores[top_k_indices]])
                })

                # --- 极致内存释放 ---
                del data, outputs, logits, attribute_score, current_scores, coords_h5, batch_data

                # 每 10 张图强制清理一次，平衡速度与内存
                if i % 10 == 0:
                    gc.collect()
                    torch.cuda.empty_cache()

                # 每 50 张图保存一次中间结果，防止中途 Killed 导致白跑
                if i % 50 == 0:
                    pd.DataFrame(extraction_results).to_csv(output_csv, index=False)

            except Exception as e:
                print(f"\n[Error] Skipping {slide_id}: {e}")
                continue

    # 最终保存
    df = pd.DataFrame(extraction_results)
    df.to_csv(output_csv, index=False)
    print(f"\n[Success] All results saved to: {output_csv}")

    return model

def get_strong_paper_cmap():
    """
    创建一个颜色浓郁、对比度强的自定义色带。
    高分区域使用完全不透明的红色，低分区域快速淡化。
    """
    # 颜色节点：(位置, (R, G, B, Alpha))
    # Alpha 增加得更快，让颜色看起来更“厚实”
    colors = [
        (0.00, (0.0, 0.0, 0.8, 0.0)),  # 底部：蓝色且完全透明
        (0.25, (0.0, 0.5, 1.0, 0.3)),  # 较低：青色，半透明
        (0.50, (0.0, 1.0, 0.0, 0.6)),  # 中间：绿色，较清晰
        (0.75, (1.0, 1.0, 0.0, 0.8)),  # 较高：黄色，很清晰
        (1.00, (1.0, 0.0, 0.0, 1.0))   # 顶部：纯红，完全不透明
    ]
    return mcolors.LinearSegmentedColormap.from_list("strong_cmap", colors)

def create_smooth_overlay_v2(scores, coords, patch_size, scale, region_size):
    """
    核心平滑逻辑：基于归一化卷积
    """
    # 转换坐标到缩略图空间
    patch_size_scaled = np.ceil(np.array([patch_size, patch_size]) * scale).astype(int)
    coords_scaled = np.floor(coords * scale).astype(int)

    h, w = region_size[1], region_size[0]
    overlay = np.zeros((h, w), dtype=float)
    counter = np.zeros((h, w), dtype=float)

    # 填充基础得分矩阵
    for idx, (x, y) in enumerate(coords_scaled):
        pw, ph = patch_size_scaled[0], patch_size_scaled[1]
        # 边界修剪确保不溢出矩阵
        y_end, x_end = min(y + ph, h), min(x + pw, w)
        overlay[y:y_end, x:x_end] += scores[idx]
        counter[y:y_end, x:x_end] += 1

    # 计算初始平均值
    valid_mask = counter > 0
    overlay[valid_mask] /= counter[valid_mask]

    # --- CLAM 风格平滑 ---
    # sigma 取 patch 缩放后尺寸的 1.2 倍左右，达到云雾感
    sigma = patch_size_scaled[0] * 1.2

    # 归一化卷积平滑：解决 patch 缺失导致的边缘塌陷
    smoothed_score = cv2.GaussianBlur(overlay, (0, 0), sigmaX=sigma, sigmaY=sigma)
    smoothed_counter = cv2.GaussianBlur(counter.astype(float), (0, 0), sigmaX=sigma, sigmaY=sigma)

    # 避免除以 0，计算最终平滑图
    final_overlay = np.divide(smoothed_score, smoothed_counter,
                              out=np.zeros_like(smoothed_score),
                              where=smoothed_counter > 1e-5)

    return np.clip(final_overlay, 0, 1), smoothed_counter

def generate_trident_clam_fusion_heatmaps(dataset, model, svs_dir_dict, save_path, data_dir_dict):
    heatmap_save_dir = os.path.join(save_path, 'heatmaps_trident_fusion')
    os.makedirs(heatmap_save_dir, exist_ok=True)

    all_svs_files = {os.path.basename(f).replace('.svs', ''): f for d in svs_dir_dict.values() for f in glob.glob(os.path.join(d, "*.svs"))}
    all_h5_files = {os.path.basename(f).replace('.h5', ''): f for d in data_dir_dict.values() for f in glob.glob(os.path.join(d, "features_uni_v1", "*.h5"))}

    device = next(model.parameters()).device
    model.eval()
    cmap = plt.get_cmap('RdBu_r')
    patch_size_level0 = 256

    with torch.no_grad():
        for i in tqdm(range(len(dataset))):
            # 将 try 放在循环的最开始，确保任何一步出错都能跳到下一个 slide
            try:
                # A. 获取 Slide ID
                if isinstance(dataset, torch.utils.data.Subset):
                    real_index = dataset.indices[i]
                    slide_id = str(dataset.dataset.slide_data['slide_id'].iloc[real_index]).strip()
                else:
                    slide_id = str(dataset.slide_data['slide_id'].iloc[i]).strip()

                save_name = os.path.join(heatmap_save_dir, f"{slide_id}_final.pdf")
                if os.path.exists(save_name): continue

                h5_path, svs_path = all_h5_files.get(slide_id), all_svs_files.get(slide_id)
                if not h5_path or not svs_path: continue

                # B. 获取数据与推理
                batch_data = dataset[i]
                data, label = batch_data[0].to(device), int(batch_data[1])
                outputs = model(data)
                # 立即将结果转为 numpy 并释放 Tensor
                attr_scores = outputs[3].detach().cpu()
                if len(attr_scores.shape) == 3:
                    raw_scores = attr_scores[0, label].numpy()
                else:
                    raw_scores = attr_scores[label].numpy()
                scores_norm = rankdata(raw_scores, 'average') / len(raw_scores)
                del data, outputs, attr_scores # 显式删除推理中间件

                with h5py.File(h5_path, 'r') as f:
                    coords = f['coords'][:]

                # C. 读取 Slide 信息
                slide = openslide.OpenSlide(svs_path)
                vis_level = slide.get_best_level_for_downsample(32)
                downsample = slide.level_downsamples[vis_level]
                scale = 1.0 / downsample
                t_w, t_h = slide.level_dimensions[vis_level]

                # D. 平滑处理与融合
                overlay, density_mask = create_smooth_overlay_v2(
                    scores_norm, coords, patch_size_level0, scale, (t_w, t_h)
                )

                img_rgb = np.array(slide.read_region((0, 0), vis_level, (t_w, t_h)).convert("RGB"))

                # 组织掩码逻辑
                gray = cv2.cvtColor(img_rgb, cv2.COLOR_RGB2GRAY)
                _, tissue_mask_raw = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
                kernel = np.ones((5,5), np.uint8)
                tissue_mask_bin = cv2.morphologyEx(tissue_mask_raw, cv2.MORPH_OPEN, kernel)
                tissue_mask_soft = cv2.GaussianBlur(tissue_mask_bin.astype(float), (15, 15), 0) / 255.0

                heatmap_color = (cmap(overlay)[:, :, :3] * 255).astype(np.uint8)
                valid_region_mask = (density_mask > 0.05) & (tissue_mask_bin > 0)
                combined_alpha = valid_region_mask[:, :, np.newaxis] * tissue_mask_soft[:, :, np.newaxis] * 0.6

                final_img = (img_rgb * (1 - combined_alpha) + heatmap_color * combined_alpha).astype(np.uint8)
                final_img[tissue_mask_bin == 0] = 255

                # E. 绘图与保存
                fig, ax = plt.subplots(figsize=(12, 12))
                ax.imshow(final_img)
                ax.axis('off')
                sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1))
                # 增加 colorbar 对象的句柄，方便设置参数
                cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
                cbar.ax.tick_params(labelsize=12) # AI 中易于识别的字体大小
                # 可以在 AI 中作为独立文本编辑
                cbar.set_label('Attribute score percentile', rotation=270, labelpad=15, fontsize=14)
                plt.savefig(save_name, format='pdf', dpi=300, bbox_inches='tight', transparent=True,pad_inches=0)

                # --- 极致释放内存 (这些 del 必须在 try 块内，且在 slide 还在时) ---
                plt.cla()
                plt.clf()
                plt.close('all') # 彻底关闭所有画布

                if 'slide' in locals():
                    slide.close()

                # 删除所有本轮循环产生的局部变量
                del final_img, img_rgb, overlay, density_mask, tissue_mask_bin, gray, tissue_mask_soft, scores_norm, raw_scores

                # 每一张图都强制回收
                gc.collect()
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()

            except Exception as e:
                print(f"\n[Error] Skipping slide {i} due to: {e}")
                # 清理逻辑也放入 except 保证容错
                gc.collect()
                continue

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Exploratory mining of PAS-associated WSI regions from a saved AttriMIL checkpoint")
    parser.add_argument('--metadata-csv', required=True, help='CSV with slide_id, case_id, group and optional Patho/source/project')
    parser.add_argument('--luad-features', required=True, help='LUAD directory containing features_uni_v1/')
    parser.add_argument('--lusc-features', required=True, help='LUSC directory containing features_uni_v1/')
    parser.add_argument('--checkpoint', type=Path, default=Path(__file__).with_name('full_model_epoch_78.pt'))
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).with_name('results'))
    parser.add_argument('--luad-slides', help='Directory containing LUAD .svs files')
    parser.add_argument('--lusc-slides', help='Directory containing LUSC .svs files')
    args = parser.parse_args()

    data_dir_dict = {'LUAD': args.luad_features, 'LUSC': args.lusc_features}
    dataset = Generic_MIL_Dataset(
        csv_path=args.metadata_csv,
        data_dir=data_dir_dict,
        shuffle=False,
        seed=1314,
        label_dict={'Low': 0, 'Medium': 1, 'High': 2},
        label_col='group',
    )
    dataset.load_from_h5(True)
    model = mine_regions_from_checkpoint(
        dataset=dataset,
        save_path=str(args.output_dir),
        checkpoint_path=args.checkpoint,
        feature_dim=1024,
    )

    if args.luad_slides and args.lusc_slides:
        generate_trident_clam_fusion_heatmaps(
            dataset=dataset,
            model=model,
            svs_dir_dict={'LUAD': args.luad_slides, 'LUSC': args.lusc_slides},
            save_path=str(args.output_dir),
            data_dir_dict=data_dir_dict,
        )
    elif args.luad_slides or args.lusc_slides:
        parser.error('Provide both --luad-slides and --lusc-slides to generate heatmaps')
