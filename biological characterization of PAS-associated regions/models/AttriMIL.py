# Adapted from MedCAI/AttriMIL (https://github.com/MedCAI/AttriMIL), Apache-2.0.
# Modified for three PAS groups and the saved checkpoint bundled here;
# see ../ATTRIMIL_LICENSE for the upstream license.
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np

class Attn_Net_Gated(nn.Module):
    def __init__(self, L=1024, D=256, dropout=True, n_classes=1):
        super(Attn_Net_Gated, self).__init__()
        self.attention_a = [
            nn.Linear(L, D),
            nn.Tanh()]

        self.attention_b = [
            nn.Linear(L, D),
            nn.Sigmoid()]

        # 在注意力分支中加入 Dropout，防止模型只关注极少数 Patch
        if dropout:
            self.attention_a.append(nn.Dropout(0.60))
            self.attention_b.append(nn.Dropout(0.60))

        self.attention_a = nn.Sequential(*self.attention_a)
        self.attention_b = nn.Sequential(*self.attention_b)
        self.attention_c = nn.Linear(D, n_classes)

    def forward(self, x):
        a = self.attention_a(x)
        b = self.attention_b(x)
        A = a.mul(b)
        A = self.attention_c(A)  # N x n_classes
        return A, x

class AttriMIL(nn.Module):
    """
    Multi-Branch ABMIL with reinforced constraints and dropout regularization.
    """
    def __init__(self, n_classes=3, dim=512):
        super().__init__()

        # 1. Adaptor 层：增加 Dropout 强力防止特征层面的过拟合
        self.adaptor = nn.Sequential(
            nn.Linear(dim, dim // 2),
            nn.ReLU(),
            nn.Dropout(p=0.6),
            nn.Linear(dim // 2, dim),
            nn.Dropout(p=0.6)
        )

        # 2. 注意力网络：显式开启 dropout=True
        self.attention_nets = nn.ModuleList([
            Attn_Net_Gated(L=dim, D=dim // 2, dropout=True, n_classes=1)
            for _ in range(n_classes)
        ])

        # 3. 实例分类器：在分类前加入 Dropout，这是解决 Class 1 泛化差的关键
        self.classifiers = nn.ModuleList([
            nn.Sequential(
                nn.Dropout(p=0.6),
                nn.Linear(dim, 1)
            ) for _ in range(n_classes)
        ])

        self.n_classes = n_classes
        self.bias = nn.Parameter(torch.zeros(n_classes), requires_grad=True)

    def forward(self, h):
        # h: [N, dim] - N 是 Patch 的数量

        # 特征自适应映射 (Residual Connection)
        h = h + self.adaptor(h)

        # 初始化存储容器，直接在 h 所在的设备上创建
        device = h.device
        n_patches = h.size(0)

        # 用于存储每个类别的注意力权重和实例得分
        # A_raw: [n_classes, n_patches]
        # instance_score: [n_classes, n_patches]
        A_raw = []
        instance_scores = []

        for c in range(self.n_classes):
            # 计算注意力权重
            A, _ = self.attention_nets[c](h) # A: [N, 1]
            A_raw.append(A.view(1, -1))     # 变为 [1, N]

            # 计算实例级别得分
            score = self.classifiers[c](h)   # score: [N, 1]
            instance_scores.append(score.view(1, -1)) # 变为 [1, N]

        # 拼接各个类别的结果
        A_raw = torch.cat(A_raw, dim=0)           # [n_classes, N]
        instance_scores = torch.cat(instance_scores, dim=0) # [n_classes, N]

        # 计算 Attribute Score (instance_score * exp(A))
        # 使用 exp 放大注意力差异是 AttriMIL 的特征
        # attribute_score: [1, n_classes, N]
        attribute_score = (instance_scores * torch.exp(A_raw)).unsqueeze(0)

        # 计算 Bag Logits
        logits = torch.empty(1, self.n_classes).to(device)
        for c in range(self.n_classes):
            # 标准 MIL 聚合逻辑：加权平均得分 + 偏置
            # 分母是该类别所有 patch 的注意力指数和
            denom = torch.sum(torch.exp(A_raw[c])) + 1e-8
            logits[0, c] = torch.sum(attribute_score[0, c]) / denom + self.bias[c]

        # 预测结果处理
        Y_prob = F.softmax(logits, dim=1)
        Y_hat = torch.argmax(logits, dim=1)

        results_dict = {} # 可根据需要扩展

        return logits, Y_prob, Y_hat, attribute_score, results_dict
