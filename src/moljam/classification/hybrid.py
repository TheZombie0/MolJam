"""
混合分类器：统计特征 + 列名关键词
结合两者优势，提高准确率和鲁棒性

决策逻辑：
1. 列名关键词提供强烈信号 → 高置信度决策
2. 统计特征提供辅助验证 → 中等置信度决策
3. 两者一致 → 提高置信度
4. 两者矛盾 → 降低置信度，标记为需要人工审核
"""

from .keyword import KeywordColumnClassifier
from .statistical import StatisticalColumnClassifier


class HybridColumnClassifier:
    """
    混合分类器：结合列名关键词和统计特征
    """

    def __init__(self, keyword_weight=0.6, statistical_weight=0.4):
        """
        初始化混合分类器

        Args:
            keyword_weight: 关键词方法的权重
            statistical_weight: 统计方法的权重
        """
        self.keyword_classifier = KeywordColumnClassifier()
        self.statistical_classifier = StatisticalColumnClassifier()

        self.keyword_weight = keyword_weight
        self.statistical_weight = statistical_weight

    def classify_column(self, col_name, col_data):
        """
        混合分类单个列

        Returns:
            dict: {
                'category': str,
                'reason': str,
                'confidence': float,
                'keyword_vote': (category, reason, conf),
                'statistical_vote': (category, reason, conf),
                'agreement': bool,
                'details': dict
            }
        """
        # 方法1: 关键词方法
        kw_result = self.keyword_classifier.classify_column(col_name, col_data)
        kw_category = kw_result['category']
        kw_reason = kw_result['reason']
        kw_conf = kw_result['confidence']

        # 方法2: 统计方法
        stat_features = self.statistical_classifier.extract_features(col_name, col_data)
        stat_category, stat_reason, stat_conf = self.statistical_classifier.classify_by_features(
            col_name, col_data, stat_features
        )

        # 融合决策
        final_category, final_reason, final_conf, agreement = self._merge_votes(
            (kw_category, kw_reason, kw_conf),
            (stat_category, stat_reason, stat_conf),
            col_name, stat_features
        )

        return {
            'category': final_category,
            'reason': final_reason,
            'confidence': final_conf,
            'keyword_vote': (kw_category, kw_reason, kw_conf),
            'statistical_vote': (stat_category, stat_reason, stat_conf),
            'agreement': agreement,
            'details': kw_result['details']
        }

    def _merge_votes(self, kw_vote, stat_vote, col_name, stat_features):
        """
        融合两个分类器的投票

        Args:
            kw_vote: (category, reason, confidence)
            stat_vote: (category, reason, confidence)
            col_name: 列名
            stat_features: 统计特征

        Returns:
            (category, reason, confidence, agreement)
        """
        kw_cat, kw_reason, kw_conf = kw_vote
        stat_cat, stat_reason, stat_conf = stat_vote

        # 情况1: 两者一致 → 提高置信度
        if kw_cat == stat_cat:
            final_conf = kw_conf * self.keyword_weight + stat_conf * self.statistical_weight
            # 一致性加成
            final_conf = min(1.0, final_conf * 1.1)
            final_reason = f"{kw_reason} | Confirmed by statistics: {stat_reason}"
            return kw_cat, final_reason, final_conf, True

        # 情况2: 关键词有强烈信号（高置信度）→ 以关键词为准
        if kw_conf >= 0.8:
            # 特殊处理：可计算描述符
            if 'Calculable descriptor' in kw_reason:
                # 统计上可能看起来像测量值，但关键词明确指出是可计算的
                final_reason = f"{kw_reason} (overrides statistical: {stat_reason})"
                return 'excluded', final_reason, kw_conf, False

            # 实验信息关键词
            if 'Experimental information' in kw_reason or 'assay' in str(col_name).lower():
                final_reason = f"{kw_reason} (overrides statistical: {stat_reason})"
                return 'useful', final_reason, kw_conf * 0.9, False

        # 情况3: 统计有强烈信号但关键词说unknown → 以统计为准
        if stat_conf >= 0.8 and kw_cat == 'unknown':
            final_reason = f"{stat_reason} | No keyword match"
            return stat_cat, final_reason, stat_conf * 0.8, False

        # 情况4: 特殊规则 - 处理稀疏标签（MUV问题）
        # 如果统计说是"Near-constant"但关键词识别为有用，且唯一值=2，缺失率高
        if (stat_cat == 'excluded' and 'constant' in stat_reason.lower() and
            kw_cat == 'useful' and
            stat_features['n_unique'] == 2 and
            stat_features['missing_rate'] > 0.5):
            # 这很可能是稀疏二分类标签
            final_reason = f"{kw_reason} | Sparse binary label (high missing rate)"
            return 'useful', final_reason, kw_conf * 0.85, False

        # 情况5: 两者矛盾，无强烈信号 → 加权平均，降低置信度
        if kw_cat != 'unknown' and stat_cat != 'unknown':
            # 按权重选择
            if kw_conf * self.keyword_weight > stat_conf * self.statistical_weight:
                final_cat = kw_cat
                final_reason = f"{kw_reason} | Conflicts with statistical: {stat_reason}"
                final_conf = kw_conf * self.keyword_weight * 0.7  # 冲突惩罚
            else:
                final_cat = stat_cat
                final_reason = f"{stat_reason} | Conflicts with keyword: {kw_reason}"
                final_conf = stat_conf * self.statistical_weight * 0.7
            return final_cat, final_reason, final_conf, False

        # 情况6: 一方说unknown → 采用另一方
        if kw_cat == 'unknown' and stat_cat != 'unknown':
            return stat_cat, f"{stat_reason} | No keyword signal", stat_conf * 0.8, False
        if stat_cat == 'unknown' and kw_cat != 'unknown':
            return kw_cat, f"{kw_reason} | Weak statistical signal", kw_conf * 0.8, False

        # 默认: 都是unknown
        return 'unknown', 'No clear signal from either method', 0.0, True

    def classify_columns(self, df, smiles_col='Smiles'):
        """
        对所有列进行分类

        Returns:
            dict: 分类结果
        """
        results = {
            'useful': [],
            'excluded': [],
            'unknown': []
        }

        # 查找SMILES列
        smiles_col_actual = None
        for col in df.columns:
            if col.lower() == smiles_col.lower():
                smiles_col_actual = col
                break

        for col in df.columns:
            if smiles_col_actual and col == smiles_col_actual:
                continue

            col_data = df[col]
            result = self.classify_column(col, col_data)

            category = result['category']
            entry = (
                col,
                result['reason'],
                result['confidence'],
                result['details'],
                result['agreement']
            )
            results[category].append(entry)

        return results
