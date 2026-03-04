"""
纯统计特征的列分类器
不依赖列名关键词，仅通过数据的统计特征判断

核心思想：
- 唯一标识符: 高唯一性 + 特定模式
- 分子名称: 高唯一性 + 文本特征
- 可计算描述符: 数值型 + 特定分布（难以区分）
- 测量值: 数值型 + 高变异 + 小数
- 分类标签: 低唯一性 + 明确类别
"""

import pandas as pd
import numpy as np
import re


class StatisticalColumnClassifier:
    """
    基于统计特征的列分类器
    """

    def __init__(self):
        self.features_cache = {}

    def extract_features(self, col_name, col_data):
        """
        提取列的统计特征

        Returns:
            dict: 特征字典
        """
        n_total = len(col_data)
        n_valid = col_data.notna().sum()
        n_missing = n_total - n_valid

        features = {
            # 基本统计
            'n_total': n_total,
            'n_valid': n_valid,
            'n_missing': n_missing,
            'missing_rate': n_missing / n_total if n_total > 0 else 0,

            # 唯一性特征
            'n_unique': col_data.nunique(),
            'uniqueness_ratio': col_data.nunique() / n_total if n_total > 0 else 0,

            # 数据类型
            'dtype': str(col_data.dtype),
            'is_numeric': pd.api.types.is_numeric_dtype(col_data),
            'is_integer': pd.api.types.is_integer_dtype(col_data),
            'is_float': pd.api.types.is_float_dtype(col_data),
            'is_object': pd.api.types.is_object_dtype(col_data),
        }

        # 数值型特征
        if features['is_numeric'] and n_valid > 0:
            valid_data = col_data.dropna()

            features['mean'] = valid_data.mean()
            features['median'] = valid_data.median()
            features['std'] = valid_data.std()
            features['min'] = valid_data.min()
            features['max'] = valid_data.max()
            features['range'] = features['max'] - features['min']

            # 变异系数
            features['cv'] = features['std'] / abs(features['mean']) if features['mean'] != 0 else 0

            # 偏度和峰度
            try:
                features['skewness'] = valid_data.skew()
                features['kurtosis'] = valid_data.kurtosis()
            except Exception:
                features['skewness'] = 0
                features['kurtosis'] = 0

            # 小数特征（针对浮点数）
            if features['is_float']:
                features['has_decimal'] = (valid_data % 1 != 0).any()
                if features['has_decimal']:
                    # 计算平均小数位数
                    decimal_counts = valid_data.apply(lambda x: len(str(x).split('.')[-1]) if '.' in str(x) else 0)
                    features['avg_decimal_places'] = decimal_counts.mean()
                else:
                    features['avg_decimal_places'] = 0
            else:
                features['has_decimal'] = False
                features['avg_decimal_places'] = 0

            # 是否有负值
            features['has_negative'] = (valid_data < 0).any()

            # 是否单调递增
            if features['is_integer'] and n_valid > 1:
                features['is_monotonic_increasing'] = valid_data.is_monotonic_increasing
            else:
                features['is_monotonic_increasing'] = False

        else:
            # 非数值型默认值
            features.update({
                'mean': None, 'median': None, 'std': None,
                'min': None, 'max': None, 'range': None,
                'cv': None, 'skewness': None, 'kurtosis': None,
                'has_decimal': False, 'avg_decimal_places': 0,
                'has_negative': False, 'is_monotonic_increasing': False
            })

        # 文本特征（针对字符串列）
        if features['is_object'] and n_valid > 0:
            valid_data = col_data.dropna().astype(str)

            # 字符串长度统计
            str_lengths = valid_data.str.len()
            features['avg_str_length'] = str_lengths.mean()
            features['std_str_length'] = str_lengths.std()
            features['min_str_length'] = str_lengths.min()
            features['max_str_length'] = str_lengths.max()

            # 是否固定长度
            features['is_fixed_length'] = (features['std_str_length'] == 0)

            # 是否有特定模式（如ID格式）
            # 检查是否是 PREFIX-NUMBER 格式
            prefix_number_pattern = r'^[A-Za-z]+[-_]?\d+$'
            features['has_id_pattern'] = valid_data.str.match(prefix_number_pattern).mean() > 0.8

            # 检查是否是UUID格式
            uuid_pattern = r'^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$'
            features['is_uuid'] = valid_data.str.match(uuid_pattern, flags=re.IGNORECASE).mean() > 0.8

        else:
            features.update({
                'avg_str_length': None, 'std_str_length': None,
                'min_str_length': None, 'max_str_length': None,
                'is_fixed_length': False, 'has_id_pattern': False,
                'is_uuid': False
            })

        # 类别分布特征（针对低唯一性列）
        if features['n_unique'] <= 100:
            value_counts = col_data.value_counts(dropna=True)
            if len(value_counts) > 0:
                features['most_common_ratio'] = value_counts.iloc[0] / n_valid if n_valid > 0 else 0

                # 类别平衡度（使用熵）
                probs = value_counts / n_valid
                entropy = -np.sum(probs * np.log2(probs + 1e-10))
                max_entropy = np.log2(len(value_counts)) if len(value_counts) > 0 else 1
                features['category_balance'] = entropy / max_entropy if max_entropy > 0 else 0
            else:
                features['most_common_ratio'] = 0
                features['category_balance'] = 0
        else:
            features['most_common_ratio'] = None
            features['category_balance'] = None

        return features

    def classify_by_features(self, col_name, col_data, features=None):
        """
        基于统计特征分类

        Returns:
            (str, str, float): (category, reason, confidence)
        """
        if features is None:
            features = self.extract_features(col_name, col_data)

        # 规则1: 完全缺失 -> 排除
        if features['missing_rate'] >= 1.0:
            return 'excluded', 'All values missing', 1.0

        # 规则2: 唯一标识符 -> 排除
        # 特征: 高唯一性(>95%) + (UUID格式 或 ID模式 或 单调递增整数)
        if features['uniqueness_ratio'] > 0.95:
            if features['is_uuid']:
                return 'excluded', 'UUID identifier', 0.95

            if features['has_id_pattern']:
                return 'excluded', 'ID pattern (PREFIX-NUMBER)', 0.90

            if features['is_integer'] and features['is_monotonic_increasing']:
                return 'excluded', 'Sequential integer ID', 0.85

            # 如果是整数且唯一性>98%，很可能是ID
            if features['is_integer'] and features['uniqueness_ratio'] > 0.98:
                return 'excluded', 'High cardinality integer (likely ID)', 0.80

            # 如果是浮点数且有小数，可能是测量值（不排除）
            if features['is_float'] and features['has_decimal']:
                # 继续后面的判断
                pass
            else:
                # 其他高唯一性情况，可能是名称或ID
                if features['is_object'] and features['avg_str_length'] > 10:
                    return 'excluded', 'High cardinality text (likely name/description)', 0.75

        # 规则3: 名称列 -> 排除
        # 特征: 唯一性50-95% + 字符串 + 长度>8 + 长度变化大
        if (features['is_object'] and
            0.5 < features['uniqueness_ratio'] < 0.95 and
            features['avg_str_length'] and features['avg_str_length'] > 8 and
            features['std_str_length'] and features['std_str_length'] > 3):
            return 'excluded', 'Text field (likely name/description)', 0.70

        # 规则4: 常量列 -> 排除
        # 特征: 唯一值=1 或 最频繁值占比>99%
        if features['n_unique'] == 1:
            return 'excluded', 'Constant column', 0.95
        if features['most_common_ratio'] and features['most_common_ratio'] > 0.99:
            return 'excluded', 'Near-constant column', 0.90

        # === 以下是有用的列 ===

        # 规则5: 二分类/多分类标签 -> 有用
        # 特征: 2-20个类别 + 类别平衡度适中
        if 2 <= features['n_unique'] <= 20:
            if features['category_balance'] and features['category_balance'] > 0.3:
                return 'useful', f'Classification label ({features["n_unique"]} balanced categories)', 0.85
            else:
                # 不平衡但也可能是标签
                return 'useful', f'Classification label ({features["n_unique"]} categories)', 0.70

        # 规则6: 分类实验条件 -> 有用
        # 特征: 20-100个类别（比标签多，但不是唯一ID）
        if 20 < features['n_unique'] <= 100:
            return 'useful', f'Categorical experimental info ({features["n_unique"]} categories)', 0.65

        # 规则7: 数值测量值 -> 有用
        # 特征: 数值型 + 有小数 + CV>0.1 + 唯一性较高
        if (features['is_numeric'] and
            features['cv'] is not None and features['cv'] > 0.1):

            # 高唯一性浮点数 -> 很可能是测量值
            if features['is_float'] and features['has_decimal']:
                if features['uniqueness_ratio'] > 0.7:
                    return 'useful', f'Continuous measurement (CV={features["cv"]:.2f})', 0.85
                else:
                    return 'useful', f'Numeric measurement (CV={features["cv"]:.2f})', 0.75

            # 整数但变异较大 -> 可能是计数或离散测量
            if features['is_integer'] and features['uniqueness_ratio'] > 0.3:
                return 'useful', f'Discrete measurement (CV={features["cv"]:.2f})', 0.70

        # 规则8: 低变异数值 -> 可能是可计算描述符或批次信息
        # 这部分很难区分，保守判断为"unknown"或低置信度"useful"
        if (features['is_numeric'] and
            features['cv'] is not None and 0.01 < features['cv'] <= 0.1 and
            features['uniqueness_ratio'] > 0.5):
            return 'useful', f'Numeric feature (low variance, CV={features["cv"]:.2f})', 0.50

        # 默认: 未知
        return 'unknown', 'No clear pattern in statistical features', 0.0

    def classify_columns(self, df, smiles_col='Smiles'):
        """
        对所有列进行分类

        Args:
            df: DataFrame
            smiles_col: SMILES列名

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
            features = self.extract_features(col, col_data)
            category, reason, confidence = self.classify_by_features(col, col_data, features)

            entry = (col, reason, confidence, features)
            results[category].append(entry)

        return results
