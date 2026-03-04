"""
Column Classifier - Keyword Version (关键词版/文字版)
基于列名关键词匹配的列分类器

与 column_classifier_statistical.py（统计版）对应
两者结合使用于 column_classifier_hybrid.py（混合版）

分类方法：通过列名中的关键词判断列的类型
- id_keywords: 标识符关键词（id, key, index等）
- name_keywords: 名称关键词（name, title等）
- activity_keywords: 活性/测量值关键词（ic50, activity等）
- label_keywords: 分类标签关键词（class, type, label等）
- experimental_keywords: 实验信息关键词（assay, target等）
- prediction_keywords: 预测值关键词（calc, predicted等）

主要特点：
1. 依赖列名语义进行分类
2. 支持词边界匹配，避免误判
3. 多阶段决策机制
"""

import pandas as pd
import re


class KeywordColumnClassifier:
    """
    关键词版列分类器（文字版）
    通过列名关键词匹配进行分类
    """

    def __init__(self, uniqueness_threshold=0.98):
        """
        初始化分类器

        Args:
            uniqueness_threshold: 唯一性阈值（提高到98%）
        """
        self.uniqueness_threshold = uniqueness_threshold
        self.min_categories = 2
        self.max_categories = 50
        self.cv_threshold = 0.1

        self._init_keywords()

    def _init_keywords(self):
        """初始化关键词"""
        # 标识符关键词（需要词边界匹配）
        self.id_keywords = ['id', 'key', 'index', 'number', 'code']
        self.id_exceptions = ['assay', 'experiment', 'target', 'document', 'cell', 'bao', 'format']

        # 名称关键词
        self.name_keywords = ['name', 'title']
        self.name_exceptions = ['assay', 'target', 'cell', 'tissue']

        # 可计算的分子描述符
        self.calculable_descriptors = [
            'molecular weight', 'mw', 'exact mass',
            'logp', 'alogp', 'xlogp', 'mlogp',
            'tpsa', 'polar surface area',
            'hbd', 'h-bond donor', 'hydrogen bond donor',
            'hba', 'h-bond acceptor', 'hydrogen bond acceptor',
            'rotatable bond', 'rotation',
            'ring', 'aromatic',
            'ro5', 'lipinski', 'violation',
            'degree', 'atom count', 'heavy atom',
            'ligand efficiency'
        ]

        # 关系符号
        self.relation_keywords = ['relation', 'operator', 'comparator']

        # 元数据
        self.metadata_keywords = ['journal', 'year', 'source description', 'document', 'comment', 'properties']

        # 实验相关（排除'ph'，改为'buffer'）
        self.experimental_keywords = [
            'assay', 'experiment', 'method', 'protocol', 'technique',
            'organism', 'tissue', 'cell', 'target', 'enzyme',
            'condition', 'temperature', 'buffer',  # 移除'ph'
            'bao', 'bioassay', 'screening',
            'parameter', 'variant', 'mutation', 'subcellular', 'fraction'
        ]

        # 生物活性/测量值
        self.activity_keywords = [
            'activity', 'ic50', 'ec50', 'ki', 'kd', 'ed50',
            'potency', 'inhibition', 'binding',
            'solubility', 'measured', 'observed', 'expt',  # 'expt'是实验测量
            'concentration', 'dose', 'response',
            'pchembl', 'value', 'result', 'standard',
            'approved', 'toxic', 'tox',
            'muv-', 'phase',
            'units', 'validity', 'action'
        ]

        # 预测值/计算值（新增）
        self.prediction_keywords = [
            'calc', 'calculated', 'predicted', 'pred',
            'estimate', 'estimated', 'computed', 'simulation', 'model'
        ]

        # 分类标签
        self.label_keywords = [
            'class', 'type', 'category', 'group', 'status',
            'label', 'flag', 'disorder', 'disease', 'issue'
        ]

        # pH相关（精确匹配）
        self.ph_keywords = ['ph', 'ph value', 'ph level']

    def _keyword_match(self, text, keywords, use_word_boundary=True, exceptions=None):
        """
        改进的关键词匹配

        Args:
            text: 待匹配文本
            keywords: 关键词列表
            use_word_boundary: 是否使用词边界
            exceptions: 例外关键词列表

        Returns:
            (bool, str): (是否匹配, 匹配的关键词)
        """
        text_lower = text.lower()

        for kw in keywords:
            if use_word_boundary:
                # 使用词边界匹配
                # 注意：对于'ph'这种短词，需要特别处理
                if len(kw) <= 2:
                    # 短关键词：精确匹配或前后有分隔符
                    pattern = r'(^|[\s_-])' + re.escape(kw) + r'($|[\s_-])'
                else:
                    # 长关键词：词边界匹配
                    pattern = r'\b' + re.escape(kw) + r'\b'

                if re.search(pattern, text_lower):
                    # 检查是否有例外
                    if exceptions:
                        has_exception = any(exc in text_lower for exc in exceptions)
                        if has_exception:
                            continue
                    return True, kw
            else:
                if kw in text_lower:
                    if exceptions:
                        has_exception = any(exc in text_lower for exc in exceptions)
                        if has_exception:
                            continue
                    return True, kw

        return False, None

    def _is_float_measurement(self, col_data):
        """
        判断是否是浮点数测量值

        Args:
            col_data: 列数据

        Returns:
            bool: 是否是浮点数测量值
        """
        if not pd.api.types.is_numeric_dtype(col_data):
            return False

        valid_data = col_data.dropna()
        if len(valid_data) == 0:
            return False

        # 检查是否有小数部分
        try:
            # 如果是浮点数且有小数部分
            if pd.api.types.is_float_dtype(col_data):
                has_decimal = (valid_data % 1 != 0).any()
                return has_decimal
        except Exception:
            pass

        return False

    def is_excluded_column(self, col_name, col_data):
        """
        判断列是否应该被排除

        Args:
            col_name: 列名
            col_data: 列数据

        Returns:
            (bool, str, float): (是否排除, 原因, 置信度)
        """
        # 处理缺失值
        valid_data = col_data.dropna()
        if len(valid_data) == 0:
            return True, "All values are missing", 1.0

        col_lower = col_name.lower()

        # A. 唯一标识符（高唯一性）
        uniqueness_ratio = col_data.nunique() / len(col_data)
        if uniqueness_ratio > 0.95:
            # 特殊处理：如果是浮点数测量值，提高阈值
            if self._is_float_measurement(col_data):
                if uniqueness_ratio <= self.uniqueness_threshold:
                    # 可能是高精度测量值，不排除
                    return False, None, 0.0
            else:
                # 非浮点数，使用原阈值
                if uniqueness_ratio > 0.95:
                    return True, f"Unique identifier (cardinality: {uniqueness_ratio:.1%})", 0.9

        # B. 列名包含标识符关键词
        matched, kw = self._keyword_match(col_lower, self.id_keywords,
                                          use_word_boundary=True,
                                          exceptions=self.id_exceptions)
        if matched:
            return True, f"Identifier column ('{kw}')", 0.8

        # C. 名称列
        matched, kw = self._keyword_match(col_lower, self.name_keywords,
                                          use_word_boundary=True,
                                          exceptions=self.name_exceptions)
        if matched:
            # 进一步检查：如果唯一性不高，可能不是真正的名称列
            if uniqueness_ratio < 0.5:
                return False, None, 0.0
            return True, f"Name column ('{kw}')", 0.7

        # D. 可计算的分子描述符
        for desc in self.calculable_descriptors:
            if desc in col_lower:
                return True, f"Calculable descriptor ('{desc}')", 0.95

        # E. 关系符号
        matched, kw = self._keyword_match(col_lower, self.relation_keywords, use_word_boundary=True)
        if matched:
            return True, f"Relation/Operator column ('{kw}')", 0.9

        # F. 元数据
        matched, kw = self._keyword_match(col_lower, self.metadata_keywords, use_word_boundary=False)
        if matched:
            return True, f"Metadata column ('{kw}')", 0.8

        # G. 重复标记
        if 'duplicate' in col_lower:
            return True, "Duplicate flag column", 0.9

        return False, None, 0.0

    def is_useful_column(self, col_name, col_data):
        """
        判断列是否是有用的实验信息

        Args:
            col_name: 列名
            col_data: 列数据

        Returns:
            (bool, str, float): (是否有用, 原因, 置信度)
        """
        col_lower = col_name.lower()

        # A. pH值（精确匹配）
        matched, kw = self._keyword_match(col_lower, self.ph_keywords, use_word_boundary=True)
        if matched:
            return True, f"Experimental condition ('{kw}')", 0.9

        # B. 预测值/计算值（新增，高优先级）
        matched, kw = self._keyword_match(col_lower, self.prediction_keywords, use_word_boundary=True)
        if matched:
            return True, f"Predicted/Calculated value ('{kw}')", 0.9

        # C. 实验相关信息
        matched, kw = self._keyword_match(col_lower, self.experimental_keywords, use_word_boundary=False)
        if matched:
            return True, f"Experimental information ('{kw}')", 0.9

        # D. 生物活性/测量值
        matched, kw = self._keyword_match(col_lower, self.activity_keywords, use_word_boundary=False)
        if matched:
            return True, f"Activity/Measurement ('{kw}')", 0.85

        # E. 分类标签（2-50个唯一值）
        n_unique = col_data.nunique()
        if self.min_categories <= n_unique <= self.max_categories:
            matched, kw = self._keyword_match(col_lower, self.label_keywords, use_word_boundary=True)
            if matched:
                return True, f"Classification label ('{kw}', {n_unique} categories)", 0.8

            # 2-10个分类，即使没有关键词也可能是标签
            if 2 <= n_unique <= 10:
                return True, f"Potential label ({n_unique} categories)", 0.6

        # F. 数值型且变异较大
        if pd.api.types.is_numeric_dtype(col_data):
            valid_data = col_data.dropna()
            if len(valid_data) > 0:
                mean_val = valid_data.mean()
                std_val = valid_data.std()
                cv = std_val / abs(mean_val) if mean_val != 0 else 0
                if cv > self.cv_threshold:
                    return True, f"Numeric measurement (CV={cv:.2f})", 0.7

        return False, None, 0.0

    def classify_column(self, col_name, col_data):
        """
        多阶段分类单个列

        Args:
            col_name: 列名
            col_data: 列数据

        Returns:
            dict: {
                'category': 'useful' / 'excluded' / 'unknown',
                'reason': str,
                'confidence': float,
                'details': dict
            }
        """
        # 阶段1：检查排除
        is_excl, excl_reason, excl_conf = self.is_excluded_column(col_name, col_data)
        if is_excl:
            return {
                'category': 'excluded',
                'reason': excl_reason,
                'confidence': excl_conf,
                'details': self._get_column_details(col_name, col_data)
            }

        # 阶段2：检查有用
        is_useful, useful_reason, useful_conf = self.is_useful_column(col_name, col_data)
        if is_useful:
            return {
                'category': 'useful',
                'reason': useful_reason,
                'confidence': useful_conf,
                'details': self._get_column_details(col_name, col_data)
            }

        # 阶段3：未知
        return {
            'category': 'unknown',
            'reason': 'No matching rules',
            'confidence': 0.0,
            'details': self._get_column_details(col_name, col_data)
        }

    def _get_column_details(self, col_name, col_data):
        """获取列的详细信息"""
        n_unique = col_data.nunique()
        n_total = len(col_data)
        if n_total == 0:
            return {
                'unique_values': 0,
                'uniqueness_ratio': 0.0,
                'missing_rate': 0.0,
                'coverage': 0.0,
                'dtype': str(col_data.dtype),
                'is_numeric': pd.api.types.is_numeric_dtype(col_data),
                'is_float': pd.api.types.is_float_dtype(col_data)
            }

        missing_rate = col_data.isnull().sum() / n_total
        coverage = 1 - missing_rate

        details = {
            'unique_values': n_unique,
            'uniqueness_ratio': n_unique / n_total,
            'missing_rate': missing_rate,
            'coverage': coverage,
            'dtype': str(col_data.dtype),
            'is_numeric': pd.api.types.is_numeric_dtype(col_data),
            'is_float': pd.api.types.is_float_dtype(col_data)
        }

        return details

    def classify_columns(self, df, smiles_col='Smiles'):
        """
        对所有列进行分类

        Args:
            df: DataFrame
            smiles_col: SMILES列名

        Returns:
            dict: {
                'useful': [(col, reason, confidence, details), ...],
                'excluded': [(col, reason, confidence, details), ...],
                'unknown': [(col, reason, confidence, details), ...]
            }
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
            entry = (col, result['reason'], result['confidence'], result['details'])
            results[category].append(entry)

        return results
