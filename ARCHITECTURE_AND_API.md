# moljam：整体架构与对外 API（当前版本）

> 本文档用于后续维护/阅读：说明项目分块、模块关系，以及所有对外可调用的函数/类/方法签名。

## 1) 整体架构（Architecture）

### 1.1 项目形态

- 标准 Python 库工程
- 源码路径：`moljam/src/moljam/`
- 建议本地开发使用：`PYTHONPATH=moljam/src` 或 `pip install -e moljam`

### 1.2 顶层入口（Top-level API）

文件：`moljam/src/moljam/__init__.py`

- `moljam.MoleculeDBScorer`：评分引擎门面（scorer façade）
- `moljam.score_database`：读 CSV → 评分 → 返回结果（helper API）

### 1.3 scoring（评分层）

- CSV 入口函数：`moljam/src/moljam/scoring/api.py`
  - 负责：读取 CSV、校验列是否存在、实例化 scorer、执行全量检查、产出报告。
- 评分门面类：`moljam/src/moljam/scoring/scorer.py`
  - 负责：通过 mixin 组合所有检查项；对外暴露统一方法（例如 `run_all_checks()`）。
- 检查项分块：`moljam/src/moljam/scoring/checks/`
  - 设计：每个领域一个“组合壳（*ChecksMixin）”，实现拆到同目录多个小文件（`*_*.py`）。
  - 主要领域：
    - `structural*`：SMILES 有效性/一致性/手性
    - `data_quality*`：标签一致性 + 数据一致性/异常值
    - `annotation*`：实验信息质量（`include_experimental_info=True` 时才计入总分）
    - `chemical_space*`：多样性/类药性
    - `distribution*`：数据规模 + 分布/均衡性
    - `orchestration*`：汇总、归一化、低分惩罚（barrel effect）
    - `cleaning*`：按规则清洗并可保存
    - `report.py`：产出文本报告
- 公共依赖/工具：
  - `moljam/src/moljam/scoring/_common.py`：公共 imports/并行/警告控制等
  - `moljam/src/moljam/scoring/chem.py`：RDKit 相关 helper（指纹、scaffold、QED、并行 worker）

### 1.4 viz（可视化层）

- 可视化门面类：`moljam/src/moljam/viz/visualizer.py`
  - 通过 mixin 组合各类 plot；自身几乎不包含逻辑。
- 核心流程：`MoleculeDBVisualizer.add_database()`
  - 内部调用 `moljam.scoring.api.score_database()`。
  - 将 `scorer/score/report/output_dir/activity_cols/class_cols` 缓存到 `self.scoring_results`，供所有 plot 使用。
- mixins 分块：`moljam/src/moljam/viz/mixins/`
  - `chemical*`：scaffold / properties / QED / t-SNE / activity consistency 等
  - `structural*`：chirality / duplication / contradictions / smiles quality 等
  - `sankey*`：echarts 版本 + matplotlib fallback
  - `distribution*`：balance / sizes / activity / class 等
  - `bars.py / radar.py / quality.py / reporting.py / useful_columns.py / annotation_quality.py`
- 绘图依赖懒加载：`moljam/src/moljam/viz/plotting.py`（matplotlib/seaborn）
- Sankey 可选后端封装：`moljam/src/moljam/viz/sankey_backend.py`

### 1.5 classification（列分类器）

- 路径：`moljam/src/moljam/classification/`
- 用途：主要被 `scoring` 的 annotation 质量检查调用，用于识别“useful/excluded/unknown columns”。

---

## 2) 对外 API 全量签名（All Public Signatures）

> 说明：实例方法签名已去掉 `self`；以下为当前代码实际可调用的 public 方法集合。

### 2.1 函数：`score_database`

```
score_database(
    db_path,
    smiles_col='smiles',
    activity_cols=None,
    class_cols=None,
    experimental_method_cols=None,
    id_col=None,
    name_col=None,
    time_col=None,
    use_parallel=True,
    include_experimental_info=False,
)
```

### 2.2 类：`MoleculeDBScorer`

构造：
```
MoleculeDBScorer(
    df,
    smiles_col='smiles',
    activity_cols=None,
    class_cols=None,
    experimental_method_cols=None,
    id_col=None,
    name_col=None,
    time_col=None,
    use_parallel=True,
    include_experimental_info=False,
)
```

公开方法：
```
analyze_chemical_diversity()
analyze_data_balance_and_distribution()
analyze_druglikeness()
apply_low_score_penalty()
calculate_quality_score(error_rate, max_score=10, threshold_low=10, threshold_high=50)
calculate_total_scores()
check_annotation_quality()
check_data_consistency_and_reliability()
check_data_size()
check_experimental_methods()
check_label_consistency()
check_representation_consistency()
check_stereochemistry()
check_time_label_availability()
clean_database(
    remove_invalid_smiles=True,
    remove_undefined_stereochemistry=True,
    remove_conflicting_labels=True,
    remove_consistent_duplicates=True,
    verbose=True,
)
compute_fingerprints_batch(mols)
get_detailed_report()
run_all_checks()
run_missing_checks()
save_cleaned_database(
    output_path,
    remove_invalid_smiles=True,
    remove_undefined_stereochemistry=True,
    remove_conflicting_labels=True,
    remove_consistent_duplicates=True,
    verbose=True,
)
score_count_based_issues(issue_count, dataset_size, max_score=10, severity='medium')
score_low_count_issues(
    issue_count,
    dataset_size,
    max_score=10,
    count_thresholds=[0, 1, 5, 10, 20, 50],
    rate_threshold=0.01,
)
validate_smiles()
```

### 2.3 类：`MoleculeDBVisualizer`

构造：
```
MoleculeDBVisualizer(output_dir='plots_bbbp_test', use_parallel=True)
```

公开方法：
```
add_database(
    name,
    file_path,
    smiles_col='smiles',
    activity_cols=None,
    class_cols=None,
    experimental_method_cols=None,
    id_col=None,
    name_col=None,
    time_col=None,
    include_experimental_info=False,
)

generate_all_plots(
    show_values=True,
    save_prefix=None,
    plot_chirality_for=None,
    top_scaffolds_n=10,
    activity_col=None,
    class_col=None,
    tsne_samples=None,
    include_tsne=True,
)

generate_comprehensive_report(output_pdf='comprehensive_report.pdf')

plot_activity_boxplot_comparison(figsize=(14, 8), save_path=None)
plot_activity_kde_distribution(figsize=(14, 8), save_path=None)
plot_activity_raincloud(activity_col=None, figsize=(14, 8), save_path=None)
plot_all_databases_combined_tsne(n_samples_per_db=None, figsize=(14, 12), save_path=None, plot_individual=True)
plot_annotation_quality_breakdown(db_name, figsize=(16, 8), save_path=None)
plot_bars_multi_canvas(show_values=True, figsize=(20, 12), save_path=None)
plot_bars_single_canvas(show_values=True, figsize=(20, 8), save_path=None)
plot_before_after_cleaning(figsize=(14, 8), save_path=None)
plot_chemical_space_tsne(n_samples=None, figsize=(12, 10), save_path=None)
plot_class_distribution(db_name, class_col=None, figsize=(10, 8), save_path=None)
plot_class_entropy_comparison(figsize=(12, 6), save_path=None)
plot_contradictory_label_molecules(db_name, save_path=None)
plot_contradictory_label_molecules_individual(db_name, output_dir=None)
plot_data_imbalance_lorenz(figsize=(16, 16), save_path=None)
plot_data_quality_metrics(figsize=(20, 10), save_path=None)
plot_database_size_comparison(figsize=(12, 6), save_path=None)
plot_duplicate_activity_consistency(figsize=(14, 8), save_path=None)
plot_invalid_smiles_comparison(mode='both', figsize=(14, 8), save_path=None)
plot_label_conflict_heatmap(figsize=(12, 10), save_path=None)
plot_lipinski_violations(figsize=(16, 8), save_path=None)
plot_molecular_properties_distribution(figsize=(20, 13), save_path=None)
plot_non_standardized_smiles_comparison(mode='both', figsize=(14, 8), save_path=None)
plot_qed_distribution(figsize=(12, 6), save_path=None)
plot_quality_heatmap(figsize=(14, 10), save_path=None)
plot_quality_waterfall(figsize=(14, 8), save_path=None)
plot_radar(figsize=(11, 11), save_path=None)
plot_sankey_both_versions(db_name, save_path=None)
plot_sankey_for_database(db_name, save_path=None)
plot_sankey_for_database_matplotlib_fallback(db_name, figsize=(16, 12), save_path=None)
plot_structural_duplication_molecules(db_name, save_path=None)
plot_structural_duplication_molecules_individual(db_name, output_dir=None)
plot_top_scaffolds(db_name, top_n=10, figsize=(16, 10), save_path=None)
plot_undefined_chirality_comparison(mode='both', use_log=False, add_inset=False, figsize=(14, 8), save_path=None)
plot_undefined_chirality_molecules(db_name, figsize=(16, 16), save_path=None)
plot_useful_columns_analysis(db_name, figsize=(18, 10), save_path=None)
```

### 2.4 列分类器（`moljam.classification`）

#### `KeywordColumnClassifier`
构造：
```
KeywordColumnClassifier(uniqueness_threshold=0.98)
```
公开方法：
```
classify_column(col_name, col_data)
classify_columns(df, smiles_col='Smiles')
is_excluded_column(col_name, col_data)
is_useful_column(col_name, col_data)
```

#### `StatisticalColumnClassifier`
构造：
```
StatisticalColumnClassifier()
```
公开方法：
```
classify_by_features(col_name, col_data, features=None)
classify_columns(df, smiles_col='Smiles')
extract_features(col_name, col_data)
```

#### `HybridColumnClassifier`
构造：
```
HybridColumnClassifier(keyword_weight=0.6, statistical_weight=0.4)
```
公开方法：
```
classify_column(col_name, col_data)
classify_columns(df, smiles_col='Smiles')
```

