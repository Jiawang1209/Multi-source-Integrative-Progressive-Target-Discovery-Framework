# 单 CAS 的 Figure 1 分析流程说明

## 文档目的

这份文档用于定义当前 `Figure 1` 的分析流程，并明确一个核心约束：

- 一个 `CAS` 编号就是一个独立项目
- 一个 `CAS` 对应一个 `case_dir`
- 一个 `CAS` 对应一个命名为 `CAS_xxxxx_Date` 的项目目录
- 一个 `CAS` 对应一条完整的 `a -> b -> c -> d -> e` 分析链
- 一个 `CAS` 对应一套独立的最终结果

后续如果要把这套流程开发成包，这份文档可以作为流程说明和设计基线。

## 核心原则

这套流程不能把多个化合物混在一个 case 结果里。

正确的项目单位应该是：

- `CAS_491-71-4_2026-03-21` 是一个独立 case
- `CAS_117-02-2_2026-03-21` 是另一个独立 case
- 每个 case 都有自己的原始输入、中间文件、图、模型和结果表

如果一个下游结果表或者预测文件里同时出现多个分子，那么它不能作为单个 CAS case 的正式最终结果。它最多只能算开发过程中的临时产物，不应该作为最终交付。

## Figure 1 的标准定义

对于单个 CAS 分子，`Figure 1` 的分析结构应为：

1. `a`：`ChEMBL / PPB2 / Swiss / SEA` 四个来源的人类靶标搜索流程示意图
2. `b`：四个数据库人类靶标的 Venn 图
3. `c`：对四库靶标并集做 `GO + KEGG enrichment`，再筛选肝脏疾病相关结果并可视化
4. `d`：对关键靶标基因再次进行 `KEGG` 分析，并用 `circlize` 画和弦图
5. `e`：严格从 `d` 图的关键 KEGG 基因出发，构建 `Chemprop multi-task` 模型并进行靶标排序和过滤

## 当前脚本对应关系

当前流程主要由 [`scripts`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts) 下这些脚本支撑：

- [`build_case_sources.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/build_case_sources.py)
- [`fetch_chembl_targets.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/fetch_chembl_targets.py)
- [`fetch_ppb2_targets.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/fetch_ppb2_targets.py)
- [`run_figure1_bcd.R`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/run_figure1_bcd.R)
- [`build_chemprop_targets_from_kegg.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/build_chemprop_targets_from_kegg.py)
- [`prepare_chemprop_data.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/prepare_chemprop_data.py)
- [`train_chemprop_multitask.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/train_chemprop_multitask.py)
- [`evaluate_chemprop_model.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/evaluate_chemprop_model.py)
- [`analyze_chemprop_predictions.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/analyze_chemprop_predictions.py)
- [`plot_figure1e.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/plot_figure1e.py)

## 推荐的 case 目录结构

每个 CAS 都应该对应一个带日期的独立项目目录，例如：

- `CAS_491-71-4_2026-03-21`

推荐结构如下：

```text
CAS_491-71-4_2026-03-21/
  a_source_collection/
  b_venn/
  c_go_kegg/
  d_kegg_circlize/
  e_chemprop/
```

每一步的结果都只能写入自己的子目录。

### 文件夹 A：`a_source_collection`

- 原始 `Swiss` 结果
- 原始 `SEA` 结果
- 原始 `ChEMBL` 靶标结果
- 原始 `PPB2` 靶标结果
- `idmapping.tsv`
- 来源汇总文件

### 文件夹 B：`b_venn`

- `venn.rds`
- `venn_inputs.json`
- `1.p_venn.pdf`

### 文件夹 C：`c_go_kegg`

- `2.GO_result_2.xlsx`
- `2.KEGG_result_2.xlsx`
- `GO_KEGG_plot.xlsx`
- `2.GO_KEGG.pdf`

### 文件夹 D：`d_kegg_circlize`

- `key_targets_from_kegg.csv`
- `3.KEGG_circos.pdf`

### 文件夹 E：`e_chemprop`

- `*_targets_from_d.csv`
- `chemprop_long.csv`
- `chemprop_multitask.csv`
- `task_summary.csv`
- `inference_template.csv`
- `training/model_0/best.pt`
- `predictions.csv`
- `test_metrics_by_target.csv`
- `ranked_targets.csv`
- 过滤后的 `Figure 1e`

## 目录规则

这个目录约定应该被视为硬约束：

- 一个 CAS 只能对应一个 case 目录
- case 目录命名格式必须是 `CAS_xxxxx_Date`
- 每个 panel 步骤都必须有自己的子文件夹
- 子文件夹顺序必须严格是 `a -> b -> c -> d -> e`
- 不同 CAS 的文件绝不能混写

## Step A：靶标来源定义

### 目标

针对一个真实 CAS 分子，先在线解析化合物身份信息，再从四个在线来源收集人类靶标蛋白：

- `ChEMBL`
- `PPB2`
- `SwissTargetPrediction`
- `SEA`

### 输入

- 一个真实 `CAS`
- 在线自动解析化合物身份信息
- 在线自动获取：
  - 化合物名称
  - 分子式
  - SMILES
  - InChIKey

### Human-only 约束

所有靶标必须限定为 `Homo sapiens`。

### 在线优先原则

正式包的预期行为应该是：

- 输入真实 `CAS`
- 在线解析化合物身份信息
- 基于解析得到的结构信息，在线提交到四个数据库
- 将四库结果统一写入 `a_source_collection`

本地缓存文件可以作为开发调试或回退手段存在，但不应成为正式包的默认主流程。

### 当前在线身份解析来源

当前化合物身份解析基于 `PubChem PUG REST`。

解析时应优先获得：

- `Title` 或 `IUPACName`
- `MolecularFormula`
- `CanonicalSMILES` 或 `ConnectivitySMILES`
- `InChIKey`

当前实现方式：

- `ChEMBL`：保留 `Organism == "Homo sapiens"`
- `SEA`：保留 `Target ID` 以 `_HUMAN` 结尾的记录
- `PPB2`：只保留成功标准化到人类基因符号的结果
- `Swiss`：提取返回的目标名，并尽可能标准化为基因层级名称

### 相关脚本

- [`src/miptd/discovery.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/src/miptd/discovery.py)
- [`fetch_chembl_targets.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/fetch_chembl_targets.py)
- [`fetch_ppb2_targets.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/fetch_ppb2_targets.py)
- [`fetch_swiss_targets.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/fetch_swiss_targets.py)
- [`fetch_sea_targets.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/fetch_sea_targets.py)
- [`build_case_sources.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/build_case_sources.py)

### 输出

每个 case 的四库标准化靶标表，以及对应的来源计数汇总，统一写入 `a_source_collection`。

## Step B：Venn 图

### 目标

用四个数据库的人类靶标集合绘制 Venn 图：

- `Swiss`
- `SEA`
- `ChEMBL`
- `PPB2`

### 关键规则

Venn 图基于的是四个来源的靶标集合本身，而不是后面进一步筛选过的靶标。

### 相关脚本

- [`run_figure1_bcd.R`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/run_figure1_bcd.R)

### 输出

- `venn.rds`
- `venn_inputs.json`
- `1.p_venn.pdf`

### 解释

Venn 图用于展示四个数据库之间的交集和重叠情况。

但后续富集分析并不是只使用 overlap，而是使用四库靶标的**并集**。

## Step C：GO 和 KEGG 富集分析

### 目标

对四个数据库靶标的并集进行富集分析。

### 输入

四库靶标的并集：

- `Swiss`
- `SEA`
- `ChEMBL`
- `PPB2`

### 方法

执行：

- `GO enrichment`
- `KEGG enrichment`

然后根据肝脏疾病和肝脏生物学相关关键词，对结果进行筛选。

这里的筛选逻辑是：

1. 先对并集靶标做完整 GO/KEGG 富集
2. 再根据命令行传入的 `--disease-keywords` / `--disase-keywords` 构造匹配规则
3. 用该规则去匹配富集结果中的 `Description`
4. 优先保留匹配到的条目作为 `c` 图基础
5. 再从被选中的 KEGG 条目中提取关键基因用于 `d` 图

当前关键词逻辑包含这类模式：

- `liver`
- `hepatic`
- `lipid`
- `fatty`
- `cholesterol`
- `bile`
- `oxidative`
- `apopt`
- `inflamm`
- `insulin`
- `AMPK`
- `PPAR`
- `TNF`
- `NF-kappa`
- `FoxO`
- `MAPK`
- `PI3K`
- `AKT`
- `ABC`

如果用户在命令行里传入自己的关键词，例如：

- `--disease-keywords "NAFLD,liver,hepatic"`

那么程序会优先使用用户提供的关键词，而不是只依赖默认肝脏相关模式。

如果一个关键词都没有匹配到，则流程会回退到按显著性和计数排序的 top 条目，以保证 `c` 和 `d` 仍能生成。

### 相关脚本

- [`run_figure1_bcd.R`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/run_figure1_bcd.R)

### 输出

- `2.GO_result_2.xlsx`
- `2.KEGG_result_2.xlsx`
- `GO_KEGG_plot.xlsx`
- `2.GO_KEGG.pdf`

### 解释

这一步的作用是把原始靶标并集收缩到具有肝脏疾病相关意义的功能和通路层面。

## Step D：关键 KEGG 基因与和弦图

### 目标

从 Step C 选出的关键 KEGG 通路中提取关键基因，并用 `circlize` 画出通路-基因关系图。

### 输入

Step C 中筛出的关键 KEGG 通路。

### 方法

对于每一条选中的 KEGG 通路：

- 提取该通路对应的基因
- 合并所有通路中的基因
- 建立 pathway-gene 对应关系表
- 用 chord diagram 进行可视化

### 相关脚本

- [`run_figure1_bcd.R`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/run_figure1_bcd.R)

### 输出

- `key_targets_from_kegg.csv`
- `3.KEGG_circos.pdf`

### 解释

这个 `key_targets_from_kegg.csv` 才是 panel `e` 的正确上游输入。

重要规则：

- panel `e` 必须从 `key_targets_from_kegg.csv` 开始
- 不能直接从 `b` 图 overlap 结果或某个早期快捷靶标表直接进入 `Chemprop`

## Step E：Chemprop 多任务建模

### 目标

基于 panel `d` 得到的关键 KEGG 基因，构建 `Chemprop multi-task` 模型，并对单个 CAS 查询分子进行靶标打分。

### 必须遵守的规则

panel `e` 必须**严格承接 panel `d`**。

正确链条是：

- `key_targets_from_kegg.csv`
- `build_chemprop_targets_from_kegg.py`
- `prepare_chemprop_data.py`
- `train_chemprop_multitask.py`
- `evaluate_chemprop_model.py`
- `analyze_chemprop_predictions.py`
- `plot_figure1e.py`

### Step E1：把 d 图关键基因转成 Chemprop 候选靶标

相关脚本：

- [`build_chemprop_targets_from_kegg.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/build_chemprop_targets_from_kegg.py)

输出：

- `*_targets_from_d.csv`

这个表表示的是来自 panel `d` 的 Chemprop 候选靶标。
该输出应写入 `e_chemprop`。

### Step E2：准备 Chemprop 训练数据

相关脚本：

- [`prepare_chemprop_data.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/prepare_chemprop_data.py)

主要逻辑：

- 使用本地 target catalog 和 `idmapping`
- 将基因映射到 `target_chembl_id`
- 优先保留 `SINGLE PROTEIN`
- 从 ChEMBL 获取活性数据
- 再次限定为 `Homo sapiens`
- 保留标准活性类型，如 `IC50 / Ki / Kd / EC50`
- 聚合到每个 molecule-target 对应一个 `pChEMBL`
- 只保留样本量足够的可训练任务

输出：

- `chemprop_long.csv`
- `chemprop_multitask.csv`
- `task_summary.csv`
- `inference_template.csv`

这些输出应写入 `e_chemprop`。

### Step E3：训练 Chemprop 模型

相关脚本：

- [`train_chemprop_multitask.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/train_chemprop_multitask.py)

当前实现：

- 回归任务
- 数据切分方式为 `SCAFFOLD_BALANCED`
- 保存 train/validation/test 的 SMILES 切分
- 保存最佳模型和最终预测结果

输出：

- `training/model_0/best.pt`
- `predictions.csv`
- `run_metadata.json`

这些输出应写入 `e_chemprop`。

### Step E4：逐靶点评估模型质量

相关脚本：

- [`evaluate_chemprop_model.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/evaluate_chemprop_model.py)

该脚本会对保存下来的测试集重新预测，并导出逐靶点测试指标。

输出：

- `test_predictions.csv`
- `test_metrics_by_target.csv`
- `test_metrics_summary.json`

这些输出应写入 `e_chemprop`。

当前重点指标包括：

- `RMSE`
- `MAE`
- `R2`
- `n_test_molecules`

### Step E5：靶标排序与模型质量过滤

相关脚本：

- [`analyze_chemprop_predictions.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/analyze_chemprop_predictions.py)

当前结果表导出的关键列包括：

- `predicted_pchembl`
- `n_training_molecules`
- `test_rmse`
- `test_mae`
- `test_r2`
- `n_test_molecules`
- `filter_basis`
- `filter_status`

当前过滤规则是：

- `n_training_molecules >= 100`
- `n_test_molecules >= 30`
- `test_r2 >= 0.30`

其目的是把结果拆分为：

- 原始预测靶标
- 高置信靶标

### Step E6：绘制 Figure 1e

相关脚本：

- [`plot_figure1e.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/plot_figure1e.py)

当前支持两种图：

- 未过滤版 `Figure 1e`
- 过滤后版 `Figure 1e`，只显示 `filter_status = keep` 的靶标

推荐用于正式结果展示的是：

- 过滤后版 `Figure 1e`

该图应写入 `e_chemprop`。

## 当前开发过程中总结出的关键经验

### 1. 单 CAS 隔离是硬约束

最重要的流程要求就是：

- 一个 CAS
- 一个项目
- 一套独立结果

后续无论开发成脚本还是包，都应该主动防止多个分子被混进同一个预测表或排序表。

### 2. panel e 必须从 panel d 出发

这条规则现在已经确定。

Chemprop 的正确上游输入只能是：

- `key_targets_from_kegg.csv`

而不能被替换成：

- 单纯的 Venn overlap 结果
- 某个临时共识靶标快捷表

### 3. 预测分数高不等于可作为最终结论

如果某个靶标的预测值很高，但该任务本身的模型测试指标很差，那么它不应该直接作为论文最终高置信靶标。

典型剔除条件如：

- `R2 < 0.30`
- `n_test_molecules < 30`

### 4. Human-only 过滤必须做两层

人类限定至少要在两层实现：

- target catalog 层
- activity record 层

这部分已经在当前的 Chemprop 数据准备脚本中落实。

## 后续开发成包时的建议模块划分

如果后续把流程开发成包，建议至少拆成这些功能模块：

1. `resolve_case(cas)`
2. `fetch_targets(case)`
3. `standardize_targets(case)`
4. `build_venn(case)`
5. `run_enrichment(case)`
6. `extract_kegg_key_genes(case)`
7. `prepare_chemprop(case)`
8. `train_chemprop(case)`
9. `evaluate_chemprop(case)`
10. `rank_targets(case)`
11. `plot_figure1e(case)`
12. `assemble_figure1(case)`

每个函数都应该只接受一个 case 对象或一个 case 目录，并且只返回这个 case 自己的结果。

## 每个 CAS 建议的最终交付物

对于一个 CAS case，建议最终输出包括：

- source summary
- `1.p_venn.pdf`
- `2.GO_KEGG.pdf`
- `3.KEGG_circos.pdf`
- `key_targets_from_kegg.csv`
- `chemprop_multitask.csv`
- `test_metrics_by_target.csv`
- `ranked_targets.csv`
- 过滤后的 `Figure 1e`
- 一张完整拼好的 `Figure 1`

这些结果应始终保留在同一个 `CAS_xxxxx_Date` 项目目录中，并按 `a-e` 子目录分类存放。

## 最终规则

后续开发时，必须把这条原则固定下来：

- 这不是一个多化合物混合报告系统
- 也不是一个多个分子共用一张靶标排序表的流程
- 而是一个**单 CAS 分析包**

这个包可以支持批量重复运行多个 CAS，但每一次运行都必须保持结果隔离。
