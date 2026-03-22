# 故障排查

## 安装问题

### `conda env create -f environment.yml` 很慢

这个环境同时混合了：

- Python
- PyTorch
- RDKit
- Node.js
- R
- Bioconductor

这类组合本来就会让依赖求解变慢，这是正常现象。

推荐做法：

1. 先创建 `miptd` 环境。
2. 激活 `miptd`。
3. 后续所有补包操作尽量在 `miptd` 环境里用自带的 `mamba` 完成。

```bash
conda env create -f environment.yml
conda activate miptd
```

### R 里缺少 `legendry`

`legendry` 目前没有放进 conda 环境定义里。

请在 `miptd` 环境中单独安装：

```bash
Rscript -e "install.packages('legendry', repos='https://cloud.r-project.org')"
```

### 找不到 `MIPTD` 命令

通常需要确认两件事：

1. 当前已经激活 `miptd`
2. 当前包已经通过可编辑模式安装

```bash
conda activate miptd
pip install -e .
```

然后检查：

```bash
MIPTD --help
MIPTD-validate --help
```

## 运行问题

### `SwissTargetPrediction` 超时

这是最常见的外部源故障之一。

当前包的处理方式是：

- 在 `status.json` 里把 Swiss 标记为失败
- 自动写入一个空的占位文件
- 只要其余来源仍然能提供足够证据，流程就继续

优先检查：

- `run.log`
- `status.json`
- `a_source_collection/swiss_fetch/`

### 某一个来源失败了，但流程还在继续

这是预期行为。

当前包采用的是“显式降级”策略：

- 单个来源失败不会立刻终止整条流程
- 失败会记录到 `status.json`
- 只要剩余证据还足够，后续 `b/c/d/e` 就会继续

### 流程在 Step A 后停止

如果四个来源标准化之后，非空的人类靶标集合少于 `2` 个，流程会主动停止。

这是故意设计的，目的是避免在证据明显不足时继续生成看起来完整但实际不可靠的结果。

优先检查：

- `status.json`
- `a_source_collection/source_summary.json`

## 资源问题

### 找不到 `idmapping.tsv` 或 `ChEMBL_target_catalog.csv`

当前包要求这两份正式资源位于：

- [`src/miptd/resources/idmapping.tsv`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/src/miptd/resources/idmapping.tsv)
- [`src/miptd/resources/ChEMBL_target_catalog.csv`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/src/miptd/resources/ChEMBL_target_catalog.csv)

如果资源缺失，可以重新生成：

```bash
Rscript scripts/build_idmapping_tsv.R --output src/miptd/resources/idmapping.tsv
python3 scripts/build_chembl_target_catalog.py \
  --output-csv src/miptd/resources/ChEMBL_target_catalog.csv \
  --summary-json src/miptd/resources/ChEMBL_target_catalog.summary.json
```

## 校验问题

### `MIPTD-validate` 失败

校验器会检查：

- 五个步骤目录是否存在
- 关键结果文件是否存在
- `inference_template.csv` 是否只包含一个分子
- `prediction_summary.json` 是否只包含一个分子
- 各处 `CAS` 是否与 `case_manifest.json` 一致

如果失败，优先检查：

- `case_manifest.json`
- `e_chemprop/chemprop_data/inference_template.csv`
- `e_chemprop/chemprop_model/analysis/prediction_summary.json`

## 推荐优先查看的调试文件

遇到大多数问题时，优先看这几个文件：

- `run.log`
- `status.json`
- `case_manifest.json`
- `a_source_collection/source_summary.json`

