# MIPTD

`MIPTD` 的全称是 `Multi-source Integrative Progressive Target Discovery`。

`MIPTD` 是一个面向单个 `CAS` 分子的、在线优先的分析包，用于复现 Figure 1 的靶标发现流程。  
该包从一个真实 `CAS` 编号出发，在线解析化合物身份信息，查询四个靶标来源，并在独立的 case 目录中完成完整的 `a -> b -> c -> d -> e` 分析链。

## 项目标题

**MIPTD: Multi-source Integrative Progressive Target Discovery**

## 项目介绍

这个项目的目标，是把一个由 `CAS` 编号标识的单分子，转换成一条可重复、可命令行执行、结果结构清晰的靶标发现流程。

整条流程包括：

1. 通过 `PubChem` 在线解析化合物身份信息
2. 获取化合物名称、分子式、`SMILES` 和 `InChIKey`
3. 查询四个靶标来源：
   - `SwissTargetPrediction`
   - `SEA`
   - `ChEMBL`
   - `PPB2`
4. 将人类靶标统一到基因层级
5. 构建 Figure 1 的完整分析链：
   - `a`：来源收集流程
   - `b`：Venn 图
   - `c`：GO/KEGG 富集与疾病关键词筛选
   - `d`：KEGG circlize 和弦图
   - `e`：Chemprop multi-task 靶标优先级分析

这个包的基本单位是：

- 一个 `CAS`
- 一个带日期的 case 目录
- 一套独立结果

## 技术栈

项目采用 Python 和 R 的混合技术栈。

### Python

- `Python 3.11`
- `requests`
- `beautifulsoup4`
- `playwright`
- `chemprop`
- `matplotlib`（在 `environment.yml` 中以 `matplotlib-base` 形式声明，被 `scripts/plot_figure1e.py` 使用）
- `pillow`（被 `scripts/build_figure1.py` 用于图像合成）
- Python 标准库中的 API、HTML 解析和流程编排模块

### R

- `r-base`
- `tidyverse`
- `readxl`
- `writexl`
- `ggvenn`
- `clusterProfiler`
- `org.Hs.eg.db`
- `ggnewscale`
- `circlize`
- `ComplexHeatmap`
- `legendry`（由 `scripts/install_miptd.sh` 从 CRAN 安装，目前 conda-forge / bioconda 上还没有）

### 外部在线来源

- `PubChem PUG REST`
- `SwissTargetPrediction`
- `SEA`
- `ChEMBL API`
- `PPB2`

## 仓库结构

- [`src/miptd`](src/miptd)
  Python 包源码
- [`scripts`](scripts)
  分析、抓取和资源生成脚本
- [`src/miptd/resources`](src/miptd/resources)
  正式静态运行资源文件
- [`doc`](doc)
  流程与包文档

## 内置资源

当前包内置两份正式资源文件：

- [`src/miptd/resources/idmapping.tsv`](src/miptd/resources/idmapping.tsv)
  人类基因的 `UNIPROT -> SYMBOL` 映射表
- [`src/miptd/resources/ChEMBL_target_catalog.csv`](src/miptd/resources/ChEMBL_target_catalog.csv)
  人类 `ChEMBL target` 目录表，用于靶标标准化与 Chemprop 任务构建

这两份文件属于运行时正式资源，发布和分发时必须和代码一起保留。

对应的资源生成脚本也保存在 [`scripts`](scripts) 下：

- [`scripts/build_idmapping_tsv.R`](scripts/build_idmapping_tsv.R)
- [`scripts/build_chembl_target_catalog.py`](scripts/build_chembl_target_catalog.py)

## 安装方法

### Step 1：clone 仓库

```bash
mkdir -p ~/Github_repos
cd ~/Github_repos
git clone https://github.com/Jiawang1209/Multi-source-Integrative-Progressive-Target-Discovery-Framework.git
cd Multi-source-Integrative-Progressive-Target-Discovery-Framework
```

### Step 2：配置 conda 环境（A / B / C 三选一）

下面三种方式都必须**在 clone 目录里**执行（也就是你刚 `cd` 进去的那个目录）。

#### 方式 A — 推荐安装脚本

```bash
bash scripts/install_miptd.sh
```

如果要自定义环境名：

```bash
bash scripts/install_miptd.sh my_miptd_env
```

#### 方式 B — 手动安装

```bash
conda create -n miptd -y python=3.11 pip
conda install -n miptd -y -c conda-forge mamba
conda run -n miptd conda install -y -c pytorch -c bioconda -c conda-forge \
  nodejs pytorch rdkit r-base r-tidyverse r-readxl r-writexl r-ggvenn \
  r-jsonlite bioconductor-clusterprofiler bioconductor-org.hs.eg.db \
  r-ggnewscale r-circlize bioconductor-complexheatmap \
  matplotlib-base pillow
conda run -n miptd Rscript -e "install.packages('legendry', repos='https://cloud.r-project.org')"
conda run -n miptd pip install -e . chemprop==2.2.2
```

如果需要跑 `SwissTargetPrediction` 的浏览器自动化流程，需要确保本机浏览器依赖可用——抓取脚本走的是 Playwright + Google Chrome。
当前脚本使用 `http://www.swisstargetprediction.ch/` 入口，并会显式触发页面的 SMILES 输入校验，以适配 SwissTargetPrediction 在 2026-05 前后的页面更新。

#### 方式 C — 全锁定安装

```bash
conda env create -f environment.lock.yml
conda activate miptd
Rscript -e "install.packages('legendry', repos='https://cloud.r-project.org')"
pip install -e .
```

## 使用方法

### 在哪里跑 MIPTD

**MIPTD 不在源码 clone 目录里跑**。请在你为这次分析挑的任意工作目录里执行——通常是仓库之外某个"按项目/按 CAS 起名的工作目录"。pipeline 会在 `--output-root`（默认就是当前目录）下面建一个 `CAS_<编号>_<日期>/` 子目录写产物，clone 目录本身完全不参与产物写入。

具体模式：

```bash
mkdir -p ~/Desktop/analyses/paclitaxel
cd ~/Desktop/analyses/paclitaxel
conda activate miptd
MIPTD --cas 33069-62-4 --disease-keywords "cancer"
# 产物会写到 ~/Desktop/analyses/paclitaxel/CAS_33069-62-4_<今日日期>/
```

`~/Desktop/analyses/paclitaxel/` 可以随时清理或归档，MIPTD 安装本身不受影响——它住在 `~/Github_repos/.../` 下的 clone 里。

### 运行完整流程

```bash
MIPTD \
  --cas 491-71-4 \
  --disase-keywords "NAFLD,liver,hepatic,steatosis"
```

主要参数：

- `--cas`
  输入的 CAS 编号
- `--disease-keywords` / `--disase-keywords`
  逗号分隔的疾病关键词
- `--output-root`
  输出根目录
- `--run-date`
  可选的日期后缀
- `--project-root`
  可选的项目根目录
- `--dry-run`
  只打印计划流程，不真正执行

### 基于关键词的过滤

`Step C` 和 `Step D` 中的过滤，是由你传入的关键词驱动的，例如：

```bash
--disease-keywords "NAFLD,liver,hepatic"
```

程序会先对靶标并集执行完整的 GO/KEGG 富集，再把逗号分隔的关键词组合成一个正则式样的匹配规则，最后用这个规则去匹配富集结果中的 `Description` 字段。

具体逻辑是：

1. 先对并集靶标做完整 GO/KEGG 富集
2. 根据 `--disease-keywords` 构造关键词匹配规则
3. 用这个规则匹配 GO/KEGG 的 `Description`
4. 优先保留匹配到的条目用于生成 `c` 图
5. 再从筛选后的 KEGG 条目中提取关键基因用于生成 `d` 图

回退规则：

- 如果用户提供的关键词一个都没有匹配到，程序不会直接输出空图，而是回退到按显著性和计数排序的 top 条目，以保证 `c` 和 `d` 仍然可以生成

### 在线来源匹配规则

Step A 会分别查询 `SwissTargetPrediction`、`SEA`、`ChEMBL` 和 `PPB2`。这四个在线来源的页面和 API 会变化，因此当前包对关键输入做了显式约束：

- `SwissTargetPrediction`
  使用 Playwright + Google Chrome 自动提交 SMILES；脚本走 HTTP 入口，并在提交前确认页面已允许预测。
- `SEA`
  使用 SEA 当前搜索页的 `rdkit_ecfp` 指纹类型提交自定义 query；化合物标签会自动转换为无空格 ID，避免 SEA 拒收含空格的 compound id。
- `ChEMBL`
  优先用 PubChem 解析出的 `InChIKey` 和 `SMILES` 做精确匹配。如果传入了这些精确标识但 ChEMBL 没有命中，脚本会把 ChEMBL 作为失败来源处理，而不会退回到名称搜索的第一条结果，避免把别的分子错当成当前 CAS。
- `PPB2`
  使用 SMILES 查询 PPB2，再把返回的 ChEMBL target 映射到人类 UniProt/gene symbol；无法映射到 gene symbol 的行不会贡献到后续 Venn 计数。

### 校验 case 目录

```bash
MIPTD-validate \
  --case-dir CAS_491-71-4_2026-03-21
```

## 测试

本地单元测试与轻量集成测试可通过以下命令运行：

```bash
python3 -m unittest discover -s tests -v
```

## 输出结果

每次运行都会创建一个独立的 case 目录：

```text
CAS_491-71-4_YYYY-MM-DD/
  compound_identity.json
  compound_identity.csv
  status.json
  run.log
  a_source_collection/
  b_venn/
  c_go_kegg/
  d_kegg_circlize/
  e_chemprop/
```

这套结构对应一个硬约束：

- 一个 `CAS`
- 一个 case 目录
- 一套结果

## 失败处理策略

当前包对在线来源采用显式降级策略：

- 单个数据库抓取失败时，不会立刻中止整条流程
- 失败状态会写入 `status.json`
- 程序会写入一个空的占位来源文件，保证 case 目录结构完整
- 只有在标准化后仍然至少有 `2` 个数据库提供非空人类靶标集合时，流程才会继续
- 如果非空来源少于 `2` 个，流程会明确报错并停止

如果某次运行已经失败，失败来源目录里可能已经写入空的占位文件。修复脚本后请换一个 `--run-date`，或者删除旧的 `CAS_<编号>_<日期>/` case 目录后再重新运行，避免复用旧占位文件。

## 文档入口

详细文档位于 [`doc`](doc)：

- [`doc/README.md`](doc/README.md)
- [`doc/包功能与使用说明.md`](doc/包功能与使用说明.md)
- [`doc/单CAS_Figure1分析流程说明.md`](doc/单CAS_Figure1分析流程说明.md)
- [`doc/Single_CAS_Figure1_Pipeline.md`](doc/Single_CAS_Figure1_Pipeline.md)
- [`TROUBLESHOOTING.md`](TROUBLESHOOTING.md)
- [`TROUBLESHOOTING_CN.md`](TROUBLESHOOTING_CN.md)

## 当前状态

当前包已经具备：

- 单 `CAS` 工作流设计
- 在线化合物身份解析
- 四源靶标收集
- 结构化的 `a -> e` 流程编排
- 内置 `idmapping.tsv` 与 `ChEMBL target catalog`
- case 校验能力
- 流程日志能力

## 发布说明

在打包或发布前，至少应检查：

- 保留 `src/miptd/resources/` 目录
- 保留 `src/miptd/resources/idmapping.tsv`
- 保留 `src/miptd/resources/ChEMBL_target_catalog.csv`
- 确认 `MIPTD --help` 正常
- 确认 `MIPTD-validate --help` 正常
- 确认可以通过 [`environment.yml`](environment.yml) 成功创建环境


