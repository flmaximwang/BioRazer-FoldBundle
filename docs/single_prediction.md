# SinglePrediction

## 概述

`SinglePrediction` 是 BioRazer-FoldBundle 中用于表示“单个设计条目的一组预测结果”的基础分析类。一个 `SinglePrediction` 对象通常对应一个 marker 目录；该目录下可以包含多个 `seed` 和多个 `sample`，每个 `(seed, sample)` 组合对应一份结构文件和一份元数据文件。

这个类的职责不是直接执行预测，而是对已经生成并整理好的预测结果进行统一访问、汇总和后处理。它提供了以下核心能力：

- 检查预测目录是否已经整理为统一格式
- 收集所有样本到 `pandas.DataFrame`
- 读取结构、序列和 metadata
- 提取简单标量指标到样本表
- 生成 PAE 图
- 计算跨链 `ipSAE`
- 将每个样本展开为单行结果，便于后续汇总

从代码结构上看，`SinglePrediction` 是一个基础类，供不同预测后端的分析器继承。具体的 `format_output()` 一般由具体子类实现，例如 AlphaFold3 本地或服务端分析器。

## 类层次

`SinglePrediction` 的继承链如下：

```text
Entry
  -> SequenceEntry
    -> SinglePrediction
      -> 具体模型分析类，例如 SingleAF3Local / SingleAF3Server
```

因此它除了自身定义的方法外，也继承了若干通用接口：

- `from_dir()`：从目录构建对象
- `from_samples()`：从 DataFrame 构建对象
- `get_dir_path()`：返回当前 marker 对应的目录路径
- `get_samples()` / `set_samples()`：读写样本表
- `get_marker()`：返回 marker 名称

## 数据模型

### 目录级别

一个 `SinglePrediction` 对象对应一个 entry 目录，也就是：

```text
parent_dir/
  marker/
    ...
```

其中：

- `marker` 是条目标识，一般就是目录名
- `parent_dir` 是该目录的上级目录

对象内部通过以下几个核心属性描述状态：

- `marker: str`
- `parent_dir: str | None`
- `samples: pandas.DataFrame | None`
- `sequence: dict[str] | str | None`

### 样本级别

在格式化后的目录中，每个预测样本必须满足如下结构：

```text
marker/
  seed-1_sample-1/
    marker_seed-1_sample-1_model.cif
    marker_seed-1_sample-1_model.pdb
    marker_seed-1_sample-1_data.pickle
  seed-1_sample-2/
    ...
  seed-2_sample-1/
    ...
```

`collect_samples()` 会扫描这些目录并生成一个标准样本表，默认包含以下列：

- `seed`
- `sample`
- `cif`
- `pdb`
- `pickle`

之后其他方法会在这个 DataFrame 上继续追加分析列，例如：

- 自定义 metric 列
- `ipSAE`

## 对象创建

### 方式一：从目录创建

最常见的用法是通过继承类的 `from_dir()` 创建对象：

```python
from biorazer_fold_bundle.alphafold3.analysis import SingleAF3Local

pred = SingleAF3Local.from_dir("path/to/formatted/marker")
```

`from_dir()` 的行为来自父类 `Entry`：

1. 根据目录名设置 `marker`
2. 根据父目录设置 `parent_dir`
3. 若 `formatted=True`，自动执行：
   - `collect_samples()`
   - `get_basic_info()`

这意味着对于已经整理好的结果目录，通常可以一步完成初始化。

### 方式二：从样本表创建

如果你已经有整理好的样本表，也可以使用：

```python
pred = SomePredictionClass.from_samples(samples=df, marker="my_design")
```

这种方式适合：

- 从外部统计表恢复对象
- 手工构造结果集
- 在批量处理后重建单个条目

但要注意，此时对象可能没有真实目录，因此依赖文件路径的方法必须保证 DataFrame 中已有正确的文件列。

## 推荐工作流

对一个已经格式化好的预测条目，推荐的调用顺序通常如下：

```python
pred = SomePredictionClass.from_dir("path/to/marker")
samples = pred.get_samples()
keys = pred.get_metadata_keys()
pred.extract_simple_metrics(...)
structure = pred.get_structure(seed=1, sample=1)
fig, ax = pred.report_pae(seed=1, sample=1)
```

如果目录还没有被整理成统一格式，则应先调用具体子类实现的 `format_output()`，再从整理后的目录重新创建对象。

## 方法说明

## is_formatted

```python
is_formatted(self) -> bool
```

检查当前目录是否满足标准格式。

判断规则：

1. 根目录下至少存在一个 `seed-*_sample-*` 子目录
2. 每个样本目录都必须是目录而不是普通文件
3. 每个样本目录中都必须同时存在：
   - `*_model.cif`
   - `*_model.pdb`
   - `*_data.pickle`

典型用途：

- 在 `collect_samples()` 前做前置检查
- 在具体子类实现 `format_output()` 后验证结果
- 在批处理前剔除不完整条目

返回值：

- `True`：目录符合约定
- `False`：目录缺失必要文件或命名不符合要求

## collect_samples

```python
collect_samples(self, overwrite=False)
```

扫描当前 entry 目录下所有 `seed-*_sample-*` 目录，并构建 `self.samples`。

主要行为：

1. 调用 `is_formatted()` 验证目录
2. 解析每个样本目录名中的 `seed` 和 `sample`
3. 拼接出对应的 `cif`、`pdb`、`pickle` 文件路径
4. 生成并按 `seed`, `sample` 升序排序 DataFrame

注意事项：

- 当前实现中的 `overwrite` 参数没有实际被使用
- 如果目录未格式化，会抛出 `ValueError`
- 该方法会直接更新对象状态，即写入 `self.samples`

示例：

```python
pred.collect_samples()
print(pred.samples.columns)
```

## get_basic_info

```python
get_basic_info(self)
```

从第一条样本记录对应的 PDB 文件中读取序列，并写入 `self.sequence`。

当前实现使用：

```python
PDB2SEQ(self.get_samples().iloc[0]["pdb"]).read()
```

需要注意的地方：

- 该方法依赖 `self.samples` 已经存在
- 它默认用第一条样本代表整个条目的序列
- 如果不同样本的序列并不一致，这种做法不会检测冲突

## flatten_samples

```python
flatten_samples(self) -> pandas.DataFrame
```

把多行样本表压平成单行，便于：

- 合并到更大的统计表
- 导出宽表结果
- 与设计级别信息做 join

输出规则：

- 新表固定包含一列 `Marker`
- 原样本表中除 `seed` 和 `sample` 以外的所有列，都会被展开成：

```text
{column}_seed-{seed}_sample-{sample}
```

例如：

- `ipTM_seed-1_sample-1`
- `pdb_seed-2_sample-3`

该方法返回新 DataFrame，不会原地覆盖 `self.samples`。

## extract_simple_metrics

```python
extract_simple_metrics(
    self,
    metric_key_list: list[str],
    metric_id_list: list[tuple],
    metric_label_list: list[str],
) -> pandas.DataFrame
```

这是 `SinglePrediction` 最常用的方法之一，用于从每个样本的 metadata 中提取标量指标，并将结果追加到 `self.samples`。

### 参数含义

- `metric_key_list`
  - metadata 中的一级 key 列表
- `metric_id_list`
  - 对应每个 key 的索引路径
  - 当 metadata[key] 是嵌套 list 时，按这个 tuple 逐层取值
- `metric_label_list`
  - 写入样本表时使用的列名

### 工作方式

对每个样本：

1. 通过 `get_metadata(seed, sample, metric_key_list)` 读取指定 key
2. 若 `metadata[key]` 不是 list，则直接写入
3. 若是 list，则按 `metric_id` 逐层索引
4. 把最终值写到 `samples.loc[sample_i, metric_label]`

### 示例

提取一个标量字段：

```python
pred.extract_simple_metrics(
    metric_key_list=["ptm"],
    metric_id_list=[()],
    metric_label_list=["pTM"],
)
```

提取矩阵中的某个元素：

```python
pred.extract_simple_metrics(
    metric_key_list=["chain_pair_iptm"],
    metric_id_list=[(0, 1)],
    metric_label_list=["A_B_ipTM"],
)
```

### 注意事项

- 必须先执行 `collect_samples()`，否则会报错
- `metric_key_list`、`metric_id_list`、`metric_label_list` 应一一对应
- 当前实现仅在 `metadata[key]` 是 Python `list` 时走索引逻辑；若是 `numpy.ndarray`，会被视为非 list，需要调用方自行处理或先确认数据类型

## get_metadata_keys

```python
get_metadata_keys(self) -> list[str]
```

读取第一条样本对应的 pickle 文件，并返回 metadata 的所有 key。

适合用于：

- 探查一个预测后端导出了哪些字段
- 在调用 `extract_simple_metrics()` 前确认字段名

注意：

- 该方法依赖 `self.samples` 已存在
- 它默认第一条样本的 key 集合代表全部样本

## get_metadata

```python
get_metadata(self, seed, sample, metadata_keys) -> dict
```

读取指定样本的 pickle metadata，并返回一个子字典，只包含请求的 key。

参数：

- `seed`: 样本 seed
- `sample`: 样本编号
- `metadata_keys`: 需要返回的 key 列表

返回值示例：

```python
{
    "ptm": 0.81,
    "pae": ...,
}
```

错误行为：

- 如果指定 key 不存在，会抛出 `ValueError`
- 如果 `(seed, sample)` 没有匹配记录，当前实现会在后续 `iloc[0]` 处失败

因此更稳妥的使用方式是先确认该样本存在，或者通过 `get_sample()` 检查。

## get_structure

```python
get_structure(self, seed, sample, file_type="cif")
```

读取指定样本的结构对象。

支持：

- `file_type="cif"`：通过 `CIF2STRUCT(...).read()` 读取
- `file_type="pdb"`：通过 `PDB2STRUCT(...).read()` 读取

适合用于：

- 后续结构分析
- 几何测量
- 与 confidence 信息联动分析

如果 `file_type` 不是 `cif` 或 `pdb`，会抛出 `ValueError`。

## get_sample

```python
get_sample(self, seed, sample) -> pandas.DataFrame
```

返回匹配指定 `(seed, sample)` 的样本表子集，并重置索引。

这里返回的是 DataFrame，而不是单行 Series，因此后续通常会写成：

```python
sample_df = pred.get_sample(1, 1)
sample_row = sample_df.iloc[0]
```

这个接口主要用于：

- 获取某个样本的路径和分析列
- 在绘图或导出前定位特定样本

## report_atom_plddts

```python
report_atom_plddts(self, seed, sample)
```

读取指定样本 metadata 中的 `atom_plddts`，并将其写入对应 CIF 结构的 `b_factor` 注释字段，然后输出一个带 `_plddt` 后缀的新 CIF 文件。

处理流程：

1. 从 pickle 中读取 `atom_plddts`
2. 读取原始 CIF
3. 将每个原子的 pLDDT 写入 `b_factor`
4. 把新结构写到：

```text
*_model_plddt.cif
```

这一方法的意义在于：

- 让结构查看器直接按 B-factor 渲染 pLDDT
- 便于用现有结构分析工具查看置信度分布

## report_residue_plddts

```python
report_residue_plddts(self, seed, sample)
```

该接口目前尚未实现，当前代码中是 `pass`。

如果后续实现，通常会期望它：

- 把 residue 级别置信度映射到结构
- 或导出 residue 级别的统计表

在当前版本中，不应依赖这个方法。

## report_pae

```python
report_pae(
    self,
    seed,
    sample,
    figsize=(5, 5),
    vmin=0,
    vmax=31.75,
    cmap="coolwarm",
) -> tuple[Figure, Axes]
```

绘制指定样本的 PAE 矩阵图。

主要行为：

1. 通过 `get_sample()` 定位样本
2. 从 pickle 读取 `pae`
3. 使用 `matplotlib.pyplot.imshow()` 画热图
4. 添加 colorbar，并将标签设为 `PAE (Å)`

返回值：

- `fig`
- `ax`

这让调用方可以继续自定义：

- 标题
- 坐标轴标签
- 保存图片

示例：

```python
fig, ax = pred.report_pae(seed=1, sample=1)
ax.set_title("seed-1 sample-1")
fig.savefig("pae.png", dpi=300, bbox_inches="tight")
```

## calc_ipsae

```python
calc_ipsae(
    self,
    chain1,
    chain2,
    pae_cutoff=10.0,
    dist_cutoff=10.0,
    pair_type="protein",
)
```

为当前对象中的所有样本批量计算 `ipSAE`，并把结果写入 `self.samples["ipSAE"]`。

工作流程：

1. 先在样本表中新建 `ipSAE` 列并初始化为 `0.0`
2. 遍历每个样本
3. 读取该样本的 `pae`
4. 读取该样本的 CIF 结构
5. 调用 `advanced_metrics.ipsae.calc_ipsae(...)`
6. 将结果写回样本表

参数说明：

- `chain1`, `chain2`
  - 要计算界面分数的两条链
- `pae_cutoff`
  - PAE 截断阈值
- `dist_cutoff`
  - 接触距离阈值
- `pair_type`
  - 链对类型，默认是 `protein`

适用场景：

- 蛋白-蛋白界面质量评估
- 多样本排序前的统一打分
- 批量筛选高可信界面模型

## 与 format_output 的关系

`SinglePrediction` 本身依赖“格式化后的结果目录”，但并不直接提供统一可用的 `format_output()` 实现。这个职责通常由具体预测器子类完成。

例如：

- AlphaFold3 本地结果分析器会把原始输出整理成标准目录，并将 JSON confidence 合并为 pickle metadata
- AlphaFold3 服务端结果分析器会将多个 model 映射为 `seed-x_sample-y` 结构

因此你可以把 `SinglePrediction` 理解为：

- 一个统一的分析接口层
- 一个要求输入目录满足约定的数据访问基类

## 最小示例

下面是一个相对完整的典型示例：

```python
from biorazer_fold_bundle.alphafold3.analysis import SingleAF3Local

pred = SingleAF3Local.from_dir("formatted_results/design_001")

print(pred.marker)
print(pred.sequence)
print(pred.get_metadata_keys())

pred.extract_simple_metrics(
    metric_key_list=["ptm", "iptm"],
    metric_id_list=[(), ()],
    metric_label_list=["pTM", "ipTM"],
)

pred.calc_ipsae(chain1="A", chain2="B")

sample = pred.get_sample(seed=1, sample=1)
structure = pred.get_structure(seed=1, sample=1, file_type="cif")
fig, ax = pred.report_pae(seed=1, sample=1)

summary = pred.flatten_samples()
print(summary.columns)
```

## 常见注意事项

### 1. 先有格式，再谈分析

`collect_samples()` 依赖严格的目录命名和文件命名。如果原始预测结果还没有整理成标准格式，应先调用具体子类的 `format_output()`。

### 2. 多数接口依赖 self.samples

以下方法都隐含要求 `self.samples` 已存在：

- `get_basic_info()`
- `flatten_samples()`
- `extract_simple_metrics()`
- `get_metadata_keys()`
- `get_metadata()`
- `get_structure()`
- `get_sample()`
- `report_atom_plddts()`
- `report_pae()`
- `calc_ipsae()`

如果对象不是通过 `from_dir(..., formatted=True)` 创建，就需要手动先调用 `collect_samples()`。

### 3. metadata 的结构需要调用方了解

`extract_simple_metrics()` 不会自动推断复杂嵌套结构的意义。使用前最好先：

```python
keys = pred.get_metadata_keys()
meta = pred.get_metadata(1, 1, keys)
```

看清楚字段的真实结构，再决定 `metric_id_list` 如何填写。

### 4. sequence 读取基于第一条样本

若不同样本可能携带不同序列，当前实现不会做一致性检查。对于这类特殊场景，建议自行额外验证。

## 相关类

通常你不会直接使用一个“裸”的 `SinglePrediction`，而是使用它的具体子类，例如：

- `SingleAF3Local`
- `SingleAF3Server`
- 其他具体预测器对应的单条目分析器

而批量场景则对应：

- `BatchPrediction`
- 各后端自己的 batch 分析类

二者关系可以概括为：

- `SinglePrediction` 负责一个 marker
- `BatchPrediction` 负责多个 marker，并把每个 entry 的结果组织成更高层级的数据集

## 总结

`SinglePrediction` 的核心价值在于把不同预测后端输出的多样本结果统一成一套稳定的分析接口。只要具体子类能把原始输出整理成标准目录结构，后续就可以通过同样的方式去：

- 收集样本
- 读取结构与 metadata
- 提取 confidence 指标
- 绘制 PAE
- 计算界面分数
- 输出汇总表

这让上层分析代码可以尽量少关心具体预测器的原始输出细节，而把注意力集中在结果筛选和结构解释上。