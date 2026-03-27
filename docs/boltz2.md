# Boltz2

## 概述

这个模块目前的核心作用是：

- 用 Python 对象描述 Boltz2 的输入
- 批量生成 Boltz2 使用的 YAML 请求文件
- 给出基础运行命令模板

对应实现位于 [biorazer_fold_bundle/apps/boltz2/execution.py](biorazer_fold_bundle/apps/boltz2/execution.py)。

当前封装里最常用的几个类是：

- `Boltz2Sequence`：描述 protein / dna / rna / ligand 输入
- `Boltz2JobGenerator`：描述单个预测任务
- `Boltz2JobBatch`：批量写出 YAML 请求

## 单链预测示例

下面给出一个最小可用示例，目标是：

- 输入 1 条蛋白序列
- 只生成 1 个预测任务
- 不依赖外部 `.a3m` 文件
- 不额外跑 data pipeline

输入序列为：

```text
MKGDTKVINMLNKLLGLLLVLINTAFLAARMAKNMGDKLANDVLYHISINAMKMADKAIELILFLEGLPNLQDLGKLNIGSSGGSSINLMNLVLGLLLVLINQAFLIARMAKNLGDKLTNDIAYHISIEAMKNADAAIETILFMEGLPNLQDLGKLNI
```

这里把“1 个 repeat”按最直接的方式理解为：只提交 1 个 job，也就是只生成 1 份 YAML 配置。

## Python 示例

```python
from biorazer_fold_bundle.apps.boltz2.execution import (
    Boltz2Sequence,
    Boltz2JobGenerator,
    Boltz2JobBatch,
)


sequence = (
    "MKGDTKVINMLNKLLGLLLVLINTAFLAARMAKNMGDKLANDVLYHISINAMKMADKAIELILFLEGL"
    "PNLQDLGKLNIGSSGGSSINLMNLVLGLLLVLINQAFLIARMAKNLGDKLTNDIAYHISIEAMKNADAA"
    "IETILFMEGLPNLQDLGKLNI"
)

protein = Boltz2Sequence(
    entity_type="protein",
    id=["A"],
    sequence=sequence,
)

# 不依赖外部 MSA 文件，直接写入一个只有 query 的最小 MSA。
# 对这个封装来说，这是当前最直接的“不要额外准备 data pipeline 输入”的写法。
protein.set_no_msa()

job = Boltz2JobGenerator(
    name="single_repeat_1",
    sequences=[protein],
)

batch = Boltz2JobBatch(
    executions=[job],
)

batch.generate_requests("example_boltz2_jobs")
```

运行后会在目标目录下生成：

```text
example_boltz2_jobs/
  single_repeat_1.yaml
  single_repeat_1_A.a3m
```

其中：

- `single_repeat_1.yaml` 是 Boltz2 的配置文件
- `single_repeat_1_A.a3m` 是由 `set_no_msa()` 自动导出的 query-only MSA

## 生成的 YAML 示例

上面的代码会生成一个与下面等价的 YAML：

```yaml
sequences:
  - protein:
            id:
                - A
            sequence: MKGDTKVINMLNKLLGLLLVLINTAFLAARMAKNMGDKLANDVLYHISINAMKMADKAIELILFLEGLPNLQDLGKLNIGSSGGSSINLMNLVLGLLLVLINQAFLIARMAKNLGDKLTNDIAYHISIEAMKNADAAIETILFMEGLPNLQDLGKLNI
            msa: ./single_repeat_1_A.a3m
```

对应的 `.a3m` 文件内容为：

```text
>query
MKGDTKVINMLNKLLGLLLVLINTAFLAARMAKNMGDKLANDVLYHISINAMKMADKAIELILFLEGLPNLQDLGKLNIGSSGGSSINLMNLVLGLLLVLINQAFLIARMAKNLGDKLTNDIAYHISIEAMKNADAAIETILFMEGLPNLQDLGKLNI
```

## 运行命令

按照当前实现，命令模板为：

```bash
boltz2 run --config single_repeat_1.yaml --output output_dir
```

如果你已经进入请求文件所在目录，直接运行即可。

## 说明

### 关于“不需要 run data pipeline”

这个封装本身没有单独暴露一个类似“skip pipeline”的命令行开关；它做的是把输入组织成 YAML。

在当前 API 下，最稳妥的写法就是：

- 不提供外部 `.a3m` 文件路径
- 直接调用 `set_no_msa()` 生成 query-only MSA

这样你不需要先准备额外的比对结果文件，就可以把最小输入写出来。

### 关于“1 个 repeat”

当前这个模块里没有专门叫 `repeat` 的字段；最自然的映射方式就是：

- 只构造 1 个 `Boltz2JobGenerator`
- 只写出 1 个 YAML 请求

如果你后面希望的是“同一条序列做多次独立重复预测”，可以把 `executions` 改成多个 job，例如：

```python
jobs = [
    Boltz2JobGenerator(name=f"single_repeat_{i}", sequences=[protein])
    for i in range(1, 4)
]
batch = Boltz2JobBatch(executions=jobs)
```

那样会一次生成多份 YAML。
