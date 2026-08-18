# 快照来源与口径

**本目录是研究仓库某一次生产状态的复制品，不是活的工作目录。**
研究仓库仍在 `03Claim/`；两边的 `src/` 与 `data/` 会随研究推进而分叉。
**引用本包的任何数字之前，先看下面的口径栏。**

---

## 快照标识

| 项 | 值 |
|---|---|
| 生成日期 | 2026-08-18 |
| 来源 | `F:\OneDrive\03researches\05SCC\03Claim` |
| 源仓库 commit | `47b4e55`（分支 `main`） |
| ⚠️ 工作区状态 | **有 73 个未提交改动**——本快照对应的是**工作区状态**，不是那个 commit |
| 目录指纹（路径+大小 SHA-256 前 16 位） | `5fa53ec23836055c` |
| 文件数 / 体积 | 106 文件 / 185.4 MB |
| 验证环境 | Julia 1.12.6 |

> **「有 73 个未提交改动」这一行是本表最重要的一行。** 它意味着仅凭 commit 号
> **无法**重建本快照。要精确复现，须以本目录本身为准，或先在研究仓库提交后重新生成。

## 自检值

```
compute_scco2(get_model(:national; climate=:fair, scenario="ssp245", end_year=2100);
              year=2030, last_year=2100)  ==  690.4507205341098
```

本快照已实测通过（独立于 `03Claim` 运行，逐位一致）。

---

## 本快照的模型口径

数值随口径变化，且**同一个模型不同口径的数字不可混引**。本快照是：

| 口径项 | 取值 |
|---|---|
| 情景 | **以实测 2024 年排放重基准后的 ssp245**（不是纯粹的 SSP2-4.5，见下）；脉冲年 2030；损害积分至 2100 |
| 排放输入 | `fair_data/rcmip_ssp*_..._rebased_obs.csv`：**2016–2024 用实测清单**（EDGAR-2025 / CEDS v_2025_03_18 / GFED5.1），2025–2500 按实测 2024 年水平整体重基准。原始 RCMIP 情景文件同目录留存但模型不读。见 `CLAIM_model_guide.md` §4.1.1 |
| 标准误 | year-block bootstrap，**块长 15**（非 Driscoll-Kraay 渐近） |
| 资本留存权重 $\rho_c$ | **由 PWT 11.0 实测资本存量 `rnna` 算出**（非 $(1-\delta)/(1+g)$ 代理） |
| 结构参数缺失填补 | 同洲之内、人均收入最近 3 国均值（非统一默认值） |
| Ramsey 折现 | FUND 口径：人均产出增长率、一阶折现核、逐国折现后跨国求和 |
| 排放脉冲 | FUND 惯例，1 吨 CO₂ 平摊到 10 年 |
| 损害核 | 单核 $\zeta_s=Ae^{-Bs}+Cs^{\nu}e^{-Ds}$，$\nu=3$ 为校准常数 |

### 到达本口径的变更

| 变更 | 日期 | SCC |
|---|---|---|
| 起点 | 2026-08-08 | 679.3368902160977 |
| $\alpha$ 缺失国改用同洲收入近邻填补 | 2026-08-08 | +0.733% → 684.3140 |
| Ramsey 折现因子期次对齐 FUND | 2026-08-08 | −2.482% → 667.3285 |
| $\rho$ 换成实测资本存量口径 | 2026-08-08 | +3.00% → 687.3172 |
| **排放输入的观测校正**（2016–2024 换实测、2025+ 重基准） | 2026-08-18 | +0.456% → **690.4507** |
| 全量去掉 `_bkw` 命名（纯改名，见下） | 2026-08-18 | **不变**（逐位 690.4507205341098） |

### 2026-08-18 的命名变更

模型已把 Bilal & Känzig 的损害量化方法**内化为一部分**，但只是一部分；到处挂 `bkw`
会让使用者误以为整个模型是他们研究的集成。故全量改名，**纯重命名、零数值影响**
（改后 SCC 逐位不变，全部测试保持绿色）：

| 旧 | 新 |
|---|---|
| `MimiCLAIM_bkw.jl` / `.ipynb` | `MimiCLAIM.jl` / `.ipynb` |
| `helper_bkw.jl` | `helper_damage.jl` |
| `components_damage/bkw_damage.jl` | `components_damage/claim_damage.jl` |
| 组件名 `:bkw_damage` | `:claim_damage` |
| `add_bkw_components!` / `add_bkw_shared_params!` / `connect_bkw_components!` / `setup_bkw_dimensions!` | `add_claim_*` / `connect_claim_*` / `setup_claim_*` |
| `update_bkw_damage_params!` | `update_claim_damage_params!` |
| `load_bkw_state_centers` | `load_damage_state_centers` |

**`:BK_global` / `:BK_region` 两个开关保留原名**：它们是论文对比开关，字面含义就是
"复现 Bilal & Känzig 的那个设定"，改名反而失真。

那次排放校正同时把模型 2024 年的温室气体浓度拉近观测：CH₄ 误差由 −21.33 ppb
变成 +10.85 ppb，三气体相对误差合计由 1.5124% 降到 1.0183%。CO₂ 略微变差
（+0.76 → +0.95 ppm），残差在 FaIR 的碳循环参数上，不在排放清单上。

---

## 复制了什么、没复制什么

**复制**（逐位复制，未做任何修改）：

- `src/` 全部 21 个文件
- `data/shared_parameters/`、`data/unshared_parameters/`、`data/fair_data/`、`data/shapefile/`
- `CLAIM_model_guide.md`

**未复制**：

| 未复制 | 原因 |
|---|---|
| `data/mesmer_data/`（32.8 GB） | 标定阶段输入，`get_model` 不读；含第三方 ESM 原始数据 |
| `scripts/`、`test/`、`results/` | 研究仓库内容 |
| `open_issues.md`、`coupling_plan.md` 等 | 研究过程文档 |

## 两档发行

本目录含**完整包**（185.4 MB）。分发时分两档，界线是**一个文件**：

| 档 | 体积 | 内容 |
|---|---|---|
| 基础包 | **75.5 MB** | 除下面那一个文件之外的全部 |
| 可选组件 | +109.9 MB | `data/unshared_parameters/climate_grid_innov_chol.csv` |

```
climate_grid_innov_chol.csv
  SHA-256  4265616099C20059E13CE981677F15784A3751EE6F9305EAD41CD3F0EE3CA71F
  字节数    115,241,260
```

`06CLAIM/.gitignore` 已把它排除，所以 `git add .` 得到的就是基础包。
要随包分发它，用 Git LFS 或单独下载——**不要直接 `git add`**。

### 这条界线是实测出来的，不是按文件名猜的

| 测试 | 结果 |
|---|---|
| 删掉该文件后跑 `get_model` + `run` + `compute_scco2` | ✅ SCC = `690.4507205341098`，与完整包**逐位一致**（2026-08-18 重新实测） |
| 删掉后调用 `cached_mesmer_grid`（外层 MC 的入口） | ❌ `missing MESMER grid parameter file: ...climate_grid_innov_chol.csv` |

调用链：`load_innov_cholesky` → `load_mesmer_grid` → `cached_mesmer_grid`，
两个终点分别是 `helper_mc.jl` 的 `run_scc_mc`（且仅当 `sample_mesmer=true`）
与 `regenerate_national_mesmer_inputs`。**`get_model` 不在链上。**

缺它不会静默降级：外层 MC 直接报错并指名缺哪个文件。

---

## 如何重新生成本快照

```powershell
python 03Claim/scripts/sync_06claim.py          # 在研究仓库根目录运行
```

该脚本会重新复制并**重新做自检**（跑一次 SCC 与本文件记录的值比对），
不一致就报错而不是静默出一个过期的包。**不要手工复制**——
手工复制没有自检，而这个仓库已经多次栽在「产物看着新、其实是几个口径之前的」上。
