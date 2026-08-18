# CLAIM 模型完整说明文档

**CLAIM: Climate-Linked Assessment of Impact Model —— 基于 Bilal & Känzig 损害量化方法的碳社会成本评估模型**

> 版本：2026-07-30
> 适用读者：希望学习、使用或扩展 CLAIM 模型的研究者
> 配套文件：`coupling_plan.md`（耦合技术方案）、`open_issues.md`（未解决问题清单）、`02regression/damage_inversion_plan.md`（反演方法论）、**`damage_uncertainty_guide.md`（损害不确定性的构建、传播与评估；含闭式反演与 bootstrap 接入的完整推导，面向初学者）**
> 主体语言为中文；变量名、文件名、组件名保留英文原名。

---

## 目录

1. [模型概述](#1-模型概述)
2. [总体架构与数据流](#2-总体架构与数据流)
3. [理论方法](#3-理论方法)
4. [Mimi 组件模块详解](#4-mimi-组件模块详解)
5. [文件逐一说明](#5-文件逐一说明)
6. [SCC 评估的实现路径](#6-scc-评估的实现路径)
7. [运行指南](#7-运行指南)
8. [当前状态与已知问题](#8-当前状态与已知问题)
9. [参考文献](#9-参考文献)

---

## 1. 模型概述

### 1.1 CLAIM 是什么

CLAIM 是一个基于 Julia **Mimi 框架**构建的集成评估模型（Integrated Assessment Model, IAM），用于评估气候变化的宏观经济损害与**碳社会成本**（Social Cost of Carbon, SCC）。它的完整因果链条为：

```
温室气体排放 → 大气浓度（FaIR v2.0 气体循环）→ 辐射强迫 → 全球温度
    → 区域/国家温度（MESMER 统计降尺度）→ 经济损害（损害函数）
    → 收入路径改变 → 折现加总 → SCC
```

### 1.2 本项目的核心创新：引入 Känzig 损害量化方法

传统 IAM（DICE、FUND、PAGE）的损害函数将**当期温度水平**映射到**当期 GDP 损失**，忽略了损害的动态传播。Bilal & Känzig（2026, QJE，下文简称 **B&K**）的核心贡献是用**局部投影（Local Projection, LP）**估计温度冲击对 GDP 的**动态脉冲响应（IRF）**，发现：

1. **全球（海洋含）温度冲击的损害远大于局部温度冲击**——1°C 全球温度冲击导致世界 GDP 峰值下降约 12%，是既有文献的 6 倍；
2. 损害具有**持续性**——冲击后损害持续累积多年而非即时体现；
3. 由此推出的 SCC 高达每吨 CO₂ 千美元量级，远超 EPA 官方估计。

本项目在 B&K 基础上做了四项扩展，然后把损害动态嵌入 CLAIM：

| 扩展 | 内容 |
|------|------|
| 状态依赖 | 借鉴 Nath et al. 思路，让温度响应系数 θ 依赖于国家气候态 barT 与初始收入 lny0，随升温动态演化 |
| 逐国损害核 | 185 国各自反演结构损害核 ζ（而非单一全球核） |
| Solow 资本动态 | 在本项目自有的"非平稳资本动态框架"（Solow 线性化）下重新回收 ζ，避免直接借用 B&K 在 Ramsey 动态下回收的参数造成核错配 |
| 不确定性传播 | LP 系数的 6×6 联合协方差 + **秩-6 因子化平滑** MVN 抽样 + 外层气候 Monte Carlo。**不是完整不确定性区间**：形状参数 $B/D$ 与校准常数 $
u$ 被冻结、气候×损害联合嵌套 MC 尚未实现，见 `open_issues.md` §6.2 R23/R24 与 P1-8 |

### 1.3 当前模型范围

当前 `03Claim` 已整理为**国家层级 B&K CLAIM 模型**。旧 FUND-like / 省级运行入口和 `impactaggregation` 增长率损害模块已移除，避免与新版水平 TFP 损害核混用。论文模型对比新增了 BHM 二次温度损害组件与 B&K 三开关接口，但默认仍为完整 B&K/CLAIM 机制。

| 维度 | 当前实现 |
|------|----------|
| 组件 | `claim_damage`（默认）/ `bhm_damage`（论文对比） |
| 主入口 | `src/MimiCLAIM.ipynb`（用户交互/SCC/出图，内联 `get_model`）+ `src/MimiCLAIM.jl`（测试/脚本库入口） |
| 函数形式 | B&K ζ_s 水平损害核卷积；BHM 二次绝对温度增长损害用于对比 |
| 空间分辨率 | 国家级 185 国（ISO3） |
| 时间结构 | TFP 核 11 年记忆 + Solow 资本无限记忆 |
| 不确定性 | 因子化 MVN + 嵌套 MC |

### 1.4 快速开始

**我只想理解模型**：读本文档 §1–3，重点抓住 FaIR → MESMER → damage → SCC 的数据流（§2.2 流程图）。

**我想跑 national 模式的 smoke test**（须用装有 Mimi 的 Julia；1.12 与 1.8.1 都可）：

```powershell
julia 03Claim/test/run_all_tests.jl
```

若要钉住某个 Julia 版本，用该版本可执行文件的完整路径代替上面的 `julia`，例如
`& "<你的 Julia 安装路径>\bin\julia.exe" 03Claim/test/run_all_tests.jl`。

**我想重新生成 B&K 输入数据**：按 §6 阶段 1 依次运行 `prepare_national_inputs.py` 与三个 `precompute_*.py`。

**我想计算 SCC**：`get_model(:national; climate=:fair, scenario="ssp245")` → `compute_scco2(m; year=2025)`，或打开 `src/MimiCLAIM.ipynb` 分块跑 SCC / 写 CSV / 出图。`scenario` 仅限 `ssp119`/`ssp126`/`ssp245`/`ssp370`/`ssp585`；默认损害起始 `bau_ref_year=2025`，积分至 2150。默认损害机制为 `damage_function=:BK_region, income_state=:dynamic, capital_response=true`。论文五方案对比可改用 `(:BHM_global,:none,false)`、`(:BK_global,:none,false)`、`(:BK_global,:none,true)`、`(:BK_region,:none,true)`、`(:BK_region,:dynamic,true)`。

---

## 2. 总体架构与数据流

### 2.1 三大目录

```
05SCC/
├── 01data/ + 01regression/     省级面板的早期探索管线（已被 02regression 取代）
├── 02regression/               计量估计端：LP 回归 + 损害核反演
│   ├── 01scripts/              01_merge_raw.py → 02_lp_regression.do
│   │                           → 03_invert_damage.py → 04_plot_slowB.py
│   ├── 02data/                 pwt110.xlsx, raw_merged.dta, analysis.dta,
│   │                           berkeley_full.dta, struct_params.csv,
│   │                           country_temperature_berkeley_grid_pw.csv
│   ├── 03table/                lp_main.csv, lp_imputed.csv,
│   │                           damage_params.csv, temp_irf.csv
│   └── 04figure/               lp_irf.png, fit_examples.png 等诊断图
└── 03Claim/                    模型端：IAM + SCC
    ├── src/                    Mimi 组件与主模型
    ├── scripts/                Python 预处理（回归输出 → CSV 参数表）
    ├── data/                   模型输入数据（FaIR、MESMER、shared/unshared 参数）
    └── test/                   Julia 单元测试
```

### 2.2 端到端数据流

```
[PWT 11.0 (`rgdpna/pop`) + Berkeley Earth + precipitation + macro controls]
        │  02regression/01scripts/01_merge_raw.py
        ▼
raw_merged.dta ──► 02_lp_regression.do（Stata）
        │   产出: lp_main.csv    （pooled θ 系数 + 21 个 VCV 元素）
        │         lp_imputed.csv （185 国 θ_c(h) + 国家级 SE，Driscoll–Kraay 渐近）
        │         temp_irf.csv   （温度自身 IRF φ^T）
        ▼
02d_bootstrap_lp.do（年份块 bootstrap，R=300，块长 10）
        │   产出: lp_bootstrap_bl10.csv（300 × 11 replicate）
        ▼
02e_lp_imputed_bootstrap.py
        │   产出: lp_imputed_boot.csv（θ 两列逐位不变，只换 SE 两列）
        │   ★ 2026-08-06 起这是**生产口径**，下游一律读它
        ▼
03_invert_damage.py（单核加权 NLS 反演，ν=3；权重 1/SE_boot²）
        │   产出: damage_params.csv（每国每渠道核类型 + 形状参数 B, D/ω）
        │         struct_params.csv（每国 α, δ, g）
        ▼
03Claim/scripts/precompute_*.py（预处理三件套）
        │   产出: data/unshared_parameters/damage_basis.csv
        │         data/unshared_parameters/damage_structural_params.csv
        │         data/unshared_parameters/damage_weights.csv
        │         data/unshared_parameters/damage_theta_base.csv
        │         data/unshared_parameters/damage_theta_draws.csv
        │         data/unshared_parameters/damage_cholesky.csv
        │   另: prepare_national_inputs.py → shared_parameters/ 与 unshared_parameters/ CSV
        ▼
[CMIP6 tas] ─► mesmer_data/export_grid_params.py（MESMER 训练，网格级 β/γ/AR 参数）
        │              │
        │   [GlobPOP + world_countries.shp]
        │              ▼
        │   scripts/build_mesmer_country_weights.py → climate_grid_weights.csv（185 × 5940）
        │              ▼
        │   src/helper/helper_mesmer.jl（网格模拟变率 → 人口加权聚合）
        │              产出: climate_ltcoef.csv / climate_ltinter.csv /
        │                    climate_region_lv{,_local}.csv / climate_global_lv.csv
        ▼
03Claim/src/MimiCLAIM.jl（get_model → run → compute_scco2）
        │
        ▼
BAU 损失路径 + SCC（全球 gsc / 分国 rsc）
```

---

## 3. 理论方法

本节自洽地给出从计量识别到 SCC 的全部理论。记号约定：c 为国家，t 为年份，h/s 为 horizon（0..10），G/L 为全球/局部渠道。

### 3.0 变量字典（速查表）

公式与代码中反复出现的核心变量，**尤其注意区分"绝对量"与"中心化量"**——历史上 `ypc0`/`temp1990` 误写成中心化值造成过全链路失真，当前 `open_issues.md` P0-A 的绝对温度/距平混装也属同类：

| 变量 | 含义 | 口径 | 主要来源/文件 |
|------|------|------|--------------|
| `barT_s` | 国家 1950–1979 气候态 | **绝对温度**（°C） | `analysis.dta` |
| `barT_c` | 中心化气候态 | `barT_s − barT_W` | `struct_params.csv` / LP 端 |
| `barT_center` | LP 温度中心化常数（= barT_W） | 跨国**无加权**均值 ≈ 19.53 | `damage_state_centers.csv` |
| `lny0_s` | 初始 log 收入 | **绝对** log 人均 GDP | `analysis.dta` |
| `lny0_c` | 中心化初始收入 | `lny0_s − lny0_W` | `struct_params.csv` / LP 端 |
| `lny0_center` | LP 收入中心化常数（= lny0_W） | 跨国无加权均值 ≈ 8.33 | `damage_state_centers.csv` |
| `regtmp` | 区域/国家绝对温度 | MESMER 分解式输出 | `climateregional` |
| `region_lv`（ν） | 局部天气变率（完整） | 近似零均值平稳噪声，**不减参考年**；仅供 `climateregional` 生成国别温度 | MESMER emulator |
| `region_lv_local` | ν 中纯特异性的局地分量 | 损害组件局地渠道 ζ_l 的驱动项（P1-19） | MESMER emulator |
| `global_lv`（gv） | 全球年际变率 | 并入损害组件全球渠道：ζ_g·(T_global − T0 + gv)（P1-19） | MESMER emulator |
| `T_global` | 全球温度路径 | 默认由 FaIR 内生生成；`climate=:csv` 时读取外生 CSV | `climate_T_global.csv` / FaIR |
| `T0_global` | 全球损害参考温度 | `bau_ref_year` 年的 T_global，Base/Pulse 共用 | `MimiCLAIM.jl` |
| `T_state` | 动态温度状态 | 过去 30 年 `regtmp` 滚动均值减 `barT_W`（**含天气噪声**，与 `barT_s` 同源） | `components_damage/claim_damage.jl` §3.4.1 |
| `Y_state` | 动态收入状态 | 滞后一期含损害收入的 log 减 `lny0_W` | `components_damage/claim_damage.jl` §3.4.1 |
| `theta_g0/g1/g2` | 全球渠道 LP 系数（θ₀,θ₁,θ₂ per horizon） | pooled draw；国家异质性由运行时状态投影产生 | `damage_theta_base.csv` / `damage_theta_draws.csv` |
| `theta_l0/l1/l2` | 局部渠道 LP 系数 | 同上 | `damage_theta_base.csv` / `damage_theta_draws.csv` |
| `Z1/Z2` | TFP 损害核基底 | 纯 ζ 核（形状冻结） | `damage_basis.csv` / `damage_structural_params.csv` / `damage_weights.csv` |
| `Y1/Y2` | 产出响应基底 | ζ 经温度 IRF 卷积 + Solow 传播后 | `damage_basis.csv` / `damage_structural_params.csv` / `damage_weights.csv` |
| `z_tfp`（z） | TFP 对数偏离 | 全球 + 局部两渠道卷积之和 | `components_damage/claim_damage.jl` |
| `khat`（k̂） | 资本对数偏离 | Solow 递归状态变量，无限记忆 | `components_damage/claim_damage.jl` |
| `yhat`（ŷ） | GDP 对数偏离 | `z_tfp + alpha_c × khat` | `components_damage/claim_damage.jl` |
| `ypc_baseline` | 无损害基准人均收入 | 外生社会经济路径（**绝对**美元） | `socioeconomic` |
| `ypc_new` | 含损害人均收入 | `ypc_baseline × exp(yhat)` | `components_damage/claim_damage.jl` |
| `loss` | 损失金额 | `(ypc_baseline − ypc_new) × population` | `components_damage/claim_damage.jl` |

### 3.1 温度冲击识别（Hamilton 滤波）

B&K 用 Hamilton (2018) 回归残差法构造温度冲击。对全球温度（Berkeley Earth 陆海合并距平 `berk_anom`，1850–2024 全样本）：

$$T^G_t = \beta_0 + \beta_1 T^G_{t-2} + \beta_2 T^G_{t-3} + \varepsilon_t, \qquad \hat\varepsilon_t \equiv g\_shock_t$$

即"当期温度中无法被 2–3 年前温度预测的部分"。**海洋分量至关重要**：只用陆地温度构造的全球冲击无法识别全球损害渠道（这正是 B&K 的中心论点——海洋含全球温度冲击代表气候系统的持续状态偏移，其损害远大于局部天气波动）。

局部冲击对每国温度 `lctmp_bkly_pw` 做同样的 Hamilton 滤波（容许国家专属截距与斜率），再**正交化**于全球冲击：

$$l\_shock_{c,t} = \gamma_c + \lambda_c \cdot g\_shock_t + u_{c,t}, \qquad \hat u_{c,t} \equiv l\_shock\_orth_{c,t}$$

正交化保证两渠道可加、系数为偏效应：`g_shock` 捕捉全球气候状态偏移的**总**效应（含其在各国的温度放大），`l_shock_orth` 捕捉纯特异性局部天气。

### 3.2 双渠道状态依赖局部投影

对每个 horizon h = 0..10 估计一个面板 LP：

$$\Delta_h y_{c,t} \equiv \ln y_{c,t+h} - \ln y_{c,t-1} = \theta^G_h(c) \, g\_shock_t + \theta^L_h(c) \, l\_shock\_orth_{c,t} + \Gamma' X_{c,t} + \alpha_c + \epsilon_{c,t+h}$$

其中响应系数是**状态依赖**的：

$$\theta^j_h(c) = \theta^j_{0,h} + \theta^j_{1,h} \cdot barT_c + \theta^j_{2,h} \cdot lny0_c, \qquad j \in \{G, L\}$$

- $barT_c = barT_s - barT_W$：国家 1950–1979 气候态减去跨国无加权均值（**中心化**）；
- $lny0_c = lny0_s - lny0_W$：初始 log 收入的中心化值；
- 二者时不变，被国家固定效应吸收，故回归中只出现 shock 交互项。

**关键识别设定**：

1. **不放年份固定效应**——年份 FE 与全体国家共享的 g_shock 完全共线，会吸收掉全球渠道；
2. 控制变量对齐 B&K：`recession L(0/2)`、`dlnpoil / treasury1y / dlnrgdpcw L(1/2)`（仅滞后项，避免坏控制）、降水 `L(1/2)`、地区专属线性趋势；
3. 标准误：Driscoll–Kraay（`xtscc`，lag = max(4, h+1)），对全球共同冲击的截面相关稳健；
4. **无加权（unweighted）回归，"一国一票"**——历史识别阶段坚决不用 GDP 加权。若按 GDP 加权，全球约 80% 的权重将集中在欧美等温带富裕国，会掩盖低收入热带国家的气候脆弱性，并压缩状态变量 $(barT_c, lny0_c)$ 的截面识别变异——而状态交互项正是靠"穷热国 vs 富冷国"的对比识别的。GDP 的规模影响只在最后一步还原：SCC 加总时按 $loss = (y^{base} - y^{base}e^{\hat y})\times pop$ 以绝对体量计价（见 §3.4.3、§3.6）。"识别用等权、计价用体量"的分工是本管线的刻意设计（`02_lp_regression.do` §10 标注 unweighted）。

每个 h 的 6 参数向量 $\vec\theta_h = [\theta^G_{0,h}, \theta^G_{1,h}, \theta^G_{2,h}, \theta^L_{0,h}, \theta^L_{1,h}, \theta^L_{2,h}]^\top$ 与完整 6×6 协方差矩阵 $\Sigma_h$（含 9 个 G–L 跨渠道协方差）一并导出至 `lp_main.csv`。

**国家 θ_c 的插补**：直接逐国回归会过拟合（全球冲击的有效时间点仅 ~19–30 个），因此每国的 $\theta_c(h)$ 由 pooled 系数与该国 $(barT_c, lny0_c)$ 线性投影得到（`lp_imputed.csv`），其标准误由完整 VC 矩阵的 delta 方法给出：

$$SE^2_{\theta_c} = V_{00} + barT_c^2 V_{11} + lny0_c^2 V_{22} + 2\,barT_c C_{01} + 2\,lny0_c C_{02} + 2\,barT_c\,lny0_c C_{12}$$

### 3.3 从 IRF 到结构损害核：Solow 反演

#### 3.3.1 为什么要反演

LP 给出的 $\theta_h$ 是**总量 GDP 响应**，混合了三层机制：(i) 温度冲击自身的持续动态；(ii) TFP 受损；(iii) 投资下降导致的资本存量伤痕。要在 IAM 中前向模拟任意温度路径的损害，必须剥离出**结构性 TFP 损害核** $\zeta_s$：温度水平偏离 1°C 在 s 年后对 TFP 的对数损害。

**为何不直接用 B&K 的参数**：B&K 在 Ramsey（最优储蓄）动态下回收核参数；本项目的 IAM 采用 Solow（固定储蓄率）资本动态。同一 $\theta_h$ 在不同资本动态下隐含不同的 ζ，直接移植会造成核错配。因此用本项目的"非平稳资本动态框架"重新反演。

#### 3.3.2 模型结构（正向映射）

三层线性算子的复合：

**第一层——温度动态**。1 单位 Hamilton 冲击触发温度水平路径 $\phi^T_h$（由辅助 LP 估计，存 `temp_irf.csv`）：

$$\phi^{T,j}_h = \frac{\partial T^j_{t+h}}{\partial shock^j_t}, \qquad j \in \{G, L\}$$

**第二层——TFP 卷积**。温度水平偏离通过核 ζ 转为 TFP 对数偏离：

$$z_h = \sum_{s=0}^{h} \zeta_s \, \phi^T_{h-s}$$

**第三层——Solow 资本传播**。TFP 偏离进入线性化 Solow 增长模型：

$$\hat y_h = z_h + \alpha_c \hat k_h, \qquad \hat k_{h+1} = \rho_c \hat k_h + (1-\rho_c)\hat y_h, \qquad \hat k_0 = 0$$

其中 $\rho_c = \dfrac{1-\delta_c}{1+g_c}$；$\alpha_c$（资本份额 = 1−labsh）、$\delta_c$（折旧率）、$g_c$（人均增长率）逐国取自 PWT 11.0（`struct_params.csv`）。

于是理论预测 $\theta_h = \text{Solow}\big[(\zeta * \phi^T)_h\big]$，反演即在此约束下由经验 $\hat\theta_{c,h}$ 求 $\zeta$。

#### 3.3.3 损害核参数化（2026-08-04 起为单核；此前为双核竞争）

$$\boxed{\;\zeta_s = A e^{-Bs} + C\,s^{\nu} e^{-Ds}\;}\qquad \nu = 3\ \text{（校准常数，非估计量）}$$

即期打击 $Ae^{-Bs}$ 以速率 $B$ 衰减；延迟驼峰 $Cs^{\nu}e^{-Ds}$ 峰点在 $s^*=\nu/D$、峰高 $|C|(\nu/D)^{\nu}e^{-\nu}$，正比于 Erlang($\nu{+}1$, $D$) 密度。外推单调归零，不震荡。

> **原 M2（阻尼振荡核 $e^{-Bs}[A\cos\omega s + C\sin\omega s]$）已删除**（`open_issues.md` P1-12 第二次追加）。三条理由：
> 1. 其外推段 $\zeta$ 符号变化中位 **23 次**（最大 90），窗口外积分与窗口内**反号**（全球 −0.339），抵消掉样本内 34% 的损害。一次"先跌后回"可解释，其后反复升降没有机制对应。
> 2. M1/M2 的分族 **99.5%** 可由 $(\bar T_c,\ln y_{0c})$ 一条直线复现——$\theta$ 矩阵秩恰为 3，分族是几何位置而非经济机制。本指南此前把 M1 获胜国解释为某类经济体的说法一并作废。
> 3. 绝对判据下持平：$z_{rms}$ 中位 M2-only 0.812 vs 本核 0.810（理论值 0.798）。
>
> **$\nu$ 由 1 提到 3**：$\nu=1$ 时 $D$ 同时决定峰点与峰宽，坑要深就必须窄，拟合不出 $h=5\!-\!7$ 的宽坑（等权拟合下失拟峰值在 $h=6$，$h\ge8$ 只占缺口 25%，与权重无关）。$\nu$ 是**校准值**：交给数据它会跑到网格上限；取 3 的依据是 $z_{rms}$ 恰好停在 $\sqrt{(11-4)/11}=0.798$ 上（$\nu=1$→0.988，$\nu=4$→0.758 已在拟合噪声）。

**为何放弃 B&K 的双指数核**：B&K 原文采用 $\zeta_s = A(e^{-Bs} - e^{-Cs})$。该函数在 $s=0$ 处**数学上强制 $\zeta_0 = 0$**——冲击当期生产率零损害。由于资本存量当期前定（$\hat k_0 = 0$），这等价于把第一年的 GDP 损失**先验地**锁死为零。本项目的核用独立的即期分量显式捕捉当期打击（$\zeta_0 = A$），把"当期是否有损害"交给数据决定而非由函数形式预设，且满足 $\lim_{s\to\infty}\zeta_s = 0$（transient 损害，无永久增长率折损）。详见 `02regression/damage_inversion_plan.md` §2.1。

> **一处此前的表述错误（2026-08-01 更正）**：本节旧版写"直接违背实证数据（LP 在 h=0 已有**显著**负效应 $\theta_0 \neq 0$）"。实测**不成立**：全球渠道 $\lvert\theta_c(0)\rvert/\mathrm{SE}_c(0)$ 中位 0.96，**0/185 国达到 $\lvert t\rvert>1.96$**，pooled $b_g(0)$ 的 $t = -1.05$。局地渠道稍好（中位 1.55，63/185 国显著）。
>
> 放弃 B&K 核的正当理由是**不预设** $\zeta_0=0$（让数据说话），不是"数据已证明 $\zeta_0\ne0$"。这个区别很重要：$A$ 解锁后（见 3.3.3.2），全球渠道有 92 国的 $A$ 由负转正、$\lvert A\rvert$ 中位减半——数据在全球渠道上并不支持一个明确的即期打击。这一点应在论文中如实交代。

估计规则（`03_invert_damage.py`）：

- **$(A, B, C, D)$ 联合加权 NLS**，残差覆盖 $h=0..10$，权重为国家级逆方差 $w_h = 1/SE^2_{\theta_c(h)}$。$A$ 不钳制（极寒国家的正红利照收）。$A$ 曾由 $h=0$ 解析锁定，2026-08-01 改为联合估计，见 3.3.3.2；
- 边界按**意图**表达：$B \in [0.05, 5]$；$D \ge \nu/H$（驼峰峰点须落在 11 年估计窗口内）；$|C| \le 0.37\,e^{\nu}(D/\nu)^{\nu}$（驼峰峰高 $\le0.37$，实现为归一化 $c_{norm}\in[-1,1]$，因该界不是箱约束）；
- **不再有核选择环节**（原按 $R^2_w$ 在 M1/M2 间择优，见 3.3.3.1 的历史记录）。

**边界必须从意图推导，不可写死数字。** 驼峰峰值 $|C|(\nu/D)^{\nu}e^{-\nu} \le 0.37$ 是意图——把单次 1°C 冲击的极限延迟生产率偏离锁在 37% 以内，杜绝小国噪声把核拟合成"生产率归零"的过拟合灾难。但它到 $C$ 的换算随 $(\nu, D)$ 变化。本项目为此犯过**三次同型错误**（`open_issues.md` P1-12）：$|C|\le0.05$ 只在 $\nu=1,D=0.05$ 处等价；$D\ge0.05$ 的"20 年内见顶"只在 $\nu=1$ 成立（$\nu=3$ 时允许峰点到 60 年，RUS 因此外推到 $+115\%$）；而 20 年本身就超出 11 年窗口。$B$ 的两端另有专门论证，见 3.3.3.1。

#### 3.3.3.1 $B$ 的边界与择优判据（2026-07-31 定案；**择优部分已于 2026-08-04 作废**）

> **本节的"择优判据"部分是历史记录。** M2 已删除，不再有 M1/M2 择优；下文所有 M1=xx / M2=xx 的计数、$R^2_5$ vs $R^2_w$ 的判据比较，均描述 2026-08-04 之前的状态，保留是为了说明当时为何那样选。
>
> **$B$ 的边界论证仍然有效，但 $B$ 的角色变了**：$\nu=1$ 且 $B$ 撞上界时尾巴由驼峰项承担，$\nu=3$ 后 $B$ 跑到**下界**（全球渠道 57/185 国，GDP 45%，其中 57 国 $A>0$），尾巴变成 $A$ 项——即一个半衰期 13.9 年的持续正效应。撞界没有消失，只是换了位置；原 M2 靠振荡使窗口外积分相消把它藏起来，本核不相消于是暴露。**$B$ 仍须改为声明的校准值 + 敏感性带，未做**（`open_issues.md` P1-12）。

$B$ 的两个端点上有大量国家聚集（全球渠道 181/185），曾被登记为"参数退化到边界"（`open_issues.md` P1-12）。逐国实测后结论是**边界数值无需修改**，但两端的性质完全不同，必须分开理解。

**上界 $B \le 5$ 是可识别区的边缘，不是"太紧"。**
固定 $(C, D)$ 在最优、只扫 $B$ 的加权 SSR 剖面显示：

| $B$ | 1.0 | 2.0 | 5.0 | 10 | 20 | 35 | 50 |
|---|---|---|---|---|---|---|---|
| $\mathrm{SSR}/\min$（全球渠道中位） | 1.1530 | 1.0388 | 1.0015 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |

$B \ge 10$ 后 SSR 平到小数点后四位，$B = 5$ 相对最优只差 **0.145%**。原因是 $e^{-5\cdot 1} = 0.0067$，$Ae^{-Bs}$ 在 $s \ge 1$ 已经消失，而 $A$ 本就由 $\theta(0)$ 解析锁定——**数据对大 $B$ 没有任何信息**。实测把上界放宽到 50，全球渠道的 $B$ 散布在 26–46（p10 26.4 / 中位 34.5 / p90 46.0），那是优化器在完全平坦的似然面上的随机落点，不是估计量。

因此上界在这里起的是**把不可识别方向正则化到一个点**的作用。"撞上界"的正确读法是"该国的即期冲击是一次性的，持续性全部由驼峰项 $Cse^{-Ds}$ 承担"，而不是数值病态。

**但可识别性是逐国的，绝不能把 $B$ 从 M1 整体删掉。** 局地渠道选中 M1 的 41 国里有 40 国 $B < 10$ 且**强可识别**：把它们的 $B$ 压到 5 会让 SSR 中位上涨 **596%**（JPN +9745%、LSO +9141%、PRT +7434%）。$B$ 大时衰减在离散年度网格上一步完成、看不见；$B$ 小时衰减过程在 11 年窗口内清晰可见——这是可分离指数模型的经典性质。

**下界 $B \ge 0.05$ 是外推先验，不是拟合约束。**
这一端 $B$ 强可识别（$B = 0.05$ 比 $B = 0.02$ 的 SSR 高 9.8%，$B = 0.20$ 高 66%），数据确实要 $B \to 0$，即 11 年窗口内看不到阻尼。但那只是关于 11 年的陈述。$B = 0.05$ 时 100 年后包络残余 $e^{-5} = 0.7\%$，$B = 0.02$ 时为 13.5%，核就不再代表"暂时性水平损害"这一定性主张。实测下探到 0.02 的收益是 $R^2_5$ 中位 +0.005，代价是 100 年残余 >1% 的国家渠道数由 33 升至 191，且撞界只是原地挪了一格（95/95 撞 0.05 → 101/101 撞 0.02），故回退。

> 注意这里的外推问题是**解释性的，不是数值稳定性的**：`claim_damage.jl` 的卷积窗口硬截断在 `for s in 0:min(10, ...)`，$s > 10$ 从不进入 SCC 计算。

**择优判据由 $R^2_5$ 改为 $R^2_w$。**
此前用可信段 $R^2_5$（$h \le 5$，不加权）择优，理由是避免 M2 追逐 $h = 8..10$ 的噪声峰。三种判据的实测对比（185 国 × 2 渠道）：

| 判据 | 定义 | 全球 G | 局地 L |
|---|---|---|---|
| $R^2_5$ | $h \le 5$，不加权 | M1=90 M2=95 | M1=41 M2=144 |
| $R^2_{full}$ | $h \le 10$，不加权 | M1=18 M2=167 | M1=42 M2=143 |
| $R^2_w$ | $h \le 10$，权重 $1/SE_h^2$ | **M1=64 M2=121** | **M1=39 M2=146** |

$R^2_{full}$ 让 SE 大 2.6 倍的长端与短端等权，全球渠道 M1 塌到 18，正是当初担心的噪声追逐。$R^2_w$ 同时解决两头：用满 11 个点，又按精度自动折价长端；更关键的是它**与 NLS 的目标函数（加权 SSR）口径统一**——此前"拟合用加权 $h=1..10$、择优用不加权 $h=0..5$"是内部不一致。

改判据后 46/370 个国家渠道换核，全球渠道 $R^2$（全段）中位由 0.5672 升至 0.5908，$R^2_5$ 中位由 0.6760 降至 0.6695（不再被直接优化，属预期）。$R^2_5$ 与 AIC 继续写入 `damage_params.csv` 作诊断，只是不再参与择优。

上表三行均在 $A$ 锁定下测得。**2026-08-01 解锁 $A$ 后（见 3.3.3.2），$R^2_w$ 下的择优结果变为全球 M1=19 M2=166、局地 M1=40 M2=145。判据本身未变。**

全球渠道的迁移是**单向的：45 国 M1→M2，0 国 M2→M1**；局地渠道几乎不动（仅 1 处翻转）。机制是 $A$ 不再被钉在 $\theta(0)$ 上后，M2 的 $\zeta_0 = A$ 与振幅包络 $\sqrt{A^2+C^2}$ 解耦，可以同时容纳"$h=0$ 近零、中段下探、长端回摆"——而这正是全球渠道 IRF 的典型走势。实证支持：M2 胜出国的 $\lvert A\rvert/\sqrt{A^2+C^2}$ 中位在**全球渠道仅 0.121**（$A$ 近零、振幅全由 $C$ 承担），在局地渠道则是 **0.863**（$A$ 主导，无解耦空间）——两个渠道差一个数量级，正对应两者迁移量的差别。

**一个应当在论文中如实交代的性质**：分歧国家的判据边际很薄（$|\Delta R^2_5|$ 中位 0.0365，$|\Delta R^2_w|$ 中位 0.0219），说明 M1/M2 的二选一本身是脆弱的，不宜当作强结论。

#### 3.3.3.2 为何不再把 $A$ 锁定在 $h=0$（2026-08-01）

恒等式 $\hat y_0 = A\,\phi^T_0$ 本身正确（$\hat k_0=0$），但"等式正确"不等于"应把 $A$ 完全交给单点 $\theta(0)$ 决定"——$\theta(0)$ 恰是整条 IRF 上最不可靠的点（全球渠道 0/185 国显著）。锁死后，核的其余自由度被迫去补偿这个噪声点。

为选出正确的修法，测了四个变体，**全部对原始 $\theta$ 评分**（避免"核拟合自己平滑过的目标"的循环）：

| 变体 | $R^2_w$ 中位 | $z_{\mathrm{rms}}$ 中位 | $\lvert A\rvert$ 中位 |
|---|---|---|---|
| A. $A$ 锁定于 $\theta(0)$（旧） | 0.6829 | 0.966 | 0.0162 |
| B. $A$ 锚于平滑后的 $\theta(0)$ | 0.7569 | 0.806 | 0.0079 |
| C. 完整 Smooth LP | 0.7516 | 0.813 | 0.0079 |
| **D. $A$ 自由（现行）** | **0.7594** | **0.782** | 0.0078 |

收益分解：$A\to B$（只换锚点）**+0.0740**，$B\to C$（再跨 horizon 平滑）**−0.0053，0/185 国改善**。

**Smooth LP（Barnichon & Brownlees 2019）因此未被采纳**：全部收益来自 $A$ 的锚点，跨 horizon 平滑零收益——下游三参数核本就是平滑器，SLP 属重复劳动。而 SLP 会引入光滑度调参 $\lambda$、收缩阶 $r$ 与样条设定三处选择，其置信带自陈为启发式且无渐近理论，在短 horizon 小样本处覆盖率最差（原文报 $T=50$、$h=2$ 时 0.588，名义 90%）——恰是 $A$ 所在的位置。

变体 D 不引入任何新调参、不改估计链、保持闭式，**185/185 国改善**。代价是 $h=0$ 残差不再恒为 0（中位约 1.0 个 SE），这是有意的取舍。

**为何 SCC 几乎不动（484.57 → 484.04，−0.1%）**：SCC 由持续升温下的**累积核** $\sum_s \zeta_s$ 驱动，而非单点 $\zeta_0 = A$。$A$ 减半被 $C$ 分量补偿，累积核基本守恒——跨国相关 **0.995**，中位变化 −4.3%（全球）/ −14.4%（局地）。

#### 3.3.4 形状冻结与闭式反演（运行时提速的关键）

注意核对 $(A, C)$ **线性**。冻结形状参数 $(B, D)$ 后，预计算两套基底：

- **Z 基底**（纯 TFP 核）：$Z_1[s] = e^{-Bs}$、$Z_2[s] = s^{
u} e^{-Ds}$；
- **Y 基底**（经 φ 卷积 + Solow 传播）：$Y_k[h] = \text{Solow}\big[(Z_k * \phi^T)\big][h]$。

由线性性 $\hat y_h = A\,Y_1[h] + C\,Y_2[h]$，反演退化为一个 $2\times2$ 的**加权闭式 OLS**：

$$\begin{pmatrix} \sum_h w_h Y_1^2 & \sum_h w_h Y_1 Y_2 \\ \sum_h w_h Y_1 Y_2 & \sum_h w_h Y_2^2 \end{pmatrix}\begin{pmatrix} A \\ C \end{pmatrix} = \begin{pmatrix} \sum_h w_h Y_1 \theta_h \\ \sum_h w_h Y_2 \theta_h \end{pmatrix}, \qquad h = 0..10$$

每次反演约 11 次点积加一个 $2\times2$ 解（微秒级），使"每国×每年×每 MC 抽样"的动态重反演在计算上可行。基底缓存于 `data/unshared_parameters/damage_basis.csv`，与上游 NLS 的口径一致（两处都让 $A$ 参与拟合）。

#### 3.3.4.1 权重口径 $w_h$：为何保留 $1/SE_h^2$（2026-08-01）

$w_h$ **在两处起作用，量化时必须分开**：

1. **NLS 拟合**（`03_invert_damage.py::_fit`）——决定形状 $B, D$；
2. **运行时闭式 OLS**（`claim_damage.jl:124-132`）——每模型年 × 每国 × 每 MC 抽样各解一次上式，决定每个未来状态下的高度与驼峰配比，直接进 SCC 积分。

两处用同一向量、作用于不同回归，**替换权重时符号相反**：

| 口径 | NLS 权重 | 运行时权重 | SCC(2030, ssp245, →2100) | 相对现状 |
|---|---|---|---|---|
| **现状** | $1/SE^2$ | $1/SE^2$ | **484.04** | — |
| DL 全链条 $1/(SE^2+\hat\tau^2)$ | DL | DL | 485.71 | +0.34% |
| 等权 全链条 | $1$ | $1$ | 508.67 | **+5.09%** |
| ├ 仅运行时等权 | $1/SE^2$ | $1$ | 568.43 | +17.43% |
| └ 仅 NLS 等权 | $1$ | $1/SE^2$ | 419.32 | −13.37% |

**只覆盖运行时权重的实验是内部不自洽的**（形状仍停留在逆方差拟合的结果上），会把杠杆高估约 3 倍。端到端敏感度是 **+5.09%**。

三条候选改法全部否决：

- **状态依赖权重**（把 $w_h$ 也投影到 $(T_{state}, Y_{state})$）：整体缩放在闭式解中精确抵消，归一化权重廓形相关 $\ge 0.998$，影响在千分位。只记文档。
- **随机效应权重 DL**：$A$ 解锁后全球渠道 $I^2$ 由 0.25 降到 **0.034**、$Q/df$ 由 1.33 降到 1.035——原先的"设定误差"主要是 $A$ 被锁死在噪声点上制造的。DL 只动 +0.34%，0/370 处择优翻转，不引入。
- **GLS（用完整跨 horizon 协方差 $\Sigma_c^{-1}$）**：$\Sigma_c$ 可由 `dk_scores.csv` 的 FWL 得分直接重建（bread 退化为标量，无需新 Stata 运行；反算的 $SE_\theta(h{=}5)$ 对实际值比值中位 0.993）。重建证实跨 horizon 相关极高（平均相邻 0.827），**但也证明 $\Sigma^{-1}$ 不可估**：33/185 国 $\Sigma_c$ 非半正定，有效独立维数中位仅 **2.23**，而估计它需要用 68 年数据定 66 个参数。任何 GLS 都要靠 ridge/特征值截断，正则化参数将成为决定答案的东西。否决。

**结论**：$\Sigma^{-1}$ 不可估，**在可估的权重族里 $1/SE_h^2$ 是自然且可辩护的选择**。等权并非"中性"——它断言 $SE_h$ 跨 $h$ 不变，而实测 min/max 权重比 0.147。剩余口径任意性由**等权全链条 +5.09%** 界定，作为稳健性上界报告。详见 `open_issues.md` P1-6。

### 3.4 前向损害模块（claim_damage 的方程组）

#### 3.4.1 动态状态更新

拒绝"静态基准谬误"：极寒国家不应在升温 3°C 后仍机械享受正红利。在未来年份 t，用**滞后一期**状态重新投影 θ：

$$\hat\theta^j_{c,h}(t) = \theta^j_{0,h} + \theta^j_{1,h} \, T^{state}_{c,t-1} + \theta^j_{2,h} \, Y^{state}_{c,t-1}$$

状态变量与 LP 端严格同口径（中心化常数 $barT_W, lny0_W$ 存 `damage_state_centers.csv`）：

$$T^{state}_{c,t} = \underbrace{\frac{1}{30}\sum_{i=t-30}^{t-1} regtmp_{c,i}}_{\text{过去 30 年气候态（含天气，与 } barT_s \text{ 同源）}} - \; barT_W, \qquad Y^{state}_{c,t-1} = \ln y_{c,t-1} - lny0_W$$

**温度状态用 30 年滚动均值**（2026-07-30 修，`open_issues.md` P1-4）。此前用的是滞后一年的 $regtmp - region\_lv$，但 LP 的 $barT_s$ 是 **1950–1979 三十年**气候态均值，"一年值"与"三十年均值"语义不同——接入 FaIR + MESMER 后 $ts\_lt$ 自身的年际波动会直接进入 θ 投影。三个设计决定：

- **窗口长度 30 不是自由参数**，由 $barT_s$ 的窗口钉死（`helper_damage.jl::BKW_STATE_WINDOW`）。窗口为 $[t-30,\,t-1]$，不含当年；历史不足 30 年时用扩张窗口，绝不静默截断。
- **不扣减 $region\_lv$**。$barT_s$ 本身就是实测 `lctmp_bkly_pw` 的三十年均值（含天气），WMO 的"气候常态"定义亦然，含 ν 才是同源口径。30 年平均后残留噪声中位仅 0.077 °C，约为 $barT$ 跨国 sd（7.49 °C）的 1%。**代价**：θ 因此依赖 MESMER 随机种子，会轻微污染 P1-7 外层气候 MC 的分解；SCC 不受影响（Base/Pulse 共用同一条 ν，均值之差只剩受迫分量）。
- **θ 随气候态时变，且不引入衰减参数**。$\theta_1$ 识别的是**截面**梯度，用于时间演化属可迁移性假设（§3 caveat 1）。但"保留截面异质性、冻结其时间演化"内部不自洽：若认为 $\theta_1$ 纯属混淆，诚实的应对是根本不用它（退回 `damage_function=:BK_global`），而非横截面上用、时间上不用。若 $\theta_1$ 确有混淆成分，真实的国内梯度满足 $0 \le |\theta_1^{within}| \le |\hat\theta_1|$ 且同号，故当前设定落在**激进那一端**。

**可迁移性假设在 t=0 就已生效**：实测损害起始年（2025）各国 $T^{state}$ 均值 +0.700，而 LP 的 $barT_c$ 均值 −0.102，漂移 **+0.802 °C**（range [+0.488, +1.204]）——这是 1950–79 到 1995–2024 已实现的升温。折算到 θ 上平移 −0.0071，占 $\theta^G_0$ 峰值（−0.141）的 **5.1%**（最差 7.6%）。不是到 2100 才发生的事。

其中 $y_{c,t-1}$ 是**上一期含损害**的实际人均 GDP，对应代码中的 `ypc_new[t-1, c]`（包括下限/上限钳制后的实际输出值）；在 `damage_start_year` 当年也使用模型上一年的收入状态，只有当损害从模型第一年开始、没有可用上一期时才退回当期基准收入。这一设定使收入状态反馈内生演化：多数仍保持增长的国家会随 `ypc_new` 上升而逐步进入更富裕状态；受损严重的低收入国家则可能因 `ypc_new` 被压低而面对更大的边际损害，形成可用于研究的 climate-induced poverty trap / 气候贫困陷阱机制。随后将 $\hat\theta(t)$ 送入 §3.3.4 的闭式 OLS 得到当年的 $A_c(t), C_c(t)$，组装当年的损害核 $\zeta_s(t) = A(t) Z_1[s] + C(t) Z_2[s]$。

#### 3.4.2 状态空间递归（资本记忆截断修复）

TFP 卷积截断在 s = 10 是合理的（ζ 指数衰减归零），但**资本存量有无限记忆**——10 年前的冲击通过投资渠道在 $\hat k$ 中留下伤痕，以速率 $\rho^t$ 缓慢愈合。若直接用 Solow 传播后的 $\hat y$ 核做 11 年截断卷积，会让第 11 年 GDP 瞬间跳回基准线。正确做法是把 TFP 卷积与资本递归**分离**：

$$\begin{aligned} z_{c,t} &= \sum_{s=0}^{10} \zeta^G_s(t) \cdot \big(T^{global}_{t-s} - T^{global}_{ref}\big) + \sum_{s=0}^{10} \zeta^L_s(t) \cdot \nu_{c,t-s} && \text{（TFP：有限记忆卷积）} \\ \hat y_{c,t} &= z_{c,t} + \alpha_c \hat k_{c,t} && \text{（GDP 对数偏离）} \\ \hat k_{c,t+1} &= \rho_c \hat k_{c,t} + (1-\rho_c) \hat y_{c,t} && \text{（资本：无限记忆递归）} \end{aligned}$$

两渠道的温度输入口径**不同**，这是实现的关键细节：

| 渠道 | 输入 | 是否减参考年 | 理由 |
|------|------|-------------|------|
| 全球 G | $T^{global}_t - T^{global}_{ref}$ | **减**（ref = `bau_ref_year`，默认 2025，Base/Pulse 两轨共用） | 非平稳气候趋势，损害由状态偏移驱动 |
| 局部 L | $\nu_{c,t} = region\_lv_{c,t}$ | **不减** | MESMER 天气噪声近似平稳、均值为零；若减 t₀ 会把参考年的一次随机热浪固化为永久强迫 |

**"均值为零"现已严格成立（`open_issues.md` P1-13，2026-08-04 修复）**：MESMER 拟合出的 AR(2) 全球变率截距隐含 **+0.0192 °C** 的平稳均值，追查后确认它是**趋势估计的残渣、不是物理变率**——主因是 `train_gt_ic_OLSVOLC` 的火山回归**不带截距**（AOD 非负、系数为负，把 `gt` 整体压低 0.0263），次因是 LOWESS 在趋势曲率处的偏差。由于 `claim_damage` 把 `global_lv` 直接加到 `T_global` 上（P1-19），该残渣等价于一个永久的虚假强迫。

现由 `export_grid_params.py` 在导出 `climate_grid_gv.csv` 时置 `ar_int = 0`（`--keep_gv_mean` 保留原值仅供审计）。重建验证：除 gv 与 manifest 外全部输出逐字节一致；`loc` 逐值不变、`gv` 精确下移 μ、`total` 下移 $W(c^{gv}\mu)$，残差 9.5e-16。净影响 BAU 累计损害 **−1.368%**、SCC −0.020%、2100 国别损害中位 −1.29%。

> **仍需注意单条实现的抽样运气**：跨 41 个 seed 的 gv 均值中位已落在 0 附近，但生产用的 seed=42 在 2025–2150 窗口内是一次偏暖抽样（gv 均值 +0.0271，第 98 百分位；BAU 损害同样第 98 百分位，两者 r = 0.964）。**SCC 不受影响**（第 51 百分位，Base/Pulse 共用同一条路径），但**报告 BAU 与国别损害水平时必须注明 seed 或跨 seed 平均**——seed=42 比跨 seed 均值高 +2.264%。

#### 3.4.3 映射回真实收入

$$y_{c,t} = y^{baseline}_{c,t} \cdot \exp(\hat y_{c,t}), \qquad loss_{c,t} = \big(y^{baseline}_{c,t} - y_{c,t}\big) \cdot pop_{c,t}$$

使用严格的 $\exp$ 而非 $(1+\hat y)$ 一阶近似——当 $|\hat y|$ 达 10% 时泰勒误差约 0.5 个百分点，在跨百年 SCC 积分中不可忽略。组件另设防爆钳制（ypc ∈ [200, 2×baseline]）并输出 `clipped` 计数供诊断。

### 3.5 不确定性框架

#### 3.5.1 因子化平滑 MVN（损害参数层）

> **口径提示（2026-08-06）**：本小节描述的共享 Z 抽样**已不是生产口径**。生产的 `damage_theta_draws.csv` 现在装的是 year-block bootstrap 的 replicate（300 条），跨 horizon 相关是实测的 **+0.42**，而不是下面推导出的 ±1。本小节保留，因为它解释了为什么这个近似曾是必要的、以及它到底假设了什么——**这些是理解 bootstrap 修正了什么的前提**。完整推导与新旧对比见 `damage_uncertainty_guide.md` §7–§9。

对第 m 次抽样，生成**单个** 6 维标准正态 $Z^{(m)} \sim \mathcal N(0, I_6)$，被全部 11 个 horizon 共享：

$$\vec\theta^{(m)}_h = \vec\theta^{base}_h + L_h Z^{(m)}, \qquad L_h = \text{Cholesky}(\Sigma_h)$$

- **horizon 内**：完整保留 6 参数联合协方差（含 G–L 跨渠道相关）；
- **horizon 间**：强制完美相关（秩-6 因子近似）——保证抽出的 IRF 曲线平滑不锯齿，代价是可能低估 horizon 正交方向的不确定性（文档如实标注为"因子化平滑"，year-block bootstrap 留作稳健性检验）。

**为什么不能改成跨 horizon 独立抽样**（这条对 bootstrap 口径同样成立，因为它约束的是"能不能逐点独立抽"，不是"相关取多少"）：若每个 h 各抽一个 $Z_h$，抽出的 IRF 将呈剧烈的高频上下锯齿——数值试验显示路径粗糙度比共享 Z 差约一个数量级。这种锯齿 θ 序列喂入 §3.3.4 的闭式 OLS 后，$C$ 的反演解会极度不稳定（残差 $R_h$ 被锯齿噪声主导，投影系数在边界间跳动），进而在动态状态更新中逐年放大，可能导致整条损害路径崩溃。共享 Z 的经济含义是保留"同一场估计误差冲击下，各 horizon 的系数同涨同跌"的平滑扰动结构——它模拟的是**整条 IRF 曲线**的抽样不确定性，而非逐点独立噪声。

注意抽的是 **pooled** 系数——国家异质性由运行时的 $T^{state}, Y^{state}$ 投影内生产生，因此 `damage_theta_base.csv` / `damage_theta_draws.csv` 不含国家维度。共享 Z 时是 1000 × 11 × 6；**生产的 bootstrap 表是 300 × 11 × 6**，抽样数上限由 replicate 数 $R$ 决定，`n` 超限会报错而非静默截断。

#### 3.5.2 嵌套 Monte Carlo（气候 × 损害）

```
for e = 1:M_c                 外层：气候不确定性（FaIR 参数集合的第 e 个成员）
    覆写 temperature:q / decay_factor、radiative_forcing:co2_f、co2_cycle:r0..rA
    for m = 1:M_d             内层：损害参数不确定性
        θ^(m) = θ_base + L_h·Z^(m)
        SCC[e,m] = compute_scco2(...)
```

**外层的内容已于 2026-08-02 修订**（`open_issues.md` P1-7）：

- **MESMER emulator 维度已删除**。实测换 seed 对 SCC 的影响仅 **0.19%**（6 个 seed，sd 0.92 \$/tCO₂）——Base/Pulse 共用同一条 ν，局地渠道在边际中相消，seed 只经 30 年滚动 `T_state` 这一条二阶通道进入。外层预算应全部给气候敏感度。
- **外层就是 FaIR 参数本身**，因此「缓存气候路径以跳过 FaIR」不可行——外层每个成员都必须重跑 FaIR。原文的该项优化只对「同一气候 × 不同 θ」的内层有意义。
- **成本估计更正**：实测稳态单次 SCC 约 3.3 s，其中 `get_model` 的参数装载占 2.84 s，FaIR 链本身仅 **0.17 s**。瓶颈在参数 I/O 而非气候模块，正确的优化是**建一次模型后原地 `update_param!`**，而非缓存温度路径。
- **逐成员覆写不需要改生产代码**：`temperature:q`、`temperature:decay_factor`、`radiative_forcing:co2_f`、`co2_cycle:r0/rU/rT/rA_co2` 本就是普通 `update_param!`（`helper_fairv2.jl:325-328, 409, 445-446`），在 `get_model` 之后、`compute_scco2` 之前覆写即可，Base/Pulse 两轨自动共用。

**气候维度的量级**（ssp245，pulse 2030，积分至 2100，基准 SCC = 484.04）：TCR 扫 AR6 likely range（1.4–2.2）使 SCC 落在 **[357, 637]**；用 Leach et al. (2021) 分发的 28 个 CMIP6 逐模式参数则是 **[417, 1018]**（中位 566.9）。这比损害参数权重口径的 +5.1% 上界大一个数量级。实测 SCC 与 **TCR** 的相关为 +0.975、与 ECS 仅 +0.763——积分到 2100 时脉冲响应是瞬态而非平衡响应；换更远的截断年该结论需重测。

### 3.6 SCC 的定义与计算

SCC 定义为：在年份 τ 额外排放 1 吨 CO₂ 所造成的、按社会折现率折现到 τ 年的全部未来边际损害。实现上用 **MarginalModel**（Base/Pulse 双轨差分）：

1. **Base 轨**：BAU 排放路径下运行完整模型；
2. **Pulse 轨**：在 τ 年注入排放脉冲（默认 `pulse_size = 1e7` 吨，分摊到 10 年，单位换算 CO₂→C 乘 12/44），其余完全相同；
3. **边际损害**：$MD_{c,t} = \dfrac{loss^{pulse}_{c,t} - loss^{base}_{c,t}}{pulse}$；
4. **折现**：按含气候损害的 Base 轨人均收入增长率进行 Ramsey 贴现；当前 `eta` 与 `prtp` 为所有国家共用的标量，逐国逐年递推折现因子

$$df_{c,t+1} = \frac{df_{c,t}}{1 + prtp + \eta \, g_{c,t}}, \qquad g_{c,t} = \frac{ypc_{c,t} - ypc_{c,t-1}}{ypc_{c,t-1}}$$

5. **加总**：

$$SCC_\tau = \sum_{t=\tau}^{T_{end}} \sum_c MD_{c,t} \cdot df_{c,t} \quad (\text{全球 gsc}); \qquad SCC^c_\tau = \sum_t MD_{c,t} \, df_{c,t} \quad (\text{分国 rsc})$$

**两个重要约定**：

- $T^{global}_{ref}$（损害参考温度）取 `bau_ref_year` 的 FaIR 温度，Base/Pulse **共用**——`compute_scco2(year=τ)` 只改变脉冲年份与折现起点，绝不重置参考温度，否则不同评估年的 SCC 差异会被参考点漂移污染；
- 局部渠道在边际中近似为零（两轨共用同一条局地变率路径（全球变率 `global_lv` 同理），$\delta\nu = 0$），但 BAU 全路径损失中局部渠道**不能丢**——它贡献真实的损害方差。

---

## 4. Mimi 组件模块详解

Mimi 中每个组件（`@defcomp`）声明 Parameter（外部输入/上游连接）、Variable（本组件计算），并实现 `run_timestep`。组件通过 `connect_param!` 串成有向计算图，按时间步推进。

### 4.1 FaIR v2.0 气候模块（`src/components_climate/`）

对 FaIR v2.0（Leach et al. 2021）简化气候模型的 Mimi 移植（源自 MimiFAIRv2.jl），把排放转为全球平均温度距平。

| 组件文件 | 功能 |
|----------|------|
| `emissions.jl` | 读取排放（`data/fair_data/rcmip_ssp*_emissions_1750_to_2500_rebased_obs.csv`，**2016–2024 已换成观测值**，见 §4.1.1），输出各气体年排放 |
| `co2_cycle.jl` | 4 池碳循环，含碳吸收随累积碳与升温的饱和反馈（state-dependent uptake），输出 CO₂ 浓度 |
| `ch4_cycle.jl` / `n2o_cycle.jl` | CH₄ / N₂O 单池寿命循环 |
| `montreal_gas_cycles.jl` | 蒙特利尔议定书气体（CFCs 等）批量循环，维度 `montreal_gases` |
| `flourinated_gas_cycles.jl` | 含氟气体（HFCs/PFCs/SF₆）批量循环，维度 `flourinated_gases` |
| `aerosol_plus_gas_cycles.jl` | 气溶胶及短寿命组分 |
| `radiative_forcing.jl` | 各组分浓度 → 辐射强迫（对数/平方根/线性形式），加总为总强迫 F |
| `temperature.jl` | 三热箱能量平衡：$T_{j,t} = T_{j,t-1} d_j + F_t q_j (1-d_j)$，$T_t = \sum_j \bar T_{j,t}$，输出全球温度距平 `T` |

参数（气体循环系数、热响应参数 q/d、初始条件）由 `helper/helper_fairv2.jl::update_MimiFAIR2_params!` 从 `data/fair_data/` 批量加载，支持 ssp119/126/245/370/585 等情景。

### 4.1.1 排放输入的观测校正：2016–2024 用实测，2025+ 重基准化

> **一句话**：模型读的**不是**原始 RCMIP 情景排放，而是 `*_rebased_obs.csv`——它的
> 2016–2024 段来自排放清单实测，2025–2500 段是原情景按 2024 年实测水平整体重基准
> 后的结果。原始 `*_rebased.csv` 保留在同目录作为出处，**模型不读它**。
> 生成脚本：`03Claim/scripts/build_observed_emissions.py`（2026-08-18 起）。

**为什么要做。** SSP 情景在 **2015 年就分叉**了：1750–2014 是 CMIP6 历史，2015 起全是
投影。而本模型的损害起始年是 2025，也就是说 2015–2024 这十年在模型里本该是「已经发生
的事实」，用的却是十年前写下的预测。实测偏差：

| 事实 | 原 SSP2-4.5 | 实测清单 |
|---|---|---|
| 2020 年 CO₂（相对 2015） | **+3.8%** | **−0.2%**（EDGAR）/ **−0.8%**（CEDS） |
| 2024 年 CO₂（相对 2015） | +6.7% | +9.2% |
| 2015–2020 年 CH₄ | 一条水平线（388.07 → 388.09 Mt） | +2.8% |
| 2020 年航空 NOx（相对 2015） | 情景无变化 | **0.650**（−35%，停摆） |

新冠停摆在情景里完全没有痕迹。后果可验证：校正前模型 2024 年 CH₄ 浓度比 NOAA 观测
**低 21.3 ppb**。

**数据源与分工。** 判据是**回接检验**——把接续法倒过来用在 2015 年以前（那里 RCMIP 有
真值），直接量出方法本身的误差；2005–2014 窗口的平均绝对百分误差：

| 物种组 | 数据源 | 回接误差 | 备注 |
|---|---|---|---|
| CO₂ / CH₄ / N₂O | **EDGAR-2025** | 1.57% / 0.73% / 2.44% | CEDS 分别是 1.40 / 2.74 / 3.02%，CO₂ 打平、CH₄ 与 N₂O 输给 EDGAR；EDGAR 原生覆盖到 2024 |
| SO₂/NOx/CO/NMVOC/BC/OC/NH₃ | **CEDS v_2025_03_18 + GFED5.1** | 1.0–10.8% | 唯一来源。CEDS 止于 2023，最后一年沿用；GFED 覆盖到 2024 |
| 航空 NOx（`nox_avi`） | **CEDS 航空部门** | — | 单列处理，见下 |
| 15 个含氟物种 | **EDGAR-2025 F-gases** | 3.8–67% | 唯一来源。SF₆ 3.8%、HFC-134a 6.4% 最好；HFC-23 67% 最差，但绝对量对 SCC 无影响 |
| 18 个蒙特利尔气体 + 8 个 EDGAR 不报的含氟 | **不动** | — | 两个清单都不报。RCMIP 里它们本就来自《蒙特利尔议定书》的大气观测反演，比自下而上清单可靠 |

**接续方法：只换增长率，不换水平。** 清单口径与 RCMIP 口径不同（EDGAR 的 CO₂ 不含土地
利用、CEDS 不含生物质燃烧、CEDS 2025 版相对 CMIP6 用的 2017 版有系统性修订），直接替换
水平会引入口径断层。故只取相对变化：

$$E[y]=E^{\text{rcmip}}_{2015}\cdot\frac{\text{INV}[y]}{\text{INV}_{2015}},\qquad 2016\le y\le 2024$$

气溶胶与前体物再拆两个分量，因为 RCMIP 的总量 = 人为 + 开放生物质燃烧，而 CEDS 只有
人为源（CO 的人为份额 0.81、OC 0.68）。把人为源的增长率施加到含生物质燃烧的总量上，
等于假设森林火灾同步下降——2023 年加拿大大火令 GFED 的黑碳比常年高 40% 以上：

$$E[y]=E^{\text{rcmip}}_{2015}\Big(sh\cdot\frac{\text{CEDS}[y]}{\text{CEDS}_{2015}}+(1-sh)\cdot\frac{\text{GFED}[y]}{\text{GFED}_{2015}}\Big),\quad sh=\frac{\text{CEDS}_{2015}}{E^{\text{rcmip}}_{2015}}$$

**2025 及以后：以 2024 观测值为起点、按情景自身的增长率续接。**

$$E[y]=E^{\text{obs}}_{2024}\prod_{t=2025}^{y}\frac{E^{\text{ssp}}[t]}{E^{\text{ssp}}[t-1]}=E^{\text{ssp}}[y]\cdot\frac{E^{\text{obs}}_{2024}}{E^{\text{ssp}}[2024]}$$

连乘望远镜式相消，**等价于一个逐物种、逐情景的常数倍率**。因此 2024→2025 没有跳变，
情景的轨迹形状原样保留，只是整条路径被重基准到实测的 2024 年水平。倍率逐情景不同
（五个情景的 $E^{\text{ssp}}[2024]$ 已经分叉）：

| 情景 | CO₂ | CH₄ | SO₂ |
|---|---|---|---|
| ssp119 | 1.2974 | 1.3421 | 1.2790 |
| ssp126 | 1.1318 | 1.2993 | 1.1544 |
| ssp245 | 1.0234 | 1.0683 | 0.9921 |
| ssp370 | 0.8903 | 0.9542 | 0.7943 |
| ssp585 | 0.8844 | 0.9981 | 1.0313 |

> **引用时怎么说**：情景不再是纯粹的 SSP2-4.5，而是**「以实测 2024 年排放重基准后的
> SSP2-4.5」**。这是 IAM 谐调化（harmonisation）的标准做法，但必须写明，因为 2100 年的
> 排放水平确实被整体缩放了。

**乘法重基准的两个后果，低排放情景尤其明显。**

*净零年份被精确保留*——常数倍率不移动零点。这是选乘法而非加法偏移的实质理由：加法会把
SSP1-1.9 的净零年推后若干年，而净零时点是这些情景最核心的政策含义。

*但负排放也被同比例放大*。SSP1-1.9 假设了 2015 年起的快速脱碳，现实没发生，于是它的
2024 年情景 CO₂ 只有 8.9857 GtC，而实测是 11.658，倍率 **1.2974**。这个倍率同样作用在
它的负排放段上：2100 年由 −3.7873 变成 −4.9136 GtC。也就是说，**模型里的 SSP1-1.9 隐含
了比原情景多 30% 的碳移除**。ssp126 同理（倍率 1.1318）。ssp245 只有 1.0234，影响可忽略，
而生产结果引用的正是 ssp245——但用本包跑 ssp119/ssp126 时必须知道这一点。

各情景 2024 年的情景 CO₂（GtC）与由此得到的 carry 倍率，实测统一为 **11.6581**：

| 情景 | 2024 情景值 | 倍率 |
|---|---|---|
| ssp119 | 8.9857 | **1.2974** |
| ssp126 | 10.3003 | 1.1318 |
| ssp245 | 11.3917 | **1.0234** |
| ssp370 | 13.0944 | 0.8903 |
| ssp585 | 13.1823 | 0.8844 |

高排放情景的倍率小于 1：SSP3-7.0 与 SSP5-8.5 假设的 2024 年排放比实际**高** 12–13%。

**效果与代价（ssp245，2030 脉冲，积分至 2100，`end_year=2150`）**：

| 量 | 校正前 | 校正后 |
|---|---|---|
| 2024 年 CO₂ 浓度误差 vs NOAA | +0.76 ppm | +0.95 ppm |
| 2024 年 CH₄ 浓度误差 vs NOAA | **−21.33 ppb** | **+10.85 ppb** |
| 2024 年 N₂O 浓度误差 vs NOAA | +0.75 ppb | +0.77 ppb |
| 三气体相对误差合计 | 1.5124% | **1.0183%** |
| **SCC 点估计** | 687.3172038776413 | **690.4507205341098**（+0.456%） |

CO₂ 浓度略微变差是预期内的：FaIR 的碳循环本来就偏高，补上（正确的、更高的）实测排放
会把它推得更高。那 0.95 ppm 残差在碳循环参数上，不在排放清单上。

**分项贡献**（相对 687.3172）：航空 NOx **−0.233%**（比全部 15 个含氟气体加起来还大）；
含氟合计 −0.16%；其余为三种主要温室气体与气溶胶。

**已知局限。** 12 个物种（`ch2cl2`、`methyl_chloride`、`chcl3`、`methyl_bromide`、
`hfc365mfc`、`so2f2`、`hfc4310mee`、`c4f10`–`c8f18`）EDGAR 与 CEDS 都不报，2016–2024
**仍是情景值**，跨情景离散 24–54%。绝对量最大的是 `ch2cl2`（2024 年跨情景极差 0.615）。
生成脚本会把这份清单逐条打印出来，不做静默处理。

**闸门**（`03Claim/scripts/test_fair_emissions_integrity.py`，7 条）：前言物理行数与表头
不变（Julia 侧 `skiplines_begin=6` 是写死的行数，多一行少一行会把表头当数据读且不报错）、
1750–2015 与出处逐位相同、2016–2024 的 26 个已替换物种**跨五情景逐值相同**（观测不能在
过去分叉，这条比原来的 1750–2014 更强）、未替换物种全时段与出处相同、接缝处增长率与原
情景一致、三条方向自检、五个文件的 MD5。

### 4.2 MESMER 统计降尺度（`climateregional`）

`src/components_climate/ClimateRegionalComponent.jl`。**2026-07-30 起分两段构造**（`open_issues.md` P1-16）：

| 段 | 年份 | 国家温度 |
|---|---|---|
| 历史 | 1950–2024 | **直接使用观测** `climate_regtmp_observed.csv`（Berkeley 人口加权，与 LP 的 `lctmp_bkly_pw` 同源） |
| 前向 | 2025– | MESMER 斜率 + **观测锚定的水平**：$\beta_c(T^{global}_t - \bar T^{global}) + \bar o_c + \nu_{c,t}$ |

此前历史段也走 MESMER 重构，而观测数据本就存在。实测重构误差：逐国 RMSE 中位 0.487 °C，**89/185 国超过 0.5 °C**，系统性偏冷 −0.182 °C（FaIR 的历史全球升温低于 Berkeley 观测，加上 MESMER 的 γ 截距）。

$\bar T^{global}$ 与 $\bar o_c$ 是 FaIR 全球温度与观测国家温度在 **1995–2024** 上的均值，在拼接年现场算定。**水平与参考点必须取同一纪元**，两者的纪元偏移才相互抵消；该式与"固定斜率 β 的 OLS 截距拟合"代数上完全等价。

不用拼接年单年值做锚，是因为观测 = 趋势 + 天气：2024 是破纪录暖年、高出趋势 +0.368 °C，单年锚会把这个异常永久固化进前向趋势线的水平，而模型还会另外抽 ν 叠加——同一年的天气算两遍。实测长期偏置：30 年均值锚 0，单年锚 +0.368 °C。

连带后果：MESMER 的 `region_lt_inter`（γ）与 `climate_tempbase` 在前向段不再被使用，水平改由观测拟合（标准的 delta-change 偏差校正：保留 β 的敏感度模式、用观测校正水平），两者保留为诊断与 legacy 分支所用。

以下分解式描述的是**前向段**；MESMER（Beusch et al.）把全球温度模拟成区域温度 = 线性趋势响应 + 随机天气变率：

$$T^{reg}_{c,t} = \underbrace{\beta_c (T^G_t - T^{G}_{1950}) + \gamma_c}_{ts\_lt\text{（趋势项）}} + \underbrace{tempbase_c}_{\text{1950 局地基准温度}} + \underbrace{\nu_{c,t}}_{region\_lv\text{（天气噪声）}}$$

- `region_lt_coef`（β_c）/ `region_lt_inter`（γ_c）：国家级趋势响应系数，由网格参数经人口权重矩阵投影而来（见 §4.2.1）；
- `region_lv`：国家级天气变率实现，由 `helper_mesmer.jl` 在网格上模拟后聚合，带显式 seed；
- `T_global_ref`：当前为 1950 年 FaIR 参考温度 `0.2749925941187845`，使 MESMER 趋势响应相对于模型起点计算；
- 输出 `regtmp`（区域绝对温度）供损害模块使用。

#### 4.2.1 从网格到国家：人口加权直投（2026-07-30 重构）

旧链路把 MESMER 的陆地格点先平均到 2,388 个次国家 `RegionID`、再按**面积**加权到 185 国。两层平均都与回归端 `lctmp_bkly_pw` 的**人口**加权口径不一致，且映射不到的国家被静默补 0——实测导致 6 国 β=γ=0（前向模型里永不升温）、26 国 `region_lv` 恒零。现已取消 `RegionID` 中间层，改为：

```
mesmer_data/export_grid_params.py     训练，持久化网格参数（5,940 陆地/岛屿格点）
      │   climate_grid_lt.csv          β、γ
      │   climate_grid_lv.csv          coef_gvtas + AR(1) 三件套
      │   climate_grid_gv.csv          AR(2) 全球变率
      │   climate_grid_innov_chol.csv  新息协方差的稀疏 Cholesky 因子（约 146 MB）
      │   climate_grid_geometry.csv    land_idx / lat / lon / land_fraction
      ▼
scripts/build_mesmer_country_weights.py   人口权重矩阵 W（185 × 5,940）
      ▼
src/helper/helper_mesmer.jl               网格模拟 → W 聚合 → 三个国家级 CSV
```

**模拟域**：MESMER 的 `extract_land(threshold_land=1/3)` 是其自身的过滤而非数据限制。原始 CMIP6 场是全球 20,592 格点，含任何陆地成分的有 5,935 个（已剔南极），1/3 门槛只留 5,121 个——斐济自己那格陆地占比 0.30，恰好被卡掉。现改用 `threshold_land=0` 并补入 SYC/BMU/MUS/MDV/CPV 五个在此分辨率下完全不算陆地的孤岛格点，共 **5,940** 点。保留全球 20,592 点不可行：新息协方差是 $n^2$，会达 3.4 GB，`load_phi_gc` 更需 35×3.4 GB。

**局地变率的模拟顺序不可颠倒**。MESMER 的局地变率是 AR(1) + **空间相关新息**：

$$\nu_{i,t} = c^{gv}_i \, gv_t + a_{i,t}, \qquad a_{i,t} = \text{AR1\_int}_i + \text{AR1\_coef}_i\, a_{i,t-1} + \varepsilon_{i,t}, \qquad \varepsilon_t = L z_t$$

其中 $L$ 是局地化新息协方差 $\Sigma$ 的 Cholesky 因子（$\Sigma = LL^\top$），把独立正态 $z$ 转成空间相关的扰动——真实热浪覆盖整片区域，逐格独立抽样会得到物理上荒谬的白噪声，聚合到国家后方差还会被严重低估。必须**先在网格上模拟、再聚合**：若先把参数聚合到国家再模拟，(a) 国家间天气相关性会丢失，(b) 各格点 AR(1) 系数不同，其加权和数学上是 ARMA 而非 AR(1)，自相关口径会变。

**局地化半径的来源**。$\Sigma$ 不能直接用样本协方差：训练样本只有 251 个（$h$-ssp126 的 251 年 × 1 个 run），而维度是 5,940，样本协方差秩至多 250，远距离项几乎全是噪声。MESMER 因此把样本协方差与 Gaspari–Cohn 函数做 Schur 积 $\hat\Sigma \odot \varphi_{GC}(d_{ij}/L_{loc})$，既压掉远距离噪声，又靠 Schur 积定理保住正定性（$\varphi_{GC}$ 在 $d = 2L_{loc}$ 处严格归零）。$L_{loc}$ 由 15 折交叉验证在候选半径上按留出折的负对数似然选出，当前值 **1,250 km**，是 `750:250:2000` 窗口内的内点解（见 `open_issues.md` P1-14 与 `mesmer_data/localisation_radius_scan.csv` 的完整 $\mathrm{nll}(L)$ 曲线）。注意此处的 $L_{loc}$ 与上式中作为 Cholesky 因子的 $L$ 是两个不同的量，勿混。

**与 Python 端不可逐位复现**：numpy 的 `multivariate_normal` 默认用 SVD 分解协方差并使用自身 RNG，Julia 端用 Cholesky + `MersenneTwister`。两者分布相同但具体实现不同，旧的 300 个 emu 文件不可能被重现；等价性在分布层面检验（方差、自相关、空间相关、与观测的量级比）。复现性由显式 seed 保证，seed 写入 `climate_region_lv_manifest.csv`。

### 4.3 社会经济模块（`socioeconomic`）

`src/components_socio/SocioEconomicComponent.jl`，刻意保持极简的外生基准：

$$ypc_{c,t} = ypc_{c,t-1}(1 + ypcgrowth_{c,t}), \qquad income_{c,t} = ypc_{c,t} \cdot pop_{c,t}$$

输入 `ypc0`（初始人均收入）、`ypcgrowth`（增长率路径）、`population`。**基准路径不受气候影响**——损害以偏离基准的方式在下游损害组件中体现。

### 4.4 B&K 损害模块（`claim_damage`，本项目主线）

`src/components_damage/claim_damage.jl`，实现 §3.4 的全部方程。每个时间步、每个国家依次执行：

| 步骤 | 计算 | 依据 |
|------|------|------|
| 1 | 状态变量：$T^{state}$ 取过去 30 年 `regtmp` 滚动均值，$Y^{state}$ 取滞后一期收入（均按 LP 口径中心化） | §3.4.1 |
| 2 | 动态 θ 投影：$\hat\theta_h(t) = \theta_0 + \theta_1 T^{state} + \theta_2 Y^{state}$（G/L 各 11 个 h） | §3.4.1 |
| 3 | 闭式 OLS（`helper/helper_damage.jl::closed_form_ols`）：θ → (A, C)，用 Y 基底 | §3.3.4 |
| 4 | 组装当年损害核 $\zeta_s = A Z_1[s] + C Z_2[s]$，用 Z 基底 | §3.3.4 |
| 5 | TFP 卷积：全球渠道用 $T^{global} - T_{ref}$，局部渠道用 $\nu$ | §3.4.2 |
| 6 | Solow 资本递归 $\hat k_t$（状态变量，跨期传递，永不截断） | §3.4.2 |
| 7 | $\exp$ 映射回 ypc_new、loss；记录 clipped 诊断 | §3.4.3 |

参数分四组：气候输入（T_global, T0_global, regtmp, region_lv_local, global_lv）、经济输入（ypc_baseline, population）、预计算缓存（8 个基底矩阵 + basis_model + w_g/w_l）、LP 系数（6 个 θ 矩阵，per region × horizon）与结构参数（alpha_c, delta_c, g_c）。

### 4.5 BHM 对比损害模块（`bhm_damage`）

`src/components_damage/bhm_damage.jl` 实现 Burke/Hsiang/Miguel 二次绝对温度增长损害，仅用于论文方案对比。年度损害为 $f(T_{climate}^{abs})-f(T_{baseline}^{abs})$，其中 $f(T)=\beta_1T+\beta_2T^2$；组件累积该增长效应到 `z_tfp`，并固定 `khat=0`、`yhat=z_tfp`，不传播资本响应。参数来自 `data/unshared_parameters/bhm_quadratic_params.csv`，当前不对 BHM 系数做不确定性抽样。

### 4.6 SCC 计算模块（`helper/helper_scc.jl`）

国家级损害组件共用的 SCC 模块（取代已删除的 MimiFUND 移植文件 `new_marginaldamages.jl`）：

- `add_co2_pulse!(m, year; pulse_size, pulse_length)`：把 `emissionspulse` 加法器插在 `emissions` 与 `co2_cycle` 之间，pulse_size 吨 CO₂ 平摊 pulse_length（默认 10）年，换算为 GtC/yr（/1e9 × 12/44）；仅支持 CO₂ 通道；
- `get_marginal_model(m; year, pulse_size, pulse_length)`：`create_marginal_model` 复制出 Pulse 轨并注入脉冲；
- `compute_scco2(m; year, last_year, eta, prtp, pulse_size, pulse_length, return_details)`：先运行 Base；B&K 模型会把 `T0_global` 的 NaN 哨兵解析为 `damage_start_year` 当年 FaIR 温度并覆写两轨（Base/Pulse 共用参考温度，§3.6），BHM 模型不需要该步骤；随后对活动损害组件（`claim_damage` 或 `bhm_damage`）的双轨 `loss` 做差分（×1e6，population 为 PWT 百万人口径），并按 Ramsey 因子（prtp + η·g，g 用含气候损害的 Base 轨 `ypc_new` 增长率，DICE 式内生折现）逐国折现加总；
- 返回标量 SCC（\$/tCO₂，2017 PPP 美元）；`return_details=true` 时返回 `(scc, scc_regional, marginaldamage, df, mm)`。

### 4.7 Monte Carlo 机制

旧 FUND 式 `src/montecarlo/`（运行时分布抽样）已删除。当前 national 模式的损害参数不确定性由 `damage_theta_base.csv` / `damage_theta_draws.csv` 预抽样驱动——外部循环 `get_model(theta_id=i)` 后逐次调用 `compute_scco2` 即可。

外层不确定性由 `helper/helper_mc.jl::run_scc_mc` 驱动，一次抽样同时抽四段（①②③ 属气候层，④ 属损害层）：

| 段 | 载体 | 应用函数 | 开关 |
|---|---|---|---|
| ① 全球平均温升响应 | `climate_fair_ensemble.csv` 的成员（热响应 q/d、CO₂ 强迫 f、碳循环 r、气溶胶 ari/aci） | `helper_fairv2.jl::apply_fair_member!` | `sample_fair` |
| ② 降尺度的**升温模式** | `climate_ltcoef_ensemble.csv` 的 58 个 CMIP6 模式（国家级 β） | `helper_mc.jl::apply_region_lt!` | `sample_esm` |
| ③ 降尺度的**内部变率** | MESMER 局地变率实现，逐抽样换 seed | `helper_mc.jl::apply_region_lv!` | `sample_mesmer` |
| ④ **损害参数 θ** | `damage_theta_draws.csv` 的 pooled 系数抽样（生产：300 条 year-block bootstrap replicate） | `helper_mc.jl::apply_theta!` | `sample_theta` |

> **④ 于 2026-08-06 加入（§6.3 H2 / R24）**。此前 θ 只能靠外部循环 `get_model(theta_id=i)` 重建模型来换，于是气候层与损害层永远各跑各的，拿不到联合区间。改为原地 `update_param!` 后，逐抽样换 θ 的代价可以忽略（`get_model` 要 2.84 s，`run` 只要 0.17 s）。
>
> **四个维度各用独立 RNG，这是设计要求而非风格**：只有这样，关掉任一维时其余维的抽样序列才原样不动。由此三个口径共用同一串随机数（**共同随机数，common random numbers**），`joint − climate` 干净地只反映 θ 层、`joint − damage` 只反映气候层，抽样噪声大部分相消。独立跑三次就没有这个性质。
>
> **`mode` 标签由开关导出**（`_mc_mode_label`），不是独立参数——设成参数就会出现"标签写着 joint、开关只开了气候"这种不报错的自相矛盾状态。默认输出路径随之为 `results/scc_mc_<mode>.csv`，三个口径互不覆盖。
>
> 三个口径 + 点估计构成 2×2 因子实验，用 `julia 03Claim/scripts/run_uncertainty_matrix.jl [n] [scenario]` 一次跑齐。
>
> **④ 与 ①②③ 的实现方式刻意相反**：θ 是**预抽样**（Python 侧 numpy PCG64 生成，Julia 只按 `draw_id` 取），③ 是**现场模拟**。理由见 §3.5.1 与 `.claude/rules/claim-data.md`：θ 只有 6 个数且需要跨语言逐位复现与可审计；③ 是 185×551 的大矩阵，池化会退回 P1-3 之前的病灶。
>
> **口径限制**：④ 只扰动 6 个 **pooled** 系数，国别 $B$、$D$ 仍冻结（R23/H5）。因此"联合区间"**仍不是完整不确定性区间**。

②③ 是 `beusch2020mesmer` 明确区分的两层：模式间散布在年尺度上 "clearly exceeds" 单模式初始条件散布。逐国 β 的跨模式变异系数中位 10.8%、p90 15.3%，确实远大于 ③ 的 seed 维度。三个开关各自独立，可单独关掉任一维做归因。

> **但对全球加总的 SCC，② 只值 0.1% 的方差**（800 次三维抽样实测，`open_issues.md` 结构性 caveat #7）。原因分两段、已完全查清：**(a)** 跨模式差异主要是**再分配**——逐国 CV 10.8%，但按损害份额加权的 β 汇总量跨模式散布只有 6.4%；**(b)** SCC 对 β 的弹性只有 **0.102**（0.8×–1.5× 上完美线性）。两段相乘 0.65%，与实测的 0.57% 一致。
>
> 弹性 0.102 的根因不是"衰减"，而是 **β 根本不乘在损害冲击上**。`claim_damage.jl` 里驱动损害的是 `T_global + global_lv`（全球渠道）与 `region_lv_local`（局地渠道）——国别温度 `regtmp` 从不作为冲击进入损害函数，β 只经 `T_state` 调制弹性 θ。这是 B&K 的设定本意：宏观损害由全球温度冲击驱动，国别异质性来自"多热、多穷"的**状态**。实测：β 放大 1.5 倍时脉冲的全球温度比为 **1.0000**（冲击不变），国别温度比 1.5030（但损害不读它），TFP 冲击比只有 1.0643。也不是 Base/Pulse 相消——BAU 累计损害的 β 弹性同样只有 0.073。
>
> **结构性事实**：降尺度产物只以**状态变量**的身份进入损害（决定各国损害弹性），不以冲击的身份进入。这同时解释了 ν 为何只值 0.18%。
>
> **引用口径**：模式选择的影响分三层，实测值各不相同，不可互相替代——**全球 SCC 方差 0.1%**；**国别损害** CV 中位 0.80%、p90 2.47%（最高 IRL 13.3%）；**国别温度产出**跨模式极差中位 **1 °C**、最高 3.16 °C（ISL）。即模式选择对国别温度是一阶的，对国别损害不是。反直觉的对照：国别损害对**天气实现 ν** 的敏感度（p90 21.65%）反而大于对模式选择的敏感度——因为 ν 经 `z_l` 作为**冲击**进入损害，而 β 只作为**状态**进入。

① 的 FaIR 成员**必须随机置换后取，不能取集合的前 $n$ 行**（`fair_sampling=:random`，`fair_seed` 默认 20260805；`open_issues.md` §6.2 R7）。`climate_fair_ensemble.csv` 是 `build_fair_constrained_ensemble.py` 中 `systematic_resample` 的产物：2000 行等权，**后验权重由重复次数编码**，行序即源先验的生成顺序。因此任何连续切片都是聚簇样本而非子样本——

- 前 300 行只覆盖 **119 个不同**的源成员，随机 300 行平均覆盖 **241 个**（200 次重抽范围 226–255），**有效样本量少一半**；
- 分布上也可区分：前 300 行对其余 1700 行的两样本 KS，ECS $p=2.3\times10^{-4}$、TCR $p=1.7\times10^{-5}$、$r_T$ $p=1.2\times10^{-7}$，且三者**均偏暖**。

**按行等概率抽，不要去重**：重复正是后验权重的编码方式，去重等于把加权后验拍回等权先验。`fair_sampling=:prefix` 保留旧行为，**仅用于复现 2026-08-05 之前的产物，不可用于生产**。输出表新增 `source_member_id` 列，"抽到了几个不同的气候参数集"只能看这一列——`member_id` 是重采样后的行号，逐行互异，看不出重复。

② 的 ESM 用**独立 RNG**（`esm_seed`，默认 20260803）抽取，不复用抽样序号：这样开关 `sample_esm` 不会改变 FaIR 成员与变率 seed 的配对。①②③ 三个维度各用独立 RNG，正是为了让任一维单独开关时其余配对不变。

**默认按模式独立性权重抽样**（`esm_weighting=:independence`）。等权抽 58 个模式会放大 CMIP6 的模式家族——EC-Earth3 五个变体、E3SM 四个，CanESM5/CESM2/GISS/CNRM 各三个——等于让同一个模式投多票。权重用 Knutti et al. (2017, GRL 44:1909)：

$$w_i \;\propto\; \frac{1}{1 + \sum_{j \neq i} \exp\!\left[-\left(S_{ij}/\sigma_S\right)^2\right]}$$

分母是模式 $i$ 的**有效近邻数**：孤立的模式得 $\approx 1$，五个变体扎堆的模式各得 $\approx 1/5$。

- $S_{ij}$ = 两个模式在 2957 个 g025 陆地格点上 **β 场**的面积加权 RMS 差。用 β 而不是外部诊断量，是因为我们传播的就是 β——独立性度量与被传播的量同源，"这两个模式对我们的结论而言是不是同一个模式"才有确切含义。
- $\sigma_S$ = "同一模式不同配置"成对距离的中位数（`SAME_MODEL_GROUPS`），逻辑同 Merrifield et al. (2020)：半径应当就是"同一个模式内部"的尺度。该分组表由公开的 CMIP6 `source_id` 命名决定，且刻意保守——仅共享机构或部件的配对不入组，这会缩小 $\sigma_S$、让合并更不激进。
- **只做独立性加权，不做性能加权**。后者需要对升温模式本身的观测约束，是另一件事，见 `open_issues.md` 结构性 caveat #7。
- `esm_weighting=:equal` 保留等权配置作为对照。权重与 58×58 距离矩阵分别落在 `climate_ltcoef_ensemble_weights.csv` 与 `climate_ltcoef_ensemble_distance.csv`。

③ **每次抽样现场在网格上模拟**（`simulate_region_lv`），不使用预生成的实现池——池化会退回 P1-3 之前「读现成 emu 文件」的模式。实测约 3.3 s/次抽样（三维），结果逐次落盘。**单独跑模型时不调用 ②③**，此时用的是自训练的 IPSL `climate_ltcoef.csv` 与 `climate_region_lv.csv` 中 seed=42 的实现，可复现。

> **两套网格并存，不可按下标混用**。② 走 MESMER v1.0.0 发布的 g025 网格（2957 陆地格点），③ 仍走 `export_grid_params.py` 自训练的 IPSL 原生网格（5940 格点）。原因是逐 ESM 复现变率需要 2957×2957 的局地化新息协方差（58 模式合计 2.66 GB），已按 P1-7 的实测（该渠道对 SCC 影响上界约 0.2%）决定不下载。故 ② 与 ③ 各有自己的权重矩阵：`climate_grid_weights_g025.csv` 与 `climate_grid_weights.csv`。`land_idx` 是位置下标不是坐标，两套混用不会报错、只会静默算错。

`apply_fair_member!` 放在 `helper_fairv2.jl` 而非构建脚本里，是因为它被**两条路径共用**——构建集合时给候选打分、消费集合时写进完整模型。两边若各留一份拷贝，「集合成员」在两个阶段可能悄悄变成不同的东西且不报错。

集合为 `unshared_parameters/climate_fair_ensemble.csv`（2000 成员），由 `03Claim/scripts/build_fair_constrained_ensemble.py` 生成：CMIP6 逐模式参数 → Liu–West 平滑自助先验 → 全记录四段升温增量 + 地球能量收支约束 → 加权重采样。方法与全部诊断见同目录的 `climate_fair_ensemble_method.md`。

> **使用规定（`open_issues.md` P1-7）**：气溶胶与外部强迫的不确定性未抽样，气溶胶/敏感度简并被人为切断，后验 ECS 因而被系统性压低、TCR 的 90% 区间仅为 AR6 *likely* 区间的 0.66 倍。**传播到 SCC 时必须把本集合的区间作为气候不确定性的下界报告**，并以 28 个 CMIP6 模式的散布作为上界参照，不得单独引用。

---

## 5. 文件逐一说明

### 5.0 数据文件分类：源数据 vs 生成物

改动任何文件之前，先确认它属于哪一类——**生成物不要手工编辑**，否则下次重跑管线会被覆盖，且破坏可复现性：

| 类别 | 路径 | 是否手工编辑 | 说明 |
|------|------|-------------|------|
| 原始/外部数据 | `01data/`、`02regression/02data/`（pwt110.xlsx 等）、`03Claim/data/fair_data/`、`03Claim/data/mesmer_data/` | ❌ 不手改 | 原始观测数据或外部模型标定参数 |
| 回归产物 | `02regression/03table/*.csv`、`02data/struct_params.csv`、`analysis.dta` | ❌ 不手改 | 由 Stata/Python 管线生成，改上游脚本后重跑 |
| B&K 预计算缓存 | `03Claim/data/unshared_parameters/damage_*.csv` | ❌ 不手改 | 由 `precompute_*.py` 生成 |
| national 输入 | `03Claim/data/shared_parameters/*.csv` 与 `03Claim/data/unshared_parameters/*.csv` | ⚠️ 原则上由脚本生成 | 非 SSP 输入由 `prepare_national_inputs.py` 产出；SSP 人口与增长路径由 `01data/build_wb_2024_actuals.py`（2024 过渡年）、`01data/ssp/01_extract_ssp.py` 和 `02_build_scenarios.py` 产出，CLAIM 只消费结果 |
| 源码 | `03Claim/src/*.jl`、`03Claim/scripts/*.py`、`02regression/01scripts/` | ✅ 可修改 | 模型逻辑与预处理逻辑 |
| 文档 | `coupling_plan.md`、`open_issues.md`、本文档 | ✅ 可修改 | 方法、审计与用户说明 |

### 5.1 `03Claim/src/`

| 文件 | 说明 |
|------|------|
| `MimiCLAIM.jl` | **测试/脚本库入口**。`get_model(:national; climate=:fair/:csv, scenario, damage_function, income_state, capital_response, theta_id, start_year, end_year, bau_ref_year, TCR, RWF, F2x)`：默认 `climate=:fair`（1950–2150，FaIR 内生升温），默认损害机制为 `damage_function=:BK_region, income_state=:dynamic, capital_response=true`；`scenario` 仅限 `ssp119`/`ssp126`/`ssp245`/`ssp370`/`ssp585`，同时选定 FaIR 排放路径与 IIASA SSP v3.2 社会经济基线（ssp119/ssp126→SSP1，ssp245→SSP2，ssp370→SSP3，ssp585→SSP5）。用户交互式 SCC/出图主入口是内联完整 `get_model` 的 `src/MimiCLAIM.ipynb` |
| `components_damage/claim_damage.jl` | B&K 损害核 Mimi 组件（§4.4） |
| `components_damage/bhm_damage.jl` | BHM 二次绝对温度增长损害对比组件（§4.5） |
| `helper/helper_global.jl` | 全局层：通用 CSV 读取器（`_load_time_region_matrix` 等，注意其"向前平推最后一个值"的填充语义）、维度/区域、组件装配与连线、shared 参数注册、socio 与 climateregional 参数加载（`update_global_params!`） |
| `helper/helper_damage.jl` | 损害模块专用：闭式 OLS、basis/theta/state_centers 加载器、`theta_parameters`（把 pooled draw 复制到各国）、`update_claim_damage_params!` |
| `helper/helper_scc.jl` | SCC/边际损害计算（§4.6） |
| `helper/helper_fairv2.jl` | FaIR v2.0 参数批量加载 `update_MimiFAIR2_params!`（构造时从默认 CSV 装载全部参数并完成组件连线）；以及 `apply_fair_member!`——把 `climate_fair_ensemble.csv` 的**单个成员**写进已建好的模型。两者职责不同：前者装默认值，后者设成员值，且后者被构建脚本与外层 MC 共用（§4.7） |
| `helper/helper_mc.jl` | 外层不确定性 MC：`run_scc_mc`（FaIR 成员 × ESM 升温模式 × MESMER seed × 损害 θ 四维联合抽样，逐次落盘；`run_climate_mc` 为强制 `sample_theta=false` 的历史别名）、`apply_region_lt!` 与 `load_lt_ensemble`（覆写国家级 β，58 模式集合）、`apply_region_lv!`（内存内覆写变率，含 Base/Pulse 共用 ν 的不变量）、`apply_theta!`（内存内覆写 pooled θ）、`_mc_mode_label`（由开关导出口径标签）、`draw_region_lv`（§4.7） |
| `helper/helper_mesmer.jl` | MESMER 网格参数加载、局地变率模拟（复刻上游 `OLS_AR1_sci`）、人口加权聚合到 185 国、`regenerate_national_mesmer_inputs`，以及出图用的 `grid_warming`（§4.2.1）。旧的 `load_region_lt`/`load_region_lv` 保留在文件末尾"存档"一节，不参与运行 |
| `components_socio/*.jl` | socioeconomic 组件（§4.3） |
| `components_climate/*.jl` | FaIR 组件 + ClimateRegional/MESMER 组件（§4.1–4.2） |
| `components_damage/*.jl` | B&K 主损害组件与 BHM 论文对比组件（§4.4–4.5） |

### 5.2 `03Claim/scripts/`（预处理：生成模型输入）

> **分层判据**：`src/` 装模型**运行时**需要的东西，`scripts/` 装**生成模型输入**所需的东西；依赖方向只能 `scripts → src`。据此，只服务于集合构建的 `build_fair_only_model`、`model_ohc_change`、似然机制等留在 `scripts/`；而 `apply_fair_member!` 因构建与消费两边共用，必须下沉到 `src/`。
>
> 本目录以 Python 为主，另有三个 Julia 文件：`generate_fair_1950_initial_conditions.jl`（一次性生成 FaIR 1950 初始条件）、`fair_likelihood.jl`（给候选参数组打分，是 `build_fair_constrained_ensemble.py` 的后端）及其测试。

| 文件 | 输入 → 输出 | 说明 |
|------|------------|------|
| `build_observed_emissions.py` | `01data/emissions/{edgar,ceds,gfed}` + 原始 RCMIP → `data/fair_data/rcmip_ssp*_..._rebased_obs.csv`（5 个） | **模型排放输入的观测校正**（§4.1.1）：2016–2024 换实测、2025+ 按 2024 实测水平重基准。`--fair_dir` 把产物导向别处（clean-output 测试用），出处始终读生产目录的原始 RCMIP。改了上游清单版本就要重跑，并更新 `test_fair_emissions_integrity.py` 的 `EXPECTED_MD5_OBS` |
| `precompute_basis.py` | damage_params + struct_params + temp_irf + **lp_imputed_boot** → `data/unshared_parameters/damage_basis.csv` | 每国每渠道算 Z₁/Z₂（纯核）、Ztilde（φ 卷积后）、Y₁/Y₂（Solow 传播后）双基底 + 逆方差权重。**`--lp_imputed` 必须与拟合 `damage_params.csv` 时用的是同一份**——它决定 `damage_weights.csv`，给错不会报错，只会让基底与参数不同源 |
| `precompute_cholesky.py` | lp_main → `data/unshared_parameters/damage_cholesky.csv` | 组装 11 个 6×6 Σ_h（含跨渠道块），robust Cholesky（对称化 + jitter） |
| `precompute_mvn_draws.py` | lp_main → `data/unshared_parameters/damage_theta_base.csv` + `damage_theta_draws.csv` | 因子化平滑抽样：1000 × Z⁽ᵐ⁾ 共享展开为 (1000, 11, 6) θ 数组。**历史口径**——生产的 `damage_theta_draws.csv` 现在由下一行的脚本产出 |
| `build_theta_draws_bootstrap.py` | `lp_bootstrap_bl10.csv` → `damage_theta_draws_bootstrap.csv` | **生产口径的 Σ_θ**：300 条 replicate 重新居中到点估计（只取离散度不取位置），块长只允许 10/15。复制成 `damage_theta_draws.csv` 才生效，这一步是人做的决定 |
| `prepare_national_inputs.py` | struct_params + analysis.dta + berkeley_full.dta + **`--temperature_cache`（必需）** → `shared_parameters/*.csv` + `unshared_parameters/*.csv` | 生成 185 国维度的非 SSP 模型输入；`climate_T_global.csv` 直接取 `berkeley_full.dta::berk_anom`；`climate_regtmp_observed.csv`（1950–2024 观测国家温度，§4.2）由 `build_regtmp_observed` 从 `country_temperature_berkeley_grid_pw.csv` + `TEMPERATURE_PROXY` 生成。**`--temperature_cache` 自 2026-08-06 起为必需参数**（`open_issues.md` R19）：漏传曾只是不写观测文件，模型端 `@warn` 一句就退回 P1-15/P1-16 已否决的单年 1950 锚。合成数据测试用 `--allow_missing_temperature_cache` 显式豁免。人口与人均增长率路径由 `01data/ssp/` 生成。**不再生成** MESMER 的三个国家级 climate 文件，改为 `check_national_climate_inputs` 只校验：β/γ 为零、`region_lv`/`region_lv_local` 恒零、`global_lv` 恒零或缺国家一律报错 |
| `build_mesmer_country_weights.py` | `climate_grid_geometry.csv` + GlobPOP + `world_countries.shp` → `unshared_parameters/climate_grid_weights.csv` + `_manifest.csv` | 185×N 稀疏人口权重矩阵（§4.2.1）。GlobPOP 读取、point-in-polygon 与 `TEMPERATURE_PROXY` 均从 `02regression/01scripts/build_country_temperature.py` 导入，结构上保证与 `lctmp_bkly_pw` 同源。同时写出**坐标指纹侧车**（`open_issues.md` R20）：`land_idx` 是位置下标，只比格点数挡不住"同维度、坐标错位"，`helper_mesmer.jl` 加载 W 前会比对该指纹 |
| `precompute_mvn_draws.py` 等的 `test_*.py` | — | pytest 单测：CSV 结构、权重、国家键一致性、MESMER 聚合、状态中心 |
| `test_theta_draws_stats.py` | 生产 CSV → — | §9 的 MVN/VCV 段：Σ_h 对称正定、Cholesky 重构、抽样均值/协方差、秩-6 因子结构、seed 复现（P1-9） |
| `test_clean_output_integration.py` | 阶段 0 产物 → 临时目录 | clean-output 集成测试（P1-10）：阶段 1 四个脚本重跑到空目录、与生产文件逐字节比对、孤儿输入检测、用重建目录跑通模型并对齐 SCC |

### 5.3 `03Claim/data/`

| 路径 | 内容 |
|------|------|
| `fair_data/` | FaIR v2.0 默认参数（气体循环、热响应、初始条件） |
| `fair_data/rcmip_ssp*_emissions_1750_to_2500{,_rebased}.csv` | RCMIP SSP 原始排放情景。**源数据，模型不读**——是下一行的出处，P1-18 的守卫看着它 |
| `fair_data/rcmip_ssp*_emissions_1750_to_2500_rebased_obs.csv` | **模型实际读取的排放**：2016–2024 换成实测清单、2025+ 按 2024 实测水平重基准（§4.1.1）。生成物，由 `build_observed_emissions.py` 产出 |
| `01data/emissions/{edgar,ceds,gfed}/` | 上述校正的清单源数据：EDGAR-2025（4 个 xlsx）、CEDS v_2025_03_18 全球部门汇总（10 个物种）、GFED5.1 全球汇总表（9 个物种）。见 `DATA_SOURCES.md` |
| `mesmer_data/cmip6-ng/` | MESMER 训练用的 CMIP6 ESM `tas` 存档（当前仅 IPSL 的 historical/ssp119/ssp126） |
| `mesmer_data/export_grid_params.py` | 训练 MESMER 并导出网格级参数（§4.2.1），须在 `mesmer_env` 中运行 |
| `mesmer_data/mesmer_results/` | **存档**：旧链路按 2,388 个次国家 `RegionID` 聚合的 302 个文件，不参与运行 |
| `mesmer_data/legacy_national_aggregation/` | **存档**：P1-3 修复前的国家级 `climate_ltcoef/ltinter/region_lv`，是旧 SCC 基准唯一可回溯的来源（旧 `regions.csv` 已不在仓库，无法重建） |
| `shared_parameters/regions_national.csv` | **国家 CLAIM 模型实际读取的区域维度文件**：185 个 ISO3、国家名与回归国家代码 |
| `unshared_parameters/climate_grid_*.csv` | MESMER 网格级参数、几何与人口权重矩阵（§4.2.1）。其中 `climate_grid_innov_chol.csv` 约 146 MB，按 Tier C 不进仓库，由 `export_grid_params.py` 重建 |
| `bootstrap.csv` | 旧损害函数 12 系数 × 1000 bootstrap |
| `shared_parameters/` | 国家级共享参数：`regions_national.csv`、`population_ssp{1..5}.csv`、`climate_T_global.csv` |
| `unshared_parameters/` | 组件专属参数：`socioeconomic_*`（含 `socioeconomic_ypcgrowth_ssp{1..5}.csv`）、`climate_*`、`damage_*` CSV |

#### 5.3.1 SSP 社会经济情景数据来源（population_ssp\* / socioeconomic_ypcgrowth_ssp\*）

情景文件由 `01data/` 下的三个脚本生成（模型端把它们当作原始输入，不含处理过程）：

- **数据来源**：IIASA SSP Scenario Database **v3.2**（2025 年 5 月发布，[下载页](https://ssp.apps.ece.iiasa.ac.at/documentation/basic-drivers)，文件 `ssp_basic_drivers_release_3.2_full.xlsx`）。人口取 IIASA-WiC POP 2025（KC et al. 2024, IIASA WP-24-003），人均 GDP 取 OECD ENV-Growth 2025（Dellink et al. 2025，`GDP|PPP [per capita]`，USD_2017/yr），SSP1–SSP5 逐国 5 年间隔至 2100。2024 过渡年另取世界银行 WDI（`01data/scioeco2021-2025.csv`，`SP.POP.TOTL` 与 `NY.GDP.PCAP.KD.ZG`）。
- **处理管线**（`build_wb_2024_actuals.py` + `ssp/01_extract_ssp.py` → `ssp/02_build_scenarios.py`）：
  1. 历史段 1950–2023 用 PWT 11.0（`rgdpna/pop`）的水平与逐年增长率；
  1b. **2024 过渡年**用实测增长率乘到 PWT 2023 水平上。只取增长率不取水平值——PWT 与 WDI 的国家口径不同，2023 年人口水平差异最大达 47%（CYP 全岛 vs 政府控制区）、20%（MDA 含不含德左）、17%（CUW）、14%（ALB），直接嫁接水平值会在 2023→2024 间劈出假跳变。WDI 覆盖 182/185 国人口与 178/185 国人均增长；其余 7 国由 `build_wb_2024_actuals.py` 的 `MANUAL_2024` 按各自权威来源补齐 2024 实际 GDP 增速，人均增速取 $(1+g^{GDP})/(1+g^{pop})-1$：

  | 地区 | 2024 实际 GDP 增速 | 来源 | 人口增速来源 |
  |---|---|---|---|
  | TWN | +4.59% | DGBAS 初步统计（2025-02-26） | 内政部户籍人口 −20,222 |
  | SSD | −26.14% | IMF WEO | WDI |
  | SYR | −1.5% | World Bank Syria Economic Monitor (Spring 2024) | WDI |
  | YEM | −1.0% | IMF WEO | WDI |
  | VGB | +3.5% | UNCTAD 国别概况 | WDI |
  | AIA | +5.7% | Anguilla Statistics Department / ECCB | PWT 2022→2023 |
  | MSR | +3.56% | Montserrat Statistics Department | PWT 2022→2023 |

  最终 185/185 国的 2024 人均增速均来自实测，仅 AIA/MSR（合计约 2 万人）的人口增速仍用 PWT 2022→2023 回退。2024 是观测年，五个 SSP 情景取值**完全相同**，情景分叉从 2025 年才开始；
  2. 2025–2100 只取 IIASA 的**增长率**（5 年点位间对数线性插值 = 区间内恒定年增长率），从 2025 年起乘到 2024 年水平上（比率拼接，规避 2017-PPP 口径差）；
  3. 2101–2500：人均增长率按 DICE 式指数衰减 g(y)=g(2100)·exp(−0.005·(y−2100))（Barrage & Nordhaus 2023 的 dela 口径），人口 2100 年后恒定；
  4. OECD/WiC 未覆盖的 9 个微型经济体（AIA/BMU/CYM/DMA/KNA/MSR/SXM/TCA/VGB）：人均增长率回退 World、人口恒定；AFG/PSE/SYR/VEN 使用 xlsx 内 "Turbulent Economy Data" 表（IIASA 标注可靠性有限）。

### 5.4 `03Claim/test/`（Julia 测试）

`test_helper_mesmer.jl`（MESMER 网格加载、变率模拟与人口加权聚合，含 AR(1) 平稳方差、lag-1 自相关、空间相关传递等统计性检验）、`test_claim_damage.jl`（组件可编译）、`test_MimiCLAIM.jl` / `test_construct_claim_model.jl`（模型可装配）、`test_helper_damage.jl`（工具函数）、`test_claim_damage_run.jl` / `test_claim_full_run.jl`（合成数据端到端 smoke test）、`test_generated_caches.jl`（真实缓存文件校验）、`test_helper_scc.jl`（SCC 端到端）、`test_ssp_scenarios.jl`（SSP 情景分叉 / ypc_new 健康 / 双折现口径）、**`test_coupling_validation.jl`**（§9 的数值交叉验证：基底重算、权重索引、θ 投影恒等式、闭式 OLS vs NLS、资本记忆脉冲实验、SCC 局部静默、SCC 数量级）。统一入口 `run_all_tests.jl`。

§9 验证清单已于 2026-08-04 基本落地（`open_issues.md` P1-9）：需要数值交叉验证的 13 项全部实现，分布在 `test_coupling_validation.jl`（8 个 testset）与 `03Claim/scripts/test_theta_draws_stats.py`（9 项，MVN/VCV 段）。仍缺的三项——MC 收敛、蒙古漂移、`exp` vs `(1+z)` 敏感性——属分析而非单元测试，随论文稳健性一并做。

> **测试的参照实现刻意不调用生产代码**。基底重算、权重重算、θ 投影都按 `coupling_plan.md` 的公式
> 在测试里独立写一遍；复用生产函数会让“交叉验证”退化为自我复读，两边同时错就同时通过。
>
> `test_coupling_validation.jl` 支持 `CLAIM_TEST_UNSHARED` / `CLAIM_TEST_TABLE` 两个环境变量覆写
> 数据根目录，**唯一用途**是把测试指向人为破坏的副本以验证断言真会失败（已用五类注入验证过）。
> 日常运行不要设置这两个变量。

### 5.5 `02regression/`（估计端，模型的上游）

| 文件 | 说明 |
|------|------|
| `01scripts/01_merge_raw.py` | PWT 11.0 (`rgdpna/pop`) + Berkeley Earth grid/global series + `country_precip_popweighted_1950_2024.csv` + 油价/利率/衰退期/VIX → `raw_merged.dta`；默认不覆盖 `berkeley_full.dta` |
| `01scripts/02_lp_regression.do` | **核心估计**：构造 Hamilton 冲击（§3.1）、正交化、中心化（§8 部分定义 barT_W/lny0_W）、跑 §3.2 的 LP（h=0..10）导出 `lp_main.csv`（21 个 VCV 元素）；§10.5 辅助 LP 导出 `temp_irf.csv`（φ^T）；§11 delta 方法插补 `lp_imputed.csv` |
| `01scripts/03_invert_damage.py` | §3.3 的单核加权 NLS 反演（ν=3）→ `damage_params.csv` + `struct_params.csv` + 诊断图 |
| `01scripts/04_plot_slowB.py` | 慢衰减（小 B）国家的外推诊断图 |
| `01scripts/05_plot_boundary.py` | **B 撞界诊断**（`open_issues.md` P1-12）：全球渠道 185/185 撞界的归属、两族退化极限的核形状、代表国拟合（含 ±1 SE 与外推）。产出 `04figure/B_boundary_overview.png` 与 `B_boundary_examples.png` |
| `03table/lp_main.csv` | pooled 6 系数 + 完整 VCV，horizon h0–h10 |
| `03table/lp_imputed.csv` | 185 国 × 11 horizon 的 θ_c 与 SE_θc |
| `03table/damage_params.csv` | 每国每渠道：model(1/2), A, B, C, D/ω, R², R²₅, RMSE, AIC, α, δ, g |
| `03table/temp_irf.csv` | φ_g, φ_l（h=0..10） |
| `02data/struct_params.csv` | 每国 α, δ, g, s, q0, flag_fallback |
| `damage_inversion_plan.md` / `method_code_audit.md` | 反演方法论与计量审计 |

### 5.6 `03Claim/` 根目录文档

| 文件 | 说明 |
|------|------|
| `coupling_plan.md` | B&K 损害耦合技术方案 v3.1（本文档 §3 的实现细化版，含 Stata 修改、预处理 CLI、Mimi 组件伪码、验证清单） |
| `open_issues.md` | **问题清单与修改记录**（P0/P1/P2 + 论文 caveat + 修复顺序）；问题修复后保留原条目，在标题标注“已修改”并记录过程与验证证据 |
| `paper_model_comparison_plan.md` | 论文五方案对比的三开关设计 |
| `CLAIM_model_guide.md` | 本文档 |

---

## 6. SCC 评估的实现路径

把 §3 的理论落到代码，完整路径分五个阶段：

### 阶段 0：计量估计（Stata + Python，一次性）

```powershell
python 02regression/01scripts/01_merge_raw.py
Start-Process -FilePath "D:\24stata\Stata 18.0\StataMP-64.exe" `
    -ArgumentList '/e','do','02regression/01scripts/02_lp_regression.do' -Wait -NoNewWindow
python 02regression/01scripts/03_invert_damage.py
python 02regression/01scripts/04_plot_slowB.py
python 02regression/01scripts/05_plot_boundary.py   # B 撞界诊断，见 open_issues.md P1-12
```

产出：`lp_main.csv`、`lp_imputed.csv`、`temp_irf.csv`、`damage_params.csv`、`struct_params.csv`，以及 `02regression/04figure/` 下的诊断图。

### 阶段 1：预处理（Python，一次性）

```powershell
# ① 排放输入的观测校正（§4.1.1）。仅在换 EDGAR/CEDS/GFED 版本时需要重跑；
#    改完要同步更新 test_fair_emissions_integrity.py 的 EXPECTED_MD5_OBS。
python 03Claim/scripts/build_observed_emissions.py

python 03Claim/scripts/prepare_national_inputs.py `
    --struct_params 02regression/02data/struct_params.csv `
    --analysis_dta 02regression/02data/analysis.dta `
    --berkeley_dta 02regression/02data/berkeley_full.dta `
    --temperature_cache 02regression/02data/country_temperature_berkeley_grid_pw.csv `
    --output 03Claim/data
python 03Claim/scripts/precompute_basis.py `
    --damage_params 02regression/03table/damage_params.csv `
    --struct_params 02regression/02data/struct_params.csv `
    --temp_irf 02regression/03table/temp_irf.csv `
    --lp_imputed 02regression/03table/lp_imputed.csv --output 03Claim/data/unshared_parameters
python 03Claim/scripts/precompute_cholesky.py `
    --lp_main 02regression/03table/lp_main.csv --output 03Claim/data/unshared_parameters
python 03Claim/scripts/precompute_mvn_draws.py `
    --lp_main 02regression/03table/lp_main.csv `
    --n_draws 1000 --seed 42 --output 03Claim/data/unshared_parameters
```

### 阶段 1b：MESMER 国家级气候参数（§4.2.1）

仅在换 ESM、换训练情景或改模拟域时重跑；日常运行直接消费已有产物。

```powershell
# ① 训练并持久化网格参数（须在装有 MESMER 的 mesmer_env 中运行，约 40 分钟）
#    局地化半径由交叉验证在 --L_start:--L_interval:--L_end（默认 750:250:2000）上选出；
#    若选中值落在窗口端点会直接报错，此时先跑 localisation_radius_scan.py 重新定位窗口。
cd 03Claim/data/mesmer_data
conda run -n mesmer_env python export_grid_params.py --output ../unshared_parameters
cd ../../..

# ② 人口权重矩阵 W（185 × 网格点数，约 5 分钟）
python 03Claim/scripts/build_mesmer_country_weights.py `
    --land_grid 03Claim/data/unshared_parameters/climate_grid_geometry.csv `
    --shapefile 03Claim/data/shapefile/world_countries.shp `
    --regions 03Claim/data/shared_parameters/regions_national.csv `
    --output 03Claim/data/unshared_parameters/climate_grid_weights.csv `
    --diagnostics 03Claim/data/unshared_parameters/climate_grid_weights_diagnostics.csv
```

```julia
# ③ 投影到国家层并抽一条带 seed 的变率实现（Julia，约 16 秒）
include("03Claim/src/helper/helper_global.jl")
include("03Claim/src/helper/helper_mesmer.jl")
regions = load_national_regions("03Claim/data/shared_parameters")
regenerate_national_mesmer_inputs("03Claim/data", regions, 1950:2500; seed=42)
```

产出 `climate_ltcoef.csv`、`climate_ltinter.csv`、`climate_region_lv.csv` 与记录 seed 的 `climate_region_lv_manifest.csv`。随后跑一次 `prepare_national_inputs.py` 即可校验三者健康（零值 β/γ 或恒零变率会直接报错）。

### 阶段 2：确定性运行（Julia）

```julia
include("03Claim/src/MimiCLAIM.jl")
m = get_model(:national; climate=:fair, scenario="ssp245",  # 排放+社会经济（SSP2-4.5）
              theta_id=0)                     # 0 = 点估计（不扰动）；默认 1950–2150
run(m)
loss_bau = m[:claim_damage, :loss]              # BAU 全路径损失
sum(m[:claim_damage, :clipped])                 # 钳制诊断，主结果应为 0
```

时间轴：默认 `end_year=2150`；`climate=:fair` 下 `end_year` 可取不晚于 2500 的任何一年（RCMIP 排放情景与 SSP 社会经济情景文件均覆盖到 2500）。

### 阶段 3：SCC 点估计

```julia
scc = compute_scco2(m; year=2025, eta=1.45, prtp=0.015)   # $/tCO₂，默认积分到模型末年
res = compute_scco2(m; year=2025, return_details=true)
res.scc            # 全球 SCC（$/tCO₂）
res.scc_regional   # 分国 SCC
```

交互式演示（SCC + 温升/人均 GDP 折线 + 分国 SCC 地图 → `03Claim/results/`）：打开 `03Claim/src/MimiCLAIM.ipynb`（Julia **1.12** / Cursor 扩展；notebook 内含 `get_model` 入口代码），按 cell 依次运行。写出 CSV 用 `CSV.write` + `DataFrame`，图在 notebook 内联显示并 `savefig` 到 `results/`。默认 `pulse_year=2025`、`last_year=2150`、`plot_from=1950`。

内部机制：先运行 Base；B&K 模型将 `T0_global` 解析为 Base/Pulse 共用参考温度，BHM 模型无需该参数。随后由 `get_marginal_model` 复制模型、注入 `emissionspulse`（脉冲经 FaIR 变成温度增量 → 经 MESMER 与活动损害组件变成损害增量）→ 双轨 `loss` 差分 → Ramsey 折现加总（§3.6）。要求 `climate=:fair` 构建的模型。正式默认：脉冲年 2025、`bau_ref_year=2025`、积分至 2150。

### 阶段 4：嵌套 Monte Carlo（不确定性区间）

```julia
scc = Vector{Float64}(undef, M_d)
for md in 1:M_d                                   # 损害参数不确定性（θ 预抽样）
    m = get_model(:national; climate=:fair, scenario="ssp245", theta_id=md)
    scc[md] = compute_scco2(m; year=2030)
end
quantile(scc, [0.05, 0.5, 0.95])
```

外层气候不确定性尚未实现，且其形态已修订——外层是 FaIR **参数集合**而非缓存的温度路径（§3.5.2、`open_issues.md` P1-7）。逐成员写法：

```julia
m = get_model(:national; climate=:fair, scenario="ssp245", end_year=2150)
update_param!(m, :temperature, :q, [q1, q2, q3])
update_param!(m, :temperature, :decay_factor, exp.(-1.0 ./ [d1, d2, d3]))
update_param!(m, :radiative_forcing, :co2_f, [f1, f2, f3])
update_param!(m, :co2_cycle, :r0_co2, r0)   # rU/rT/rA 同理
scc = compute_scco2(m; year=2030, last_year=2100)
```

生产代码无需改动；覆写须在 `compute_scco2` 之前完成，Base/Pulse 两轨才会共用。

报告 SCC 的 5% / 50% / 95% 分位，并按 `coupling_plan.md` §9 做收敛与对标检验（B&K 原文、EPA 数量级）。

---

## 7. 运行指南

### 7.1 环境要求

| 软件 | 用途 | 备注 |
|------|------|------|
| Stata 18 (`D:\24stata\Stata 18.0\StataMP-64.exe`) | LP 回归 | 需 `reghdfe`、`ftools`、`xtscc`（.do 内自动 `ssc install`）；**必须经 PowerShell 启动**（Git Bash 会把 `/e` 误解析） |
| Python 3.11 (Anaconda) | 数据预处理 | pandas, numpy, scipy, openpyxl, matplotlib |
| Julia + Mimi.jl | IAM 与 SCC | 另需 DelimitedFiles, CSVFiles, DataFrames, Distributions；旧 notebook 可视化还用 PyCall/PyPlot/VegaLite |

### 7.2 常见坑（实战经验）

1. **Stata 批处理**：`/e` 模式是 GUI 进程，`&` 会立即返回——须轮询 `Get-Process StataMP-64` 等其退出再读 log；
2. **Stata 行内注释**只能用 `//`，代码行尾的 `*` 会被解析为乘法；
3. **全球冲击变量**：主规格必须用 Berkeley 陆海合并口径；自建的 CRU 陆地口径给出零全球渠道是**预期结果**（印证 B&K 海洋论点），不要"修复"它；
4. **不要放年份 FE**（吸收全球冲击）；`absorb(iso3#year)` 之类需先 `encode` 成数值；
5. **中心化口径**：`damage_state_centers.csv` 的 barT_W/lny0_W 必须与 Stata §8 的无加权均值逐位一致；前向模型中所有状态变量绝不允许二次中心化；
6. **horizon 索引**：数学的 h=0..10 对应 Julia 数组下标 1..11，`closed_form_ols` 与基底数组务必对齐；
7. **Base/Pulse 共用 T0_global**：改评估年只改脉冲与折现起点。

### 7.3 学习路径建议

1. 读本文档 §1–3 建立理论图景；
2. 跑通 `julia 03Claim/test/test_claim_full_run.jl`，确认国家层级 socioeconomic → climateregional → claim_damage 链条可运行；
3. 读 `coupling_plan.md` §1（数学）与 §4（组件设计），对照 `src/components_damage/claim_damage.jl` 逐行理解七步计算；
4. 读 `02_lp_regression.do`（变量字典在文件头注释）与 `03_invert_damage.py`；
5. 读 `open_issues.md` 了解当前缺陷，按其第 5 节顺序参与修复。

### 7.4 如何修改模型（按任务索引）

**换情景**
- 改哪里：`get_model(:national; scenario="ssp119"/"ssp126"/"ssp245"/"ssp370"/"ssp585")`，同时切换 FaIR 排放路径与对应的社会经济基线（`population_sspX.csv` / `socioeconomic_ypcgrowth_sspX.csv`）。
- 需要重跑：FaIR → MESMER → damage → SCC 全链（即重新 `get_model` + `compute_scco2`）。更新情景数据：改 `01data/ssp/` 下的 `01_extract_ssp.py` / `02_build_scenarios.py`；换 2024 过渡年的观测数据源改 `01data/build_wb_2024_actuals.py`（§5.3.1）。

**换 ESM / 换训练情景 / 换一条天气变率实现**
- 换 ESM 或训练情景：`conda run -n mesmer_env python export_grid_params.py --esm ... --scenarios h-ssp126 h-ssp245 ...`。MESMER 的模式缩放系数按设计是**情景无关**的，故只存一套、不按情景分文件；当前只在 `h-ssp126` 上训练是因为仓库仅有 IPSL 的 historical/ssp119/ssp126，补齐高强迫情景数据后改这一个参数即可重训。换 ESM 或改模拟域后必须重建人口权重矩阵（网格点数会变）。
- 换一条天气变率实现：改 `regenerate_national_mesmer_inputs(...; seed=...)` 的 seed 重新生成 `climate_region_lv.csv`。**SCC 的 Base/Pulse 两轨必须共用同一条实现**（§3.6 局地渠道相消的前提），因此只能在 `get_model` 之前生成一次，不可在两轨之间重抽。
- 不要恢复 2,388 `RegionID` 中间层或面积加权聚合（`open_issues.md` P1-3）。`prepare_national_inputs.py` 现在只校验这三个国家级文件、不生成它们，遇到零值 β/γ 或恒零变率会直接报错。

**修改 B&K 损害函数**
- 只改核形状族（M1/M2 之外的新核）：改 `03_invert_damage.py`（反演）+ `precompute_basis.py`（基底）+ `helper/helper_damage.jl::closed_form_ols`（若基底数目/线性结构变化）。
- 只改 θ 状态投影（如加二次项）：同步修改 `02_lp_regression.do`（交互项）、`theta_parameters`（参数展开）与 `src/components_damage/claim_damage.jl` 的 Step 2 投影公式。
- 改 Solow 结构参数来源：改 `03_invert_damage.py::load_struct_params`。

**修改折现参数**
- 改哪里：`compute_scco2(...; eta=, prtp=)` 调用实参；当前 API 接受所有国家共用的 `Float64` 标量。
- 注意：折现只影响 SCC，不影响 BAU 损失路径；当前未实现逐国 `eta`/`prtp` 或 equity weighting。

**增加新的不确定性维度**
- B&K 管线：损害参数不确定性走 `damage_theta_base.csv` / `damage_theta_draws.csv`（改 `precompute_mvn_draws.py`）；气候不确定性走 FaIR 参数集合逐成员覆写（§3.5.2；plan §5 的路径缓存方案已作废，见 `open_issues.md` P1-7）。

---

## 8. 当前状态与已知问题

**管线成熟度**（截至 2026-07-30）：

| 环节 | 状态 |
|------|------|
| LP 估计（θ 与 6×6 VCV） | ✅ 完成 |
| 温度 IRF φ^T | ✅ 完成（2026-07-28 重估：全样本 `berk_anom` + 温度自身滞后，φ_g[0] 精确为 1，路径单调为正） |
| 单核反演（ν=3）+ 结构参数 | ✅ 完成（全球渠道 $z_{rms}$ 中位 0.810，理论值 0.798；B 仍有 57/185 撞下界，见 `open_issues.md` P1-12） |
| Python 预处理三件套 | ✅ 完成 |
| claim_damage 组件 + 闭式 OLS | ✅ 完成（$A$ 与 $C$ 由 2×2 加权正规方程联合求解，2026-08-01 起不再锁定 $A=θ₀/Y1[0]$，见 §3.3.3.2；损害起始年参数化已修复） |
| 国家级输入数据 | ✅ 已重新生成（`climate_tempbase`/ypc0 为绝对量、ypcgrowth 用 expm1） |
| MESMER 国家参数 | ✅ 完成（2026-07-30 重构：取消 2,388 `RegionID` 中间层，网格人口加权直投 185 国；模拟域 5,121→5,940 点救回岛国；零值国家 6→0、恒零变率 26→0；§4.2.1） |
| FaIR 接入损害模型 | ✅ 完成（`climate=:fair` 默认模式，1950 spin-up，端到端可跑） |
| SCC 模块 | ✅ 完成（`helper_scc.jl::compute_scco2`，支持活动的 B&K/BHM 损害组件并按 `ypc_new` 增长率折现；数值基准以当前 SCC 测试输出为准） |
| SSP 社会经济情景 | ✅ 完成（IIASA SSP v3.2，经统一 `scenario` 参数接入，§5.3.1；2024 过渡年已换为实测增长率） |
| 嵌套 MC | ⚠️ θ 维度可用（外部循环 theta_id）；外层气候维度未做，已改为 FaIR 参数集合方案（P1-7，B2 构建中） |
| 验证测试（plan §9 的 20 项） | ❌ 仅 smoke test，数值交叉验证项全缺 |

**未解决问题**：完整清单与证据见 `open_issues.md`，当前**无 P0 级问题**。

修复时间线与对 SCC 的影响：

| 日期 | 问题 | SCC(2025, ssp245, →2150) |
|---|---|---|
| 2026-07-28 | P0-A：`gtmp_bkly_aw` 绝对温度/距平混装使 φ^G 完全失真（论证见 `open_issues.md` §5.1–5.3） | 1007.37 → 600.76 |
| 2026-07-29 | P1-11：清理 P0-A 根因，生产面板不再携带 `gtmp_bkly_aw` | 600.76（不变） |
| 2026-07-30 | P1-3：MESMER 改人口加权直投（§4.2.1） | 600.76 → 596.32（−0.7%） |
| 2026-07-30 | P1-4：`T_state` 改为 30 年滚动气候态（§3.4.1） | 见下表 |
| 2026-07-30 | P1-16：历史段改用观测国家温度（§4.2） | 见下表 |

同口径（ssp245，→2100）的连续影响：

| 顺序 | 问题 | SCC | 变化 |
|---|---|---|---|
| 起点 | 旧 MESMER 面积加权链路 | 515.32 | — |
| ① | P1-3：MESMER 人口加权直投 | 510.43 | −0.9% |
| ② | P1-4：`T_state` 改 30 年滚动均值 | 494.57 | −3.1% |
| ③ | P1-16：历史段改用观测 + 前向段观测锚 | **502.03** | +1.5% |

三项净效果 −2.6%。逐项机制：

- **P1-3 位移小**，因为全球损害渠道直接用 $T^{global}-T^{global}_{ref}$ 而非 `regtmp`，β/γ 只通过状态依赖的 θ 投影间接影响 SCC，且局地渠道在 Base/Pulse 差分中相消（§3.6）；跨国均值 β 几乎未动（1.0922→1.0897）。
- **P1-4 影响最大**，因为 `T_state` 直接驱动 θ：30 年后视均值比"滞后一年"更冷（ssp245 下约 0.3 °C），$\theta_1<0$ 故 θ 变得不那么负、损害减小。
- **P1-16 方向相反**：历史段改用观测后 2025 年的气候态由 19.822 修正到真实的 20.037（+0.215，最差国家 +1.47），θ 更负、损害增大。

**P1-14**（MESMER 局地化半径撞搜索下界）已于 2026-07-31 修复：交叉验证窗口从撞界的 `1500:250:10000` 改为 `750:250:2000`，选出内点解 $L_{loc}=1250$ km，邻国天气相关下降 0.03–0.09（FRA-DEU 0.708→0.655）而远距离对不动，SCC 502.03→501.95（−0.016%，局地渠道在边际中相消故点估计几乎不变）。

**P1-12**（损害核 $B$ 撞界）已于同日结案，结论是**边界本身无需修改**，详见 §3.3.1；同时把模型择优判据由 $R^2_5$ 改为加权全段 $R^2_w$，SCC 501.95→484.57（−3.5%）。

**P1-17**（$A$ 解析锁定）已于 2026-08-01 修复：$A$ 改为与形状参数联合估计，见 §3.3.3.2，SCC 484.57→484.04（−0.1%）。

**P1-6**（闭式 OLS 权重口径）已于同日结案，结论是**保留逆方差权重 $1/SE_h^2$**，见 §3.3.4.1。当前优先项转为局地变率均值不为零（P1-13）。

---

## 9. 参考文献

1. Bilal, A. & Känzig, D. R. (2026). *The Macroeconomic Impact of Climate Change: Global vs. Local Temperature*. Quarterly Journal of Economics.（复现数据：Harvard Dataverse doi:10.7910/DVN/O4AUDQ）
2. Hamilton, J. D. (2018). Why You Should Never Use the Hodrick-Prescott Filter. *Review of Economics and Statistics*, 100(5), 831–843.
3. Nath, I., Ramey, V. & Klenow, P. (working paper). *How Much Will Global Warming Cool Global Growth?*
4. Leach, N. J. et al. (2021). FaIRv2.0.0: a generalized impulse response model for climate uncertainty and future scenario exploration. *Geoscientific Model Development*, 14, 3007–3036.
5. Beusch, L., Gudmundsson, L. & Seneviratne, S. I. (2020). Emulating Earth system model temperatures with MESMER. *Earth System Dynamics*, 11, 139–159.
6. Feenstra, R. C., Inklaar, R. & Timmer, M. P. Penn World Table 11.0.
7. IIASA SSP Scenario Database v3.2 (May 2025). 人口：KC, S. et al. (2024). *Updating the Shared Socioeconomic Pathways (SSPs) Global Population and Human Capital Projections*. IIASA WP-24-003；GDP：Dellink, R. et al. (2025). OECD ENV-Growth 2025 projections. https://ssp.apps.ece.iiasa.ac.at/
8. Barrage, L. & Nordhaus, W. (2023). *Policies, Projections, and the Social Cost of Carbon: Results from the DICE-2023 Model*. NBER WP 31112.（2100 后增长率外推的 dela 衰减口径）
8b. 2024 过渡年数据源：World Bank, *World Development Indicators*（导出于 2026-07-13，系列 `SP.POP.TOTL`、`NY.GDP.PCAP.KD.ZG`）；WDI 未覆盖的 7 个地区另取 DGBAS 初步统计（TWN）、内政部户籍人口（TWN 人口）、IMF WEO（SSD、YEM）、World Bank *Syria Economic Monitor* Spring 2024（SYR）、UNCTAD 国别概况（VGB）、Anguilla Statistics Department / ECCB（AIA）、Montserrat Statistics Department（MSR）。逐条出处与检索日期见 `01data/build_wb_2024_actuals.py` 的 `MANUAL_2024`。
9. 本项目文档：《The Macroeconomic Damage Module: A Non-Stationary Capital Dynamics Framework》（`references/`）；`coupling_plan.md`；`02regression/damage_inversion_plan.md`；`open_issues.md`。




