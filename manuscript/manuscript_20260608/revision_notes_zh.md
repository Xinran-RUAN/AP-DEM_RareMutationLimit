# 修改说明：重构版论文更新

本次修改基于当前代码 `AP_RML_code_v17` 与原稿 `manuscript_20260518_手动修改`。核心目标是把代码中已经形成的重构思想写入论文，并把“全离散格式”部分改写为“全离散重构流程”。

## 1. 对用户总结的核对

你的总结总体正确，但论文中需要保留以下限定，避免把后处理/数值重构误写成半离散 AP 定理的一部分：

- `log P` + WENO 峰定位确实用于小 `eps` 下的 `P(theta)` 重构，主要体现在后处理脚本中；其作用是稳定确定 `theta_m` 位置，并避免直接指数化带来的下溢/上溢。
- 质量修正在求解器层面实现，用一个正的标量因子修正 `rho` 与 `W`，保持质量平衡，同时不改变 `u` 的峰位置。
- Laplace 重构主要用于 `rho(x)` 的稳定恢复，特别是小 `eps` 时由相函数局部二次结构给出指数积分近似。
- `uW` 双网格思想已经在代码中实现：`u` 在较细的 trait 网格上演化，`W` 在较粗的 trait 网格上演化；这与高维空间下的效率优势相匹配。不过当前代码中它主要用于 frozen-H / split amplitude 路径，因此论文中表述为数值实现与效率策略，而不是半离散理论证明的一部分。
- 二次拟合用于峰形/峰宽模拟，尤其服务于小 `eps` 下 `P(theta)` 的窄峰重构。

## 2. 论文主要修改

- 摘要和引言加入重构流程：`log P`、WENO 峰定位、质量修正、Laplace `rho` 重构、二次峰宽拟合和 `uW` 双网格。
- 将第 3.4 节从“Fully discrete realization and time-step restrictions”改为“Fully discrete reconstruction pipeline and time-step restrictions”。
- 第 3.4 节新增/重写以下内容：
  - `rho^n = R_rho[W^n,u^n]` 与 `P^n = R_P[W^n,u^n]` 的重构观点；
  - log-sum-exp 密度重构；
  - Laplace 密度重构；
  - frozen-H 相方程更新；
  - split amplitude 更新；
  - 可选 density-compatible integrating-factor 更新；
  - 可选 Patankar 密度层更新；
  - 质量修正公式；
  - `log P(theta)` 重构、WENO/二次峰定位和峰宽拟合；
  - gauge-invariant residual 与自适应时间步；
  - `uW` 双 trait 网格描述。
- 第 4 节数值实验设置中修正/补充与代码一致的内容：例如 `K(x)=K0-K1 cos(2*pi x)`、`W` 初值通常取 1、重构诊断不再只是直接积分。
- 重写 Test 3，使其突出小 `eps` 下的 trait-marginal 重构和双网格效率。
- 扩展 Test 4，说明 `uW` 双网格在空间高维时的效率意义。
- 附录 C 改写为重构附录，集中说明 log-stable `rho`、Laplace `rho`、`logP`、WENO 峰搜索、二次峰宽、质量归一化和双网格重构。
- 修正了原稿中两个未定义引用，替换为现有文献条目。

## 3. 编译与检查

- 已成功编译为 PDF。
- 编译日志中未发现未定义引用或未定义文献。
- 已将 PDF 渲染为页面图像并抽查重点修改页；未发现明显裁切、黑块或文字破损。
- 仍有少量常见的 LaTeX/hyperref/overfull hbox 警告，不影响正文阅读。
