# 修改说明（2026-05 draft）

## 总体处理思路

这版没有继续把 WKB remainder 条件写成“自然成立”的弱假设，而是把它降调为一个显式的 conditional one-sided stability assumption。同时，为了降低审稿风险，文中新增了一个 density-compatible WKB variant：该变体使重构密度在离散层面严格满足 conservative density balance，因此 weighted remainder 可以由周期离散分部求和直接控制。这样全文有两条线：

1. generic split WKB scheme：stationary AP 结构是 conditional 的；
2. density-compatible WKB variant：one-sided remainder bound 可验证，因此可以作为理论和数值实验中的稳妥备选。

## 主要修改

1. **Introduction / Abstract**
   - 删除了 `\cite{}` 空引用和红色占位文字。
   - 修正了基本假设中 `K(\theta)` 的写法，统一为 `K(\mathbf{x})`。
   - 摘要和引言改为更诚实的表述：强调 fixed-grid stationary AP structure 是在 explicit one-sided stability condition 下得到的，并说明 density-compatible variant 可以验证该条件。

2. **Section 3：半离散格式与 stationary AP 结构**
   - 保留原来的 split WKB semi-discrete scheme。
   - 在 reconstructed density equation 后新增 **density-compatible variant**，对应新的方程 `\eqref{eq:density_compatible_W}` 和 `\eqref{eq:density_compatible_n}`。
   - 将原来的 remainder 条件改写为正式的 Assumption：
     `One-sided weighted remainder stability`，并使用正部形式
     `((\mathcal R_j^{\varepsilon,D})_+ \le C_R m_j^\varepsilon+r_h)`。
   - 将 density reconstruction 的上下界整理为 `Reconstruction stability` assumption。
   - 精简并重写了 conditional density bound 的证明，去掉了过密的 remark，把解释性内容并入正文逻辑。
   - 保留 fixed-grid phase estimate 和 formal AP trait-selection limit，但明确其条件性与 fixed-grid 性质，避免过度宣称完整 AP convergence theorem。

3. **Appendix**
   - 删除原稿中重复的 appendix 和重复的 fully discrete realization。
   - 新增 Appendix `A sufficient condition for the one-sided remainder bound`：
     - 对 density-compatible variant，证明 weighted remainder bound 可由 discrete summation by parts 得到；
     - 对 generic split WKB scheme，给出 bounded exponential ratio 和 local defect bound 形式的充分条件；
     - 明确指出普通的 `u` 离散导数有界不足以推出 epsilon-uniform bound。
   - 整理了 monotone Hamiltonian、conservative flux、fully discrete realization、direct density solver 的 appendix。
   - 添加了简短 bibliography，避免引用全部 unresolved。

4. **Section 4：数值实验部分**
   - 将原来的实验计划扩写为可直接执行的实验设置：
     - baseline 参数：`K(x)=1+0.5 cos(2πx)`，`D(θ)=0.2+0.4(1-cos(2π(θ-θ_m)))`，`θ_m=0.35`；
     - 初值：phase peak 初始放在 `θ_0=0.70`，避免 well-prepared；
     - 明确了 WKB/direct solver 的对比指标、停止准则、log-stabilized density reconstruction、gauge normalization。
   - 设计了 5 组实验：
     1. direct density solver vs WKB solver；
     2. fixed coarse trait grid with decreasing epsilon；
     3. trait-mesh dependence；
     4. nontrivial heterogeneous landscapes；
     5. remainder and density-compatibility diagnostics。
   - 第 5 组实验专门用于回应审稿人对 remainder assumption 的疑问，建议保留。
   - 表格中数值暂用 `--` 占位，跑完数值实验后直接填入即可。

5. **LaTeX 清理**
   - 关闭了 draft label boxes / hyperlink red boxes，改为 clean output。
   - 删除了未使用的 `psfig` 输入。
   - 清理了重复 appendix。
   - 已经用 `pdflatex` 编译两遍；当前版本无 undefined references、无 overfull hbox 警告。

## 建议后续数值策略

- 如果 split WKB 的 remainder diagnostic `\Gamma_R^\varepsilon` 随 `\varepsilon` 明显爆炸，正文数值部分应主推 density-compatible variant 的 stationary results。
- 如果 `\Gamma_R^\varepsilon` 在表 5 中保持温和，则可以保留 split WKB 作为主格式，并把 density-compatible variant 作为理论保障和对照。
- 投稿前建议将表格中的 `--` 全部替换为数据，并加入至少三类图：
  1. trait marginal `P_k` 对比；
  2. phase `u_k` 与选中特征位置；
  3. `\Gamma_R^\varepsilon` 随 `\varepsilon` 的变化。
