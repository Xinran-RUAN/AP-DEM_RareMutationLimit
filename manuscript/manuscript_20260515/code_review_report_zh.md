# 代码核对报告（WKB rare-mutation dispersal model）

## 结论

当前 1D/2D 代码的核心思路与文档中的 WKB split 模型一致：使用
`n = W exp(u/eps)`，先由 `rho` 求离散主特征值 `H(theta)`，再更新相位 `u` 与振幅 `W`。但是，旧代码与当前稿件之间存在若干重要不一致，尤其是实验参数、全离散振幅方程、归一化方式、密度重构、直接求解器路径以及 2D 主程序。新的 `+src` 代码已按稿件模型重构，并保留 split WKB 与 density-compatible WKB 两个振幅更新。

## 主要不一致与风险点

1. **实验模型参数不一致**
   - 文档 Section 4 的 baseline 是
     `K(x)=1+0.5 cos(2*pi*x)`，
     `D(theta)=0.2+0.4(1-cos(2*pi*(theta-theta_m)))`，`theta_m=0.35`。
   - 旧 1D/2D 主程序仍使用 `K=1+20*(1-4*(x-0.5)^2)^8` 或其 2D 类似形式，以及指数型 `D`。
   - 旧 `initial_D(theta_m,theta)` 在 `theta_m ~= 0.5` 时并未真正把极小值放在输入的 `theta_m`。
   - 已在 `+model/default_params_1d.m` 中统一为稿件 baseline 和 heterogeneous 两套 profile。

2. **全离散 split WKB 振幅更新与文档公式不完全一致**
   - 文档的 fully discrete split 公式把 `D_theta F(W^{n+1},u^{n+1})` 和 `-2 eps delta_theta^2 u^{n+1}` 对应的线性项放在新时层。
   - 旧 `solve_w.m` / `solve_w_fem.m` 中 upwind transport 和 `2 eps W delta_theta^2 u` 是显式放到右端的。
   - 这不是半离散模型错误，但它是一个 IMEX 变体；若不写入文档，稳定性/正性结论不能直接对应。
   - 新代码默认使用 implicit upwind transport；也保留 `par.transportImplicit=false` 作为旧 IMEX 风格。

3. **文档的 density-compatible variant 旧代码未实现**
   - 稿件中 density-compatible 更新要求离散层面严格满足
     `eps n_t - D delta_xx n - eps^2 delta_thetatheta n = n(K-rho)`。
   - 旧代码没有该变体。
   - 新增 `+src/update_w_density_compatible_fd_1d.m`，通过比值
     `exp((u_{k±1}-u_k)/eps)` 构造 `E_k^{-1} delta_thetatheta(W E)_k`。

4. **相位 Hamiltonian 的实现高阶化但文档未说明**
   - 旧代码用 `WENO5_left_RXR/right_RXR` 做单边导数重构，再用 Roe/entropy-fix flux 近似 `-|u_theta|^2`。
   - 文档当前主要写 Godunov/Lax-Friedrichs monotone Hamiltonian。WENO5 不满足理论中使用的单调性假设，因此只能作为数值高阶选项。
   - 新代码默认 `phaseHamiltonian='godunov'`，并保留 `weno5` 选项；TeX 中补充了 implementation appendix。

5. **`rho` 重构方式与文档需要进一步区分**
   - 旧代码后期默认改用 Laplace 重构 `solve_rho_laplace`，但文档中的 Assumption 3.2 对 direct composite quadrature 最直接。
   - 旧 Laplace 重构在峰值靠近周期边界时使用 `theta(km)=1-dtheta` 会导致局部三点插值坐标错误；当 `u'' >= 0` 或曲率很小时也没有 fallback。
   - 新代码默认使用 log-stabilized direct quadrature；Laplace/hybrid reconstruction 作为可选诊断。周期边界采用局部坐标 `{-dtheta,0,dtheta}` 并带 fallback。

6. **归一化 gauge 不一致**
   - 文档写的是 `max_k u_k=0` 的 gauge，并相应缩放 `W` 以保持 `n` 不变。
   - 旧 `normalize_u.m` 采用 `int exp(u/eps) dtheta = 1` 的积分归一化，虽然也保持 `n` 不变，但与文档表述不一致，且使用固定 `0.001` spline 积分。
   - 新代码使用 `+src/normalize_gauge.m` 的 max-gauge。

7. **直接求解器路径/命名空间错误**
   - 旧 `run/original_bio.m` 未调用 `startup`，且调用 `initial_D`、`prepare_part_ori`、`solve_original_rho` 等未加 package 前缀，按当前目录结构会失败。
   - 新代码统一为 `+src/run_direct_1d.m`。

8. **停止准则不完整**
   - 旧 1D WKB 主程序只用 `norm(u-u_new)` 判断收敛，可能在 `W` 未收敛时停止。
   - 新代码使用 `u` 与 `W` 的联合相对残差。

9. **1D Neumann 边界离散需要在文档中说明**
   - 旧代码在包含端点的网格上使用反射 ghost（边界行系数为 `[-2,2]/dx^2`），而稿件半离散部分写成了较抽象的 ghost values。
   - 新 TeX 中建议增加 implementation note：端点网格采用反射 Neumann；理论部分仍可保留抽象记号。

10. **2D FEM 部分仍建议作为 legacy/extension 使用**
    - 2D 代码主体可运行思路与 WKB 模型一致，但还混有旧参数、旧归一化和显式 transport。
    - 新结构中将其放在 `src_2d_legacy`，提供 `main_bio_2d_restructured.m` wrapper，用新参数和 max gauge 做二维实验入口。
    - 若投稿数值重点依赖 2D，建议下一步把 FEM 质量矩阵、刚度矩阵、transport 矩阵也完全 package 化。

## 新结构中对应的关键文件

- WKB 主求解器：`+src/run_wkb_1d.m`
- direct density 求解器：`+src/run_direct_1d.m`
- phase update：`+src/step_phase_1d.m`
- split amplitude update：`+src/update_w_split_fd_1d.m`
- density-compatible amplitude update：`+src/update_w_density_compatible_fd_1d.m`
- log-stable/Laplace density reconstruction：`+src/reconstruct_rho_1d.m`，`+src/reconstruct_laplace_1d.m`
- remainder diagnostic：`+src/compute_remainder_gamma_1d.m`
- 实验脚本：`run/run_test1_*.m` 到 `run/run_test6_*.m`
