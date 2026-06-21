# Bugfix 2 notes

本次修复主要针对 MATLAB 路径和运行时稳定性：

1. 修复 `utils.ensure_dir('data/...')` 在部分 MATLAB 环境中被解释为 URL 的问题。输出目录统一使用项目根目录下的绝对路径 `data/...`。
2. `ensure_dir` 增加 MATLAB `mkdir`、Java `mkdirs` 和系统 `mkdir -p` fallback。
3. 删除 `eigs` 的 `opts.issym` 和 `opts.isreal`，避免矩阵输入时被 MATLAB 忽略并刷 warning。
4. 1D 特征值问题改为 trapezoidal-weighted symmetric similarity transform，保留 reflected Neumann 离散的特征值结构。
5. 增加 `run_smoke_test_1d.m`，可先用小网格检查路径、求解器和保存功能。
6. 1D 默认不再把 `src_2d_legacy` 加入路径；只有运行 2D wrapper 时用 `startup_AP_RML(true)` 加入 legacy 路径。

建议先运行：

```matlab
startup_AP_RML
run_smoke_test_1d
```
