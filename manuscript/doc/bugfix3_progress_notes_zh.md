# bugfix3：进度输出、中文注释和普通主程序

本版主要响应以下需求：运行主程序时需要看到进度、预计还剩多少；代码中增加中文注释；在 `run/` 下增加非数值测试的普通主程序。

## 新增文件

- `run/run_main_1d.m`：普通单组参数主程序，只跑一组 WKB 参数；可选打开 direct solver。
- `+utils/progress_init.m`、`+utils/progress_update.m`、`+utils/progress_finish.m`：统一进度输出工具。
- `+utils/format_seconds.m`：把 elapsed / ETA 格式化为易读字符串。

## 进度输出内容

WKB/direct/2D legacy wrapper 现在会打印：

- 开始参数：`eps, Nx, Ntheta, T, dt, tol`
- 预计最多步数
- 当前 step / 总 step
- 当前 t / T
- residual
- elapsed 和 ETA
- 已知的 trait 诊断量，例如 `thetaWKB`, `thetaDIR`, `gammaR`, `rhoMax`
- 最终保存路径

输出频率由以下参数控制：

```matlab
par.verbose = true;
par.progressEveryPercent = 5;
par.progressEverySeconds = 10;
par.progressEverySteps = [];
par.historyEveryTime = 0.2;
```

## 中文注释

已在以下核心文件补充中文说明：

- `+model/default_params_1d.m`
- `+model/default_params_2d.m`
- `+src/run_wkb_1d.m`
- `+src/run_direct_1d.m`
- `+src/solve_H_fd_1d.m`
- `+src/reconstruct_rho_1d.m`
- `+src/update_w_split_fd_1d.m`
- `+src/update_w_density_compatible_fd_1d.m`
- `+src/update_direct_density_fd_1d.m`
- `run/run_main_1d.m`
- `run/run_test*.m`

## 建议启动方式

```matlab
startup_AP_RML
run_smoke_test_1d
run_main_1d
```
