# AP_RML_code_v16 重写版使用说明

这个版本把一维算例的运行方式统一成一个主入口：

```matlab
run/run_case_1d.m
```

普通算例脚本只负责在文件开头填写 `cfg`，然后调用：

```matlab
runInfo = run_case_1d(cfg);
```

## 1. 推荐工作流

1. 复制 `run/run_template_1d.m`。
2. 只修改开头“用户接口”部分：算例名、求解器、网格、时间步、初值、算法、保存设置。
3. 运行脚本。
4. 用 `post/post_process_1d.m` 填入一个或多个数据目录做后处理。

## 2. 数据目录格式

一维数据统一保存为：

```text
data/算例名/wkb_eps..._dt..._Nx..._Ntheta..._t...
data/算例名/direct_eps..._dt..._Nx..._Ntheta..._t...
```

每个数据目录中常见文件：

```text
case_info.mat
result_wkb_final.mat 或 result_direct_final.mat
checkpoint_wkb_latest.mat 或 checkpoint_direct_latest.mat
history_wkb_latest.mat 或 history_direct_latest.mat
snapshot_wkb_stepXXXXXXX_tT.mat 或 snapshot_direct_stepXXXXXXX_tT.mat
```

## 3. 批量参数规则

一个主程序最多只批量跑一个变量。写法示例：

```matlab
cfg.eps = [1e-2, 1e-3];
```

或者：

```matlab
cfg.Ntheta_wkb = 16;
cfg.Ntheta_direct = [16, 32, 64];
```

`cfg.solvers = {'wkb','direct'}` 不算批量变量。

## 4. 后处理入口

统一入口是：

```matlab
post/post_process_1d.m
```

在开头填：

```matlab
dataDirs = {
    fullfile(root, 'data', '算例名', 'wkb_eps...'),
    fullfile(root, 'data', '算例名', 'direct_eps...')
};

fileNames = {};        % 为空则自动找 final；也可手动指定 snapshot 文件名
opt.jobs = {'history','snapshot'};
opt.snapshotTime = []; % [] 画 final；例如 1 画 t≈1 的 snapshot
```

底层画图函数是：

```matlab
post/plot_case_1d.m
post/plot_history_1d.m
post/plot_snapshot_1d.m
```

## 5. 入口脚本

`run/` 中的脚本现在都只是薄入口，参数集中在开头：

```text
run_main_1d.m
run_template_1d.m
run_smoke_test_1d.m
run_test1_assumption_diagnostics.m
run_test2_long_time_concentration.m
run_test3_wkb_advantage_refinement.m
run_test4_coarse_theta_AP.m
run_test5_accuracy.m
run_timeap_compare_1d.m
run_all_1d.m
run_test6_2d_legacy.m
```

其中一维脚本都调用 `run_case_1d`；二维 legacy 脚本保持单独入口。
