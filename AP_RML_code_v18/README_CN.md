# AP_RML 1D 精简版

本版本只保留当前 1D 数值实验需要的入口：

- `run/run_wkb_main_1d.m`：WKB 双 theta 网格求解；
- `run/run_direct_main_1d.m`：direct density 求解；
- `post/post_main_1d.m`：统一后处理，自动识别 WKB/direct 数据。

## 关键整理

- 删除 `cfg.profile='baseline'` 及所有 baseline/profile 参数预设。
- 删除 `+model/default_params_1d.m`、`+model/default_params_2d.m`、`+model/apply_patch_preset_1d.m`。
- 删除旧示例脚本、旧 post 脚本、2D legacy 目录和 doc 目录。
- 模型参数、初值、网格、时间步和算法选项都写在两个主程序开头。

## 使用

WKB：

```matlab
run/run_wkb_main_1d.m
```

Direct：

```matlab
run/run_direct_main_1d.m
```

后处理：

```matlab
post/post_main_1d.m
```
