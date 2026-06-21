本整理版只保留三个用户入口：

1. run/run_wkb_main_1d.m
   一维 WKB 双 theta 网格求解。
   所有模型参数、初值、网格、dt、T、算法选项和保存选项都在文件开头显式设置。
   不再使用 cfg.profile，也不再依赖 +model/default_params_1d.m。

2. run/run_direct_main_1d.m
   一维 direct density 求解。
   参数同样全部在文件开头显式设置。

3. post/post_main_1d.m
   统一后处理。自动扫描 data/main_1d 下的 WKB/direct 最终结果，
   自动识别 n 数据和 WKB 数据并构造 P(theta)。也可切换为 manual 模式手动给标签。

内部入口：
   run/run_case_1d.m
   不建议平时修改。它只负责把主程序中显式给出的 cfg 传给求解器。

已删除：
   cfg.profile / baseline 参数预设；
   +model/default_params_1d.m, +model/default_params_2d.m, +model/apply_patch_preset_1d.m；
   旧 run/archive_examples、旧 post 脚本、src_2d_legacy、doc 目录。
