% RUN_TEST6_2D_LEGACY
% -------------------------------------------------------------------------
% 旧版二维代码入口。二维部分保持 legacy，不纳入一维 run_case_1d。
% -------------------------------------------------------------------------

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_AP_RML(true);

%% ========================= 用户接口 =========================

par2d = model.default_params_2d();
par2d.outdir = utils.output_dir('test6_2d_legacy');
par2d.verbose = true;

%% ========================= 运行 =========================

utils.ensure_dir(par2d.outdir);
fprintf('输出目录: %s%s', par2d.outdir, newline);

try
    result2d = main_bio_2d_restructured(par2d); %#ok<NASGU>
    summaryFile = utils.summary_file(par2d.outdir);
    save(summaryFile, 'result2d', 'par2d');
    fprintf('TEST 6 完成。summary 保存到: %s%s', summaryFile, newline);
catch ME
    fprintf(2, 'TEST 6 失败：%s%s', ME.message, newline);
    rethrow(ME);
end
