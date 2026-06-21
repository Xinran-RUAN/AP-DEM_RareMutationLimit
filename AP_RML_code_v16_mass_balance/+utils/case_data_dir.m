function outdir = case_data_dir(caseName, solverName, par)
%CASE_DATA_DIR 统一构造一维算例的数据目录。
%
% 数据目录约定：
%   data/算例名/wkb_eps..._dt..._Nx..._Ntheta..._t...
%   data/算例名/direct_eps..._dt..._Nx..._Ntheta..._t...
%
% 其中 solverName 为 'wkb' 或 'direct'。目录名只包含后处理最常用的必要参数。

if nargin < 3
    error('utils.case_data_dir:NotEnoughInputs', ...
        'Usage: utils.case_data_dir(caseName, solverName, par).');
end

caseName = local_safe_token(caseName, 'case');
solverName = utils.normalize_solver_name(solverName);

dirName = sprintf('%s_eps%s_dt%s_Nx%d_Ntheta%d_t%s', ...
    solverName, ...
    local_num_token(local_get_num(par, 'eps', NaN)), ...
    local_num_token(local_get_num(par, 'dt', NaN)), ...
    round(local_get_num(par, 'Nx', NaN)), ...
    round(local_get_num(par, 'Ntheta', NaN)), ...
    local_num_token(local_get_num(par, 'T', NaN)));

if isfield(par, 'dirExtraTag') && ~isempty(par.dirExtraTag)
    dirName = [dirName '_' local_safe_token(par.dirExtraTag, 'tag')];
end

outdir = fullfile(utils.project_root(), 'data', caseName, dirName);
end

function v = local_get_num(s, name, defaultValue)
v = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)) && isnumeric(s.(name)) && isscalar(s.(name))
    v = double(s.(name));
end
end

function s = local_safe_token(v, defaultValue)
if nargin < 2
    defaultValue = 'token';
end
if isempty(v)
    s = defaultValue;
    return;
end
s = local_to_char(v);
s = strtrim(s);
% 只替换常见路径非法字符和空白，中文算例名会保留。
s = regexprep(s, '[\\/:*?"<>|\s]+', '_');
s = regexprep(s, '^_+|_+$', '');
if isempty(s)
    s = defaultValue;
end
end

function tf = local_isstring(v)
tf = false;
if exist('isstring', 'builtin') || exist('isstring', 'file')
    tf = isstring(v);
end
end

function s = local_to_char(v)
if ischar(v)
    s = v;
elseif local_isstring(v)
    s = char(v);
elseif isnumeric(v) && isscalar(v)
    s = num2str(v);
else
    try
        s = char(v);
    catch
        s = 'token';
    end
end
end

function s = local_num_token(x)
% 使用统一的数值命名规则，避免 2.5e-4 被四舍五入成 3e-4。
s = utils.num_token(x);
end
