function pathname = ensure_dir(pathname)
%ENSURE_DIR Create a local directory if needed and return its absolute path.
% Relative paths are resolved relative to the project root. This prevents
% MATLAB from interpreting paths such as 'data/test3_compare_n' as URLs in
% some environments and makes outputs independent of the current folder.
if nargin < 1 || isempty(pathname)
    pathname = utils.project_root();
else
    pathname = utils.resolve_path(pathname);
end

if exist(pathname, 'dir') == 7
    return;
end

ok = false;
msg = '';
msgid = '';
mkdirException = [];
try
    [ok, msg, msgid] = mkdir(pathname);
catch ME
    mkdirException = ME;
end

if (~ok || exist(pathname, 'dir') ~= 7) && usejava('jvm')
    try
        jf = javaObject('java.io.File', pathname);
        ok = jf.exists() || jf.mkdirs();
    catch
        % Keep the MATLAB mkdir error below.
    end
end

if (~ok || exist(pathname, 'dir') ~= 7)
    q = strrep(pathname, '"', '\"');
    if ispc
        cmd = sprintf('mkdir "%s"', q);
    else
        cmd = sprintf('mkdir -p "%s"', q);
    end
    [status, ~] = system(cmd);
    ok = (status == 0) && (exist(pathname, 'dir') == 7);
end

if ~ok && exist(pathname, 'dir') ~= 7
    if ~isempty(mkdirException)
        error('utils.ensure_dir:MkdirFailed', ...
              'Could not create directory "%s". MATLAB mkdir error: %s', ...
              pathname, mkdirException.message);
    else
        error('utils.ensure_dir:MkdirFailed', ...
              'Could not create directory "%s". mkdir returned msgid="%s", msg="%s".', ...
              pathname, msgid, msg);
    end
end
end
