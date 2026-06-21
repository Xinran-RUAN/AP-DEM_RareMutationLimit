function pathname = resolve_path(pathname)
%RESOLVE_PATH Convert a local path to a usable absolute path.
% Relative paths are resolved under the project root.
if nargin < 1 || isempty(pathname)
    pathname = utils.project_root();
    return;
end
if isstring(pathname)
    if numel(pathname) ~= 1
        error('utils.resolve_path:InvalidPath', 'Path must be a scalar string or char vector.');
    end
    pathname = char(pathname);
elseif iscell(pathname)
    if numel(pathname) ~= 1
        error('utils.resolve_path:InvalidPath', 'Path cell array must contain one entry.');
    end
    pathname = char(pathname{1});
else
    pathname = char(pathname);
end
pathname = strtrim(pathname);
if numel(pathname) >= 2
    if (pathname(1) == '''' && pathname(end) == '''') || ...
       (pathname(1) == '"'  && pathname(end) == '"')
        pathname = pathname(2:end-1);
    end
end
if isempty(pathname)
    pathname = utils.project_root();
elseif pathname(1) == '~'
    homeDir = char(getenv('HOME'));
    if isempty(homeDir)
        try
            homeDir = char(java.lang.System.getProperty('user.home'));
        catch
            homeDir = '';
        end
    end
    if isempty(homeDir)
        error('utils.resolve_path:NoHome', 'Cannot expand a path beginning with ~ because the home directory is unknown.');
    end
    if numel(pathname) == 1
        pathname = homeDir;
    elseif pathname(2) == filesep || pathname(2) == '/' || pathname(2) == char(92)
        pathname = fullfile(homeDir, pathname(3:end));
    else
        error('utils.resolve_path:UnsupportedTilde', 'Only paths of the form ~ or ~/folder are supported.');
    end
elseif ~utils.is_absolute_path(pathname)
    pathname = fullfile(utils.project_root(), pathname);
end
end
