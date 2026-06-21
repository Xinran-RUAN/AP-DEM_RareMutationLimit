function tf = is_absolute_path(p)
%IS_ABSOLUTE_PATH True for local absolute paths on Unix/macOS/Windows.
if isstring(p)
    p = char(p);
end
if isempty(p)
    tf = false;
    return;
end
p = char(p);
if ispc
    tf = ~isempty(regexp(p, '^[A-Za-z]:[\\/]', 'once')) || ...
         strncmp(p, '\\', 2) || strncmp(p, '//', 2);
else
    tf = p(1) == '/' || p(1) == '~';
end
end
