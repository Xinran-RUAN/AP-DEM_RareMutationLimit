function root = project_root()
%PROJECT_ROOT Return the absolute root directory of this code package.
root = fileparts(fileparts(mfilename('fullpath')));
end
