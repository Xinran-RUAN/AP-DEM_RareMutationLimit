function root = startup_AP_RML(includeLegacy2D)
%STARTUP_AP_RML Add the project folders needed by the AP-RML code.
%   startup_AP_RML() adds the package root, run/, and post/.
%   startup_AP_RML(true) also adds src_2d_legacy/ recursively.
if nargin < 1 || isempty(includeLegacy2D)
    includeLegacy2D = false;
end
root = fileparts(mfilename('fullpath'));
addpath(root);
addpath(fullfile(root, 'run'));
addpath(fullfile(root, 'post'));
if includeLegacy2D
    legacy2d = fullfile(root, 'src_2d_legacy');
    if exist(legacy2d, 'dir') == 7
        addpath(genpath(legacy2d));
    end
end
end
