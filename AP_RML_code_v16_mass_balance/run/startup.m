function root = startup()
%STARTUP 将工程根目录、run/ 和 post/ 加入 MATLAB 路径。
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_AP_RML(false);
end
