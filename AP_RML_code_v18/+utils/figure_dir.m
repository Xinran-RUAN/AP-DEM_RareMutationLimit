function outdir = figure_dir(dataDir, varargin)
%FIGURE_DIR Canonical figure directory under an experiment data directory.
%  Default: dataDir/figures. Extra arguments create subfolders under figures.
if nargin < 1 || isempty(dataDir)
    dataDir = utils.output_dir();
end
outdir = utils.ensure_dir(fullfile(char(dataDir), 'figures', varargin{:}));
end
