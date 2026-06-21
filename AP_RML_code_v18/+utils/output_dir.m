function pathname = output_dir(varargin)
%OUTPUT_DIR Build an absolute path under the project data/ directory.
pathname = fullfile(utils.project_root(), 'data', varargin{:});
end
