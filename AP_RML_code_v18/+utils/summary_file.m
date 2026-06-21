function filename = summary_file(outdir)
%SUMMARY_FILE Canonical experiment summary filename.
filename = fullfile(char(outdir), 'summary.mat');
end
