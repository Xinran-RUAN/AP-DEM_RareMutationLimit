function filename = compare_collection_file(outdir)
%COMPARE_COLLECTION_FILE Canonical paired comparison collection filename.
filename = fullfile(char(outdir), 'snapshots_compare.mat');
end
