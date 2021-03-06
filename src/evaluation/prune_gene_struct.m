function pruned_genes = prune_gene_struct(genes)

keep_fields = {...
    'name', ...
    'strand', ...
    'chr', ...
    'chr_num', ...
    'start', ...
    'stop', ...
    'transcripts', ...
    'transcript_valid', ...
    'exons', ...
    'cds_exons', ...
    'utr5_exons', ...
    'utr3_exons', ...
    'is_valid', ...
    'transcript_coding'
};

rm_fields = setdiff(fieldnames(genes), keep_fields);
%for f=1:length(rm_fields),
%  fprintf('pruning field %s from gene structure\n', rm_fields{f})
%end
pruned_genes = rmfield(genes, rm_fields);
