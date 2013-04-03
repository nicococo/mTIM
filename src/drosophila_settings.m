function CFG = drosophila_settings(CFG)

% This script prepares paths and input/output files specialized
% for drosophila test.

% 8G as defined in general_settings.m is too less!
CFG.grid_memreq = 12;
fprintf('12GB for the grid.\n');   

%CFG.PAR.bmrm.K = 250;             % number of cutting planes used -1 + aggregated plane

% parameter combinations to be used for independent training
% (typically overwritten)
CFG.train_params = { ...
  [0.01],   [0.05],   [0.05],  [1500], 'QP', [1]; ...  % max amount of train data available without val data
  [0.01],   [0.5 ],   [0.5],  [1500], 'QP', [1]; ...  % max amount of train data available without val data
  [0.01],   [5.0 ],   [0.5],  [1500], 'QP', [1]; ...  % max amount of train data available without val data
  [0.001],   [5.0],   [0.5],  [1500], 'QP', [1]; ...  % max amount of train data available without val data
};



%%%%%%%%%%%%%%%%%%%%%%%%%
% Reset features
%%%%%%%%%%%%%%%%%%%%%%%%%


CFG.max_train_chunk_len = 200000;
% number of discrete expression levels
CFG.PAR.num_levels = 5; % something like 3-10;
% number of supporting points for each scoring function
CFG.PAR.num_plif_nodes = 15; %2*CFG.PAR.num_levels; %CFG.PAR.num_levels+2 %2*CFG.PAR.num_levels;
% include mate-pair information for paired reads as a learning features
CFG.PAR.use_pair_feats = 0;
% include data on genomic repeats as a learning features
CFG.PAR.use_repeat_feats = 0;
% use splice feats
CFG.PAR.use_splice_feats = 0;
% enforce the monotonicity for some feature scoring functions
% note that this can make the training problem more difficult
% (i.e. time-consuming) to solve!
CFG.PAR.enf_monot_score_funcs = 0;
% use heuristic training procedure
CFG.PAR.constraint_margin = 5;
% options for enforcing a certain  state given a certain range of feature values
CFG.PAR.switch_features = 1;
% the maximum allowed value for the low-coverage block feature in exon
% and intron states (forcing to decode these as intergenic regions)
CFG.PAR.gene_states_low_cover_cutoff = 10; % a value of 10 corresponds to a block
                                      % of ~2kb with mean coverage of 1 or 
                                      % a 1kb-region with 0-coverage
% restrict the range of splice site predictions to be > -inf in the IF and
% IL states to effectively enforce the canonical splice site dinucleotide
% consensus
CFG.PAR.enforce_splice_site_consensus = 1;
% compute two matrices states x features specifying a range of allowed
% feature values; if a feature value outside this range is encountered,
% decoding the corresponding state will not be possible. Hence, too many
% / too strict ranges can make decoding (& training) impossible!
CFG.PAR.perm_feature_ranges = set_permitted_feature_ranges(CFG);




%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT FILES
%%%%%%%%%%%%%%%%%%%%%%%%%


% directories where input data will be loaded from
rgasp_dir = '/fml/ag-raetsch/nobackup/projects/rgasp.2/';
genome_dir = [rgasp_dir 'genomes/drosophila/drosophila.gio'];
anno_dir = [rgasp_dir 'annotations/drosophila/'];
data_dir = '/fml/ag-raetsch/nobackup/projects/sequencing_runs/D_melanogaster/reads/';

% splice site predictions learned from read alignments directly
splice_site_dir = [data_dir 'D_melanogaster_L3.6.bestpaired.filtered.bam.splice_pred/'];
acc_splice_dir = [splice_site_dir 'splice_acc.bspf/pred/'];
don_splice_dir = [splice_site_dir 'splice_don.bspf/pred/'];

% check whether these exist
if exist(genome_dir, 'dir'),
  fprintf('Located GENOME_DIR (at %s)\n', genome_dir);
else
  error('Could NOT find GENOME_DIR (at %s)!', genome_dir);
end
if exist(anno_dir, 'dir'),
  fprintf('Located ANNO_DIR (at %s)\n', anno_dir);
else
  error('Could NOT find ANNO_DIR (at %s)!', anno_dir);
end

if exist(data_dir, 'dir'),
  fprintf('Located DATA_DIR (at %s)\n', data_dir);
else
  error('Could NOT find DATA_DIR (at %s)!', data_dir);
end
if exist(splice_site_dir, 'dir'),
  fprintf('Located SPLICE_SITE_DIR (at %s)\n', splice_site_dir);
else
  error('Could NOT find SPLICE_SITE_DIR (at %s)!', splice_site_dir);
end
if exist(acc_splice_dir, 'dir'),
  fprintf('Located ACC_SPLICE_DIR (at %s)\n', acc_splice_dir);
else
  error('Could NOT find ACC_SPLICE_DIR (at %s)!', acc_splice_dir);
end
if exist(don_splice_dir, 'dir'),
  fprintf('Located DON_SPLICE_DIR (at %s)\n', don_splice_dir);
else
  error('Could NOT find DON_SPLICE_DIR (at %s)!', don_splice_dir);
end


% set current reads directory here
read_map_dir  = data_dir;
read_map_file = [read_map_dir 'D_melanogaster_L3.6.bestpaired.filtered.bam'];
if isempty(read_map_file),
  error('No READ_MAP_FILE supplied!')
else
  if exist(read_map_file, 'file'),
    fprintf('Located READ_MAP_FILE (at %s)\n', read_map_file)
  else
    error('Could NOT find READ_MAP_FILE (at %s)!', read_map_file)
  end
end


% genome details
CFG.genome_dir = genome_dir;
CFG.genome_info = sprintf('%s/genome.config', genome_dir);
info = init_genome(CFG.genome_info);
CFG.chr_names = info.contig_names;
CFG.num_chr = length(CFG.chr_names);
for c=1:CFG.num_chr,
  info.flat_fnames{c} = [info.basedir '/genome/' info.contig_names{c} '.flat'];
  d = dir(info.flat_fnames{c});
  CFG.chr_lens(c) = d.bytes;
  assert(CFG.chr_lens(c)>0);
end


% annotation details
CFG.annotation_dir = sprintf('%s', anno_dir);
CFG.gene_fn = sprintf('%sgenes.mat', anno_dir);

% splice site prediction details
CFG.splice_site_dir.acc = acc_splice_dir;
CFG.splice_site_dir.don = don_splice_dir;


CFG.read_map_dir = read_map_dir;
CFG.read_map_file = read_map_file;

% competitor
CFG.cufflinks_pred_file = sprintf('%scufflinks/drosophila/D_melanogaster_L3.%s/sanitized_out/transcripts.gtf.mat', data_dir, '6.bestpaired.filtered');
CFG.cufflinks_convert = 0;

