function CFG = elegans_settings(CFG)

% This script prepares paths and input/output files specialized
% for elegans test.


% parameter combinations to be used for independent training
% (typically overwritten)
CFG.train_params = { ...
  [0.1],   [0.5],   [0.05],  [750], 'QP', [1]; ...  % max amount of train data available without val data
%  [0.01],   [0.05],   [5],  [750], 'QP', [1]; ...  % max amount of train data available without val data
%  [0.01],   [2.5],   [0.05],  [750], 'QP', [1]; ...  % max amount of train data available without val data
%  [2.5],   [0.5],   [0.05],  [750], 'QP', [1]; ...  % max amount of train data available without val data
};



%%%%%%%%%%%%%%%%%%%%%%%%%
% Reset features
%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters for generating the intron span feature
CFG.read_intron_span_max_intron_len = 200000; % TODO sufficient for worm&fly
CFG.read_intron_span_min_score = 15; % C. elegans & D. melanogaster
CFG.read_intron_span_max_mismatches = 1; 

% parameters for extracting splice site features from spliced read alignments
CFG.read_splice_site_max_intron_len = 200000; % TODO sufficient for worm&fly
CFG.read_splice_site_min_score = 10; 
CFG.read_splice_site_max_mismatches = 1; 

% parameters for retrieving read pair features from read alignments
CFG.read_pair_cover_max_intron_len = 200000;
CFG.read_pair_cover_min_score = 15;
CFG.read_pair_cover_max_mismatches = 1;

% threshold on total read coverage used to define low-coverage regions
CFG.low_cover_threshold = 5;  


CFG.max_train_chunk_len = 100000;
% number of discrete expression levels

%warning('ONLY 1EXPRESSION SUBMODEL!');
CFG.PAR.num_levels = 1; % something like 3-10;
% number of supporting points for each scoring function
%CFG.PAR.num_plif_nodes = 2*CFG.PAR.num_levels; %CFG.PAR.num_levels+2 %2*CFG.PAR.num_levels;
CFG.PAR.num_plif_nodes = 30; %CFG.PAR.num_levels+2 %2*CFG.PAR.num_levels;
% include mate-pair information for paired reads as a learning features
CFG.PAR.use_pair_feats = 0;
% include data on genomic repeats as a learning features
CFG.PAR.use_repeat_feats = 0;
% enforce the monotonicity for some feature scoring functions
% note that this can make the training problem more difficult
% (i.e. time-consuming) to solve!
CFG.PAR.enf_monot_score_funcs = 0;
% use heuristic training procedure
CFG.PAR.constraint_margin = 5;
% options for enforcing a certain  state given a certain range of feature values

% use splice feats
CFG.PAR.use_splice_feats = 1;
% use original filtered intron span feature
CFG.PAR.use_filtered_intron_feats = 1;
% use the binned intron span features
CFG.PAR.use_binned_span_feats = 0;
% use the cufflinks feature
CFG.PAR.use_cuffl_feats = 0;


%warning('disable intron span feature in ige state.');
CFG.PAR.ige_intron_span = 0;

%fprintf('Elegans test withOUT permitted feature ranges.\n');
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
%warning('No permitted feature ranges');
%CFG.PAR = rmfield(CFG.PAR,'perm_feature_ranges');


CFG.grid_memreq = 16;
%warning('For debugging purposes turning off grid usage.');
%CFG.grid_use = 0;


% bundle method settings
warning('Extremly tiny duality gap');
CFG.PAR.bmrm.useMonotConstr = 0;  % use monotonicity constraints (0=no 1=yes)
CFG.PAR.bmrm.eps = 1e-6;          % relative gap
CFG.PAR.bmrm.useAggregation = 1;  % use aggregation (set appropriate K)
CFG.PAR.bmrm.K = 168;             % number of cutting planes used -1 + aggregated plane
CFG.PAR.bmrm.checkRealObjMod = 5; % heuristic 1
CFG.PAR.bmrm.monotMargin = 0.01;  % heuristic 2
% bmrm linesearch settings
CFG.PAR.bmrm.useLs = 1;           % use linesearch (0=no 1=yes)
CFG.PAR.bmrm.lsEps = 1e-3;        % break if parabola approx gap is smaller than lsEps
CFG.PAR.bmrm.lsSteps = 50;        % max. number of iterations
CFG.PAR.bmrm.lsTheta = 0.25;      % w_new = (1-lsTheta)*w_ls + lsTheta*w_old;

warning('USING HMM AS CLASSIFIER.')
CFG.PAR.classifier = @train_hmm;

% dont filter label noise examples out of the training set
CFG.PAR.train_filter_label_noise = 0;
CFG.PAR.train_filter_amb_label = 0;
CFG.PAR.train_filter_ige_cover = 0;




%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT FILES
%%%%%%%%%%%%%%%%%%%%%%%%%


% directories where input data will be loaded from
data_dir = '/home/nico/mtim/elegans/';
genome_dir = [data_dir 'elegans.gio'];
anno_dir = [data_dir];

% splice site predictions learned from read alignments directly
splice_site_dir = [data_dir];
acc_splice_dir = [splice_site_dir 'acc/'];
don_splice_dir = [splice_site_dir 'don/'];

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
read_map_file = [read_map_dir 'elegans.6.bam'];
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
%CFG.cufflinks_pred_file = sprintf('%selegans.%s.cufflinks.mat', data_dir, '6');
CFG.cufflinks_convert = 1;


