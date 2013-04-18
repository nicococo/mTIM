mTIM
====

transcript identification method


A) PREPARATION
--------------

- software package comes without parallelization software package 
  but training and prediction are ready for it
  - see general_settings.m : CFG.grid_use = 0; 

- for generating mex-files you need to adapt the 
  makefiles in the corresponding directories (just change 
  the path to your matlab)
   mTIM/src/utils : change Makefile, then type 'make all'
   mTIM/src/sotool/native  : change Makefile, then type 'make all'
   mTIM/src/sotool/losses  : change Makefile, then type 'make all'
   mTIM/src/model : just type 'mex convert_states2labels.cpp' in your matlab command line

B) TESTING
----------
- get the mini test example which is a driectory consisting of
  - acceptor and donor splice-sites
  - genome
  - annotation
  - rna-seq reads
- you have to manually adapt the path in mini.gio/genome.config
- you also need to adapt the paths in mTIM/src/mini_settings.m

- generate data by typing 'mTIM_prepare_data('mini')'
  - you can find the generated files under mTIM/out/mini
  
- training and prediction with 'mTIM_predict('mini')'
  - training output files are stored in mTIM/out/mini/{DATE}/
    where the trained predictor is in mTIM/out/mini/{DATE}/xval_fold{X}/sosvm_final.mat
  - predicted annotation is in mTIM/out/mini/{DATE}/prediction.mat
  
MISC
----
- feature generation files are in /src/data_preparation
- viterbi is in /src/sotool/native/best_path.cpp