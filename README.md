# CERENKOV3 

CERENKOV3 is machine learning pipeline which constructs molecular networks, creates new features from tuning Node2vec on those networks, and runs a customized, XGBoost-based classifier to identify regulatory SNPs in the non-coding human genome.

CERENKOV3 was created by Yao Yao and Stephen Ramsey at Oregon State University, and will be published at Pacific Symposium on Biocomputing (PSB) 2020.

The pipeline consists of 3 parts:

1. The construction of molecular network. It's carried out mainly by the scripts in [cerenkov3_data](https://github.com/ramseylab/cerenkov3/tree/master/cerenkov3_data) and the basic network will be stored as multiple edge list.
2. The customized classifier, [`cerenkov3_classifier`](https://github.com/ramseylab/cerenkov3/blob/master/cerenkov3_classifier/cerenkov3_classifier.py), which is responsible to:
    - accept the weights of edge types as hyperparameters and run Node2vec on the dynamically weighted network to gain new featurs;
    - accept the annotation features from CERENKOV and the new clustering-based feature, _Locus Size_
    - initialize the base classifier (XGBoost in our case) to train and test
3. Hyperparamter optimaztion, training/testing, performance comparision and visualization

Note that `cerenkov3_classifier` is not fully scikit-learn compatible because Node2vec only accept networks whose node IDs are integers and threfore to simplify the implementation, we designed `cerenkov3_classifier` to accept integer IDs instead of the whole feature matrix.
