#!/bin/bash  

# Training
rm -fr model_Za_BDT_run3.pkl
rm -fr Za-study.db

python run3_Za_BDT.py

# After decide to use the model, copy to using folder
# rm -fr ../using/model_Za_BDT_run3.pkl
# rm -fr ../using/Za-study.db

# cp model_Za_BDT_run3.pkl ../using/
# cp Za-study.db ../using/