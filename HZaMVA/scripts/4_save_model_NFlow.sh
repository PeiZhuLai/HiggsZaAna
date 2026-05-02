#!/bin/bash

# After decide to use the model, copy to using folder

rm -fr ../using/model_Za_BDT_run3_NFlow.pkl
rm -fr ../using/Za-study_NFlow.db

cp model_Za_BDT_run3_NFlow.pkl ../using/
cp Za-study_NFlow.db ../using/