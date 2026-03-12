#!/bin/bash
pwd
#------------------------------Dry Run---------------------------------------
#1
# python scripts/plot_variable_dataVmc.py -y run3 --ln -b
#2
# python3 scripts/plot_variable_dataVmc.py -y run3 -m --ln -b
#3
# python3 scripts/plot_variable_dataVmc.py -y run3 -m --region 1 --ln
#4
# python3 scripts/plot_variable_dataVmc.py -y run3 -m --region 2 --ln
#5
# python3 scripts/ALP_Optimization.py -y run3 -o ./optimize_run3UL --region 0 -p --sigVSscore -s --doOpt -c 1
#----------------------------------------------------------------------------

# 0: Full Region, 1: Signal Region, 2: Contral Region
#-----------------------------Parallel Run-----------------------------------
path_code='/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/scripts'
path_output='/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots'
export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/lib:/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/lib


if [ "$1" -eq 0 ]; then
    python3 $path_code/plot_variable_dataVmc.py -y run3 --ln -b
elif [ "$1" -eq 1 ]; then
    python3 $path_code/plot_variable_dataVmc.py -y run3 -m --ln -b
elif [ "$1" -eq 2 ]; then
    python3 $path_code/plot_variable_dataVmc.py -y run3 -m --region 1 --ln
elif [ "$1" -eq 3 ]; then
    python3 $path_code/plot_variable_dataVmc.py -y run3 -m --region 2 --ln
elif [ "$1" -eq 4 ]; then
    python3 $path_code/ALP_Optimization.py -y run3 -o $path_output/optimize_run3UL --region 1 -p --sigVSscore -s --doOpt -c 2
elif [ "$1" -eq 5 ]; then
    python3 $path_code/ALP_Optimization.py -y run3 -o $path_output/optimize_run3UL --region 1 -p --sigVSscore -s --doOpt -c 1
fi
#----------------------------------------------------------------------------

# python plot_variable_dataVmc.py --tune
# python plot_variable_dataVmc.py -m -y 2016 
# python ALP_Optimization.py -y run2 -o ./optimize_run2 --doOpt -c 5
# python plot_variable_dataVmc.py -y run2 -m --ln
# python plot_variable_dataVmc.py -y run2 -m -S #--ln #--cut --mA M30
# python plot_variable_dataVmc.py -y run2 -m 
# python ALP_plot_bkgCorr.py -y run2 -m


# Optional: Commented-out commands
# python plot_variable_dataVmc.py -y run2 -m --mu -S 
# python plot_variable_dataVmc.py -y run2 -m -S 
# elif [ "$1" -eq 6 ]; then
#     #python ALP_
#python plot_variable_dataVmc.py -y run2 -m --mu -S 
#python plot_variable_dataVmc.py -y run2 -m -S 
# elif [ $1 -eq 6 ]; then
#     #python ALP_Optimization.py -y run2 -o ./optimize_run2UL -p --sigVSscore -s --ele --doOpt -c 1
#     python ALP_plot_bkgCorr.py -y run2 -m
# elif [ $1 -eq 7 ]; then
#     python ALP_Optimization.py -y run2 -o ./optimize_run2UL -p --sigVSscore -s --mu --doOpt -c 1
#     python ALP_plot_bkgCorr.py -y run2 -m --mu
# elif [ $1 -eq 8 ]; then
#     #python ALP_Optimization.py -y run2 -o ./optimize_run2UL_ele -p --sigVSscore -s --ele --doOpt -c 2
#     python ALP_BDTSys.py -y $2 -m --mu
# elif [ $1 -eq 9 ]; then
#     #python ALP_Optimization.py -y run2 -o ./optimize_run2UL_mu -p --sigVSscore -s --mu --doOpt -c 2
#     python ALP_BDTSys.py -y $2 -m --ele
# elif [ $1 -eq 10 ]; then
#     python ALP_NormalizationSys_param.py -y run2 -m --ele
#     python ALP_NormalizationSys_param.py -y $2 -m --mu
# elif [ $1 -eq 11 ]; then
#     python ALP_NormalizationSys_param.py -y $2 -m --ele
#     #python ALP_BDTSys.py -y run2 --mu
# elif [ $1 -eq 12 ]; then
#     python plot_variable_dataVmc.py -y run2  -m --cut --mA M10 --ln
# elif [ $1 -eq 13 ]; then
#     python plot_variable_dataVmc.py -y run2  -m --cut --mA M20 --ln
# elif [ $1 -eq 14 ]; then
#     python plot_variable_dataVmc.py -y run2  -m --cut --mA M30 --ln
# fi

# export PATH=$PATH:/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin/
# source /cvmfs/cms.cern.ch/cmsset_default.sh

# export PATH=/eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/bin:$PATH
# source /eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/etc/profile.d/conda.sh

# # >>> conda initialize >>>
# # !! Contents within this block are managed by 'conda init' !!
# __conda_setup="$('/eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
# if [ $? -eq 0 ]; then
#     eval "$__conda_setup"
# else
#     if [ -f "/eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/etc/profile.d/conda.sh" ]; then
#         . "/eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/etc/profile.d/conda.sh"
#     else
#         export PATH="/eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/bin:$PATH"
#     fi
# fi
# unset __conda_setup
# # <<< conda initialize <<<
# # conda init

# conda activate higgs-alp-ana
