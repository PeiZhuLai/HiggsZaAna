#!/bin/bash
/bin/hostname
gcc -v
pwd
# export PATH=$PATH:/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin/
source /cvmfs/cms.cern.ch/cmsset_default.sh

export PATH=/eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/bin:$PATH
source /eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/etc/profile.d/conda.sh

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
conda init

conda activate higgs-alp-ana


#### job

python ALP_plot_param.py -y run3 --ln -b

# python ALP_plot_param.py --tune

##########hep_sub runjob.sh -g cms -mem 8000 -wt mid -o job.out -e job.err
# python ALP_plot_param.py -m -y 2016 

#python ALP_Optimization.py -y run2 -o ./optimize_run2 --doOpt -c 5
#python ALP_plot_param.py -y run2 -m --ln

#python ALP_plot_param.py -y run2 -m -S #--ln #--cut --mA M30
#python ALP_plot_param.py -y run2 -m 
#python ALP_plot_bkgCorr.py -y run2 -m

#parser.add_argument("--region", help="0 for full region, 1 for signal region, 2 for sideband region")
# 0 and 4 cannot run parallely

# Modify your runjob.sh to include the following:
path_code='/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot'

export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/lib:/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/lib

if [ "$1" -eq 0 ]; then
    python3 $path_code/ALP_plot_param.py -y run3 --ln
elif [ "$1" -eq 1 ]; then
    python3 $path_code/ALP_plot_param.py -y run2 -m --ln
elif [ "$1" -eq 2 ]; then
    python3 $path_code/ALP_plot_param.py -y run2 -m --region 1 --ln
elif [ "$1" -eq 3 ]; then
    python3 $path_code/ALP_plot_param.py -y run2 -m --region 2 --ln
elif [ "$1" -eq 4 ]; then
    python3 $path_code/ALP_plot_param.py -y run2 -m
elif [ "$1" -eq 5 ]; then
    python3 $path_code/ALP_Optimization.py -y run2 -o $path_code/optimize_run3UL --region 0 -p --sigVSscore -s --doOpt -c 1
fi

# Optional: Commented-out commands
# python ALP_plot_param.py -y run2 -m --mu -S 
# python ALP_plot_param.py -y run2 -m -S 
# elif [ "$1" -eq 6 ]; then
#     #python ALP_
#python ALP_plot_param.py -y run2 -m --mu -S 
#python ALP_plot_param.py -y run2 -m -S 
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
#     python ALP_plot_param.py -y run2  -m --cut --mA M10 --ln
# elif [ $1 -eq 13 ]; then
#     python ALP_plot_param.py -y run2  -m --cut --mA M20 --ln
# elif [ $1 -eq 14 ]; then
#     python ALP_plot_param.py -y run2  -m --cut --mA M30 --ln
# fi