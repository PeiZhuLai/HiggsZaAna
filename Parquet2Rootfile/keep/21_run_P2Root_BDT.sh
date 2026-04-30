#!/bin/bash                                                                                                                                                                       
echo "==============STARTED=============="

# input="/eos/user/j/jiehan/parquet/nanov9/"
# target="/eos/home-j/jiehan/root/skimmed_ntuples_run2/"
# input="/eos/user/j/jiehan/parquet/nanov12/"
# target="/eos/home-j/jiehan/root/skimmed_ntuples_run3/"
# target="./"

# Run3
input="/eos/home-p/pelai/HZa/parquet_DNA/"
target="/eos/home-p/pelai/HZa/root_P2Root/run3_mergedBDT"

# input="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/Parquet/"
# target="/eos/home-p/pelai/HZa/Root_Dataset/run3/"

# Run2
# input="/eos/home-p/pelai/HZgamma/Parquet/NanoV9/run2/"
# target="/eos/home-p/pelai/HZgamma/Root_Dataset/run2/NanoV9/"


# years=(2022preEE 2022postEE 2023preBPix 2023postBPix 2024)
years=(2022preEE)
# years=(2022preEE 2022postEE)
# years=(2023preBPix 2023postBPix)
systs=("FNUF" "Material" "Electron_scale" "Electron_smear" "Muon_pt_scale" "Muon_pt_smear" "Photon_scale" "Photon_smear")
# systs=("FNUF" "Material" "Scale" "Smearing" "JER" "JES" "MET_JES" "MET_Unclustered" "Muon_pt")

ma_list=(1 2 3 4 5 6 7 8 9 10 15 20 25 30)

# 函数定义：执行命令并处理错误
execute_command() {
    local cmd="$1"
    local max_retries=10  # 最大重试次数
    local attempt=1

    while [ $attempt -le $max_retries ]; do
        echo "Attempt $attempt: $cmd"
        $cmd && break  # 如果命令成功执行，则跳出循环
        echo "Command $cmd failed. Retrying..."
        ((attempt++))
    done

    if [ $attempt -gt $max_retries ]; then
        echo "Error: Maximum retries reached. Command failed: $cmd"
    fi
}

# 函数定义：处理样本数据
process_sample() {
    local sample="$1"
    local type="$2"
    # 原本用 if/else 設定 corr，改為更簡潔的預設值寫法
    local corr="${3:-nominal}"
    
    # 修正條件語法: 使用 || 並確保中括號空白
    if [ "$type" = "Sig_MC" ] || [ "$type" = "Data" ]; then
        for year in "${years[@]}"; do
            command="python /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile/Parque2Root_BDT.py "
            command+="-i ${input}${type}/${sample}_${year}/merged_${corr}.parquet "
            if [ "$type" = "Data" ]; then
                command+="-o ${target}Data/${year}.root"
            else
                command+="-o ${target}${sample}/${year}.root"
            fi
            
            if [ "$type" = "Sig_MC" ]; then
                # 從樣本名擷取質量值 (e.g. ALP_M5 -> 5)
                ma_from_sample="${sample#*_M}"
                command+=" --split"
                command+=" --ma ${ma_from_sample}"
            fi

            execute_command "$command" &
            pid_list+=($!)
        done
    fi

    if [ "$type" = "Bkg_MC" ]; then
        for year in "${years[@]}"; do
            for ma in "${ma_list[@]}"; do
                command="python /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile/Parque2Root_BDT.py "
                command+="-i ${input}${type}/${sample}_${year}/merged_${corr}.parquet "
                command+="-o ${target}${sample}/${year}.root"
                command+=" --ma ${ma}"

                execute_command "$command" &
                pid_list+=($!)
            done
        done
    fi

    # 等待所有后台任务完成
    for pid in "${pid_list[@]}"; do
        wait $pid
    done

    echo "Sample $sample completed successfully."
}

# 函数定义：处理样本数据
process_sample_syst() {
    local sample="$1"
    local type="$2"
    local year="$3"
    local uod="$4"
    
    for syst in "${systs[@]}"; do
        corr="${syst}_${uod}"
        base_cmd="python /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile/Parque2Root_BDT.py "
        base_cmd+="-i ${input}${type}/${sample}_${year}/merged_${corr}.parquet "
        base_cmd+="-o ${target}${sample}_${syst}_${uod}/${year}.root"

        if [ "$type" = "Sig_MC" ]; then
            # 從樣本名擷取質量值
            ma_from_sample="${sample#*_M}"
            execute_command "$base_cmd --split --ma ${ma_from_sample}" &
            pid_list+=($!)
        elif [ "$type" = "Bkg_MC" ]; then
            # 若未來啟用 Bkg_MC 系統誤差流程，逐一掃描 ma_list
            for ma in "${ma_list[@]}"; do
                execute_command "$base_cmd --ma ${ma}" &
                pid_list+=($!)
            done
        else
            execute_command "$base_cmd" &
            pid_list+=($!)
        fi
    done

    # 等待所有后台任务完成
    for pid in "${pid_list[@]}"; do
        wait $pid
    done

    echo "Sample $sample completed successfully."
}


# 处理 signal 样本

# samples=(ggH_M125 VBF_M125 WplusH_M125 WminusH_M125 ZH_M125 ttH_M125 ggH_M120 VBFH_M120 WplusH_M120 WminusH_M120 ZH_M120 ttH_M120 ggH_M130 VBFH_M130 WplusH_M130 WminusH_M130 ZH_M130 ttH_M130 ggH_mix VBF_mix ggH VBF WplusH WminusH ZH ttH)
# samples=(ALP_M5 ALP_M15 ALP_M30)

samples=(ALP_M5 ALP_M15 ALP_M30)
type="Sig_MC"
for sample in "${samples[@]}"; do
    mkdir -p "$target${sample}/"
    # 存储后台任务的进程ID列表
    pid_list=()

    # 调用函数处理样本数据
    process_sample "$sample" "$type"
done

samples=(ALP_M5 ALP_M15 ALP_M30)
type="Sig_MC"
for sample in "${samples[@]}"; do
    for sf in "up" "down"; do #  "up" "down"
        for syst in "${systs[@]}"; do
            mkdir -p "$target${sample}_${syst}_${sf}"
        done
        for year in "${years[@]}"; do
            # 存储后台任务的进程ID列表
            pid_list=()

            # 调用函数处理样本数据
            process_sample_syst "$sample" "$type" "$year" "$sf"
        done
    done
done

# ****************************
# ********** Bkg *************
# ****************************

# ****************************
# ********* Nomianl **********
# ****************************

#  Run 2
#       Z + Fake g   |   Z + g          |   tt        |   tg/ttg  |   VBS Z + g   |   tt + X  |   Multibosons              |
#       DYJetsToLL   |   ZGToLLG        |   TTtoLNu2Q |   TGJets  |   ZG2JToG2L2J |   ttWJets |   WW  / WZ  / ZZ  / WGToLNuG  |
#       EWKZ2J       |                  |   TT        |   TTGJets |               |   ttZJets |   WWG / WZG / ZZG /

#  Run 3
#  2022 DYJetsToLL   |DYGto2LG_10to50   |   TTtoLNu2Q |   Lack    |   ZG2JToG2L2J |   Lack    |   WW / WZ / ZZ
#  2022EE            |DYGto2LG_50to100  |   TT        |   Lack    |               |   Lack    |

#  2023 DYJetsToLL   |DYGto2LG_10to100  |   TTtoLNu2Q |   Lack    |   ZG2JToG2L2J |   Lack    |   WW / WZ / ZZ 
#                    |                  |   TT        |   Lack    |               |   Lack    |

# 处理 bkgmc 样本

# samples=(ZGToLLG DYJetsToLL WGToLNuG ZG2JToG2L2J EWKZ2J TT TTGJets TGJets ttWJets ttZJets WW WZ ZZ DYGto2LG_10to50 DYGto2LG_50to100)
for i in {1..3};do
    if [ "$i" = "1" ]; then
        samples=(DYGto2LG_10to50 DYGto2LG_50to100)
        years=(2022preEE 2022postEE)

    elif [ "$i" = "2" ]; then
        samples=(DYGto2LG_10to100)
        years=(2023preBPix 2023postBPix)

    elif [ "$i" = "3" ]; then
        samples=(DYJetsToLL)
        years=(2022preEE 2022postEE 2023preBPix 2023postBPix)
    fi

    # samples=(DYGto2LG_10to100)
    # samples=(DYJetsToLL)
    type="Bkg_MC"
    for sample in "${samples[@]}"; do
        mkdir -p "$target$sample"
        # 存储后台任务的进程ID列表
        pid_list=()

        # 调用函数处理样本数据
        process_sample "$sample" "$type"
    done
done 

# samples=(ZGToLLG DYJetsToLL WGToLNuG ZG2JToG2L2J EWKZ2J TT TTGJets TGJets ttWJets ttZJets WW WZ ZZ DYGto2LG_10to50 DYGto2LG_50to100)
# samples=(DYGto2LG_10to50 DYGto2LG_50to100)
# samples=(DYGto2LG_10to100)
# samples=(DYJetsToLL)
# type="Bkg_MC"
# for sample in "${samples[@]}"; do
#     for sf in "down"; do #  "up" "down"
#         for syst in "${systs[@]}"; do
#             mkdir -p "$target${sample}_${syst}_${sf}"
#         done
#         for year in "${years[@]}"; do
#             # 存储后台任务的进程ID列表
#             pid_list=()

#             # 调用函数处理样本数据
#             process_sample_syst "$sample" "$type" "$year" "$sf"
#         done
#     done
# done

# ****************************
# ********** Data ************
# ****************************

# ****************************
# ********* Nomianl **********
# ****************************

# # 处理 data 样本

# samples=(Data)
# type="Data"
# for sample in "${samples[@]}"; do
#     mkdir -p "$target$sample"
#     # 存储后台任务的进程ID列表
#     pid_list=()

#     # 调用函数处理样本数据
#     process_sample "$sample" "$type"
# done

# Use fake photon background estimation with data-driven

# mkdir -p /eos/home-j/jiehan/root/2017/skimmed_ntuples/data_med/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/data_fake/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/mc_true/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/mc_med/
# python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/apply_weight.py

# ######################
# Non prompt MC sample
# ######################

echo "==============FINISHED==========="
