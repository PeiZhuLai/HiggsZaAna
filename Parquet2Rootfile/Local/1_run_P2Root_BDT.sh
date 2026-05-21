#!/bin/bash                                                                                                                                                                       
timer_start=$(date +%s)
echo "==============TIMER STARTED: $(date '+%F %T')=============="
print_runtime() {
    local timer_end elapsed hours minutes seconds
    timer_end=$(date +%s)
    elapsed=$((timer_end - timer_start))
    hours=$((elapsed / 3600))
    minutes=$(((elapsed % 3600) / 60))
    seconds=$((elapsed % 60))
    printf "==============RUNTIME: %02d:%02d:%02d==============\n" "$hours" "$minutes" "$seconds"
}
trap print_runtime EXIT

echo "==============STARTED=============="


# Run3
# input="/eos/home-p/pelai/HZa/parquet_DNA/"
input="/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA/"
target="/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_nominal/"

# Gen Info
# input="/eos/home-p/pelai/HZa/parquet_Sig_MC_DNA/"
# target="/eos/home-p/pelai/HZa/root_Sig_MC_P2Root/run3_mergedBDT/"

# sum Cut Study
# input="/eos/home-p/pelai/HZa/parquet_sumStudy_DNA/"
# target="/eos/home-p/pelai/HZa/root_P2Root/run3_sumStudy_mergedBDT/"


# input="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/Parquet/"
# target="/eos/home-p/pelai/HZa/root_P2Root/run3_mergedBDT/"


years=(2022preEE 2022postEE 2023preBPix 2023postBPix 2024)
# years=(2023preBPix) # Sum Cut Study
# years=(2022preEE 2022postEE 2023preBPix 2023postBPix)
# years=(2022preEE)
# years=(2022preEE 2022postEE)
# years=(2023preBPix 2023postBPix)
systs=("FNUF" "Material" "Electron_scale" "Electron_smear" "Muon_scale" "Muon_smear" "Photon_scale" "Photon_smear")
# systs=("FNUF" "Material" "Scale" "Smearing" "JER" "JES" "MET_JES" "MET_Unclustered" "Muon_pt")
MAX_PARALLEL_JOBS=6

# 函数定义：执行命令并处理错误
execute_command() {
    local cmd="$1"
    local max_retries=2  # 最大重试次数
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

wait_for_job_slot() {
    while [ "$(jobs -pr | wc -l)" -ge "$MAX_PARALLEL_JOBS" ]; do
        wait -n || true
    done
}

submit_background_job() {
    local cmd="$1"
    wait_for_job_slot
    execute_command "$cmd" &
}

wait_for_background_jobs() {
    local label="$1"
    while [ "$(jobs -pr | wc -l)" -gt 0 ]; do
        wait -n || true
    done
    echo "$label completed successfully."
}

# 函数定义：处理样本数据
process_sample() {
    local sample="$1"
    local type="$2"
    if [ -z "$3" ]; then
        local corr="nominal"
    else
        local corr="$3"
    fi
    
    for year in "${years[@]}"; do
        command="python /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile/Parque2Root_BDT.py "
        command+="-i ${input}${type}/${sample}_${year}/merged_${corr}.parquet "
        if [ "$type" = "Data" ]; then
            command+="-o ${target}Data/${year}.root"
        else
            command+="-o ${target}${sample}/${year}.root"
        fi
        
        if [ "$type" = "Sig_MC" ]; then
            command+=" --split"
        fi

        # 使用函数执行命令
        submit_background_job "$command"
    done

    echo "Sample $sample submitted."
}

# 函数定义：处理样本数据
process_sample_syst() {
    local sample="$1"
    local type="$2"
    local year="$3"
    local uod="$4"
    
    for syst in "${systs[@]}"; do
        corr="${syst}_${uod}"
        command="python /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile/Parque2Root_BDT.py "
        command+="-i ${input}${type}/${sample}_${year}/merged_${corr}.parquet "
        command+="-o ${target}${sample}_${syst}_${uod}/${year}.root"
        
        if [ "$type" = "Sig_MC" ]; then
            command+=" --split"
        fi

        # 使用函数执行命令
        submit_background_job "$command"
    done

    echo "Sample $sample $year $uod systematics submitted."
}

#---------------------------------------------------------------------------------------
## 处理 signal 样本
samples=(mA_M1 mA_M2 mA_M3 mA_M4 mA_M5 mA_M6 mA_M7 mA_M8 mA_M9 mA_M10 mA_M15 mA_M20 mA_M25 mA_M30)
# samples=(mA_M7 mA_M8)
# samples=(mA_M5)
type="Sig_MC"
for sample in "${samples[@]}"; do
    mkdir -p "$target${sample}/"

    # 调用函数处理样本数据
    process_sample "$sample" "$type"
done
wait_for_background_jobs "Signal nominal jobs"

samples=(mA_M1 mA_M2 mA_M3 mA_M4 mA_M5 mA_M6 mA_M7 mA_M8 mA_M9 mA_M10 mA_M15 mA_M20 mA_M25 mA_M30)
# samples=(mA_M7 mA_M8)
type="Sig_MC"
for sample in "${samples[@]}"; do
    for sf in "up" "down"; do #  "up" "down"
        for syst in "${systs[@]}"; do
            mkdir -p "$target${sample}_${syst}_${sf}"
        done
        for year in "${years[@]}"; do
            # 调用函数处理样本数据
            process_sample_syst "$sample" "$type" "$year" "$sf"
        done
    done
done
wait_for_background_jobs "Signal systematic jobs"

# ****************************
# ********** Bkg *************
# ****************************

# ****************************
# ********* Nomianl **********
# ****************************

#  Run 3
# -------------------------------------
#  2022preEE, 2022postEE
#  DYJetsToLL   |DYGto2LG_10to50   
#               |DYGto2LG_50to100
# -------------------------------------
#  2023preBPix, 2023postBPix
#  DYJetsToLL   |DYGto2LG_10to100  
# -------------------------------------
#  2024
#  DYJetsTo2E   |DYGto2LG_10to100
#  DYJetsTo2Mu
#  DYJetsTo2Tau

# ## 处理 bkgmc 样本
for i in {1..4};do
# for i in {2..3};do # Sum Cut Study
# for i in 2; do # Compliment Study
    if [ "$i" = "1" ]; then
        samples=(DYGto2LG_10to50 DYGto2LG_50to100)
        years=(2022preEE 2022postEE)

    elif [ "$i" = "2" ]; then
        samples=(DYGto2LG_10to100)
        years=(2023preBPix 2023postBPix 2024)
        # years=(2023preBPix) # Sum Cut Study

    elif [ "$i" = "3" ]; then
        samples=(DYJetsToLL)
        years=(2022preEE 2022postEE 2023preBPix 2023postBPix)
        # years=(2023preBPix) # Sum Cut Study

    elif [ "$i" = "4" ]; then
        samples=(DYJetsTo2E DYJetsTo2Mu DYJetsTo2Tau)
        years=(2024)
    fi

    type="Bkg_MC"
    for sample in "${samples[@]}"; do
        mkdir -p "$target$sample"

        # 调用函数处理样本数据
        process_sample "$sample" "$type"
    done
done 
wait_for_background_jobs "Background jobs"

# ****************************
# ********** Data ************
# ****************************

# ****************************
# ********* Nomianl **********
# ****************************

# 处理 data 样本
samples=(Data)
years=(2022preEE 2022postEE 2023preBPix 2023postBPix 2024)
type="Data"
for sample in "${samples[@]}"; do
    mkdir -p "$target$sample"

    # 调用函数处理样本数据
    process_sample "$sample" "$type"
done
wait_for_background_jobs "Data jobs"

# Use fake photon background estimation with data-driven
# mkdir -p /eos/home-j/jiehan/root/2017/skimmed_ntuples/data_med/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/data_fake/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/mc_true/ /eos/home-j/jiehan/root/2017/skimmed_ntuples/mc_med/
# python /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/scripts/apply_weight.py

# ######################
# Non prompt MC sample
# ######################

echo "==============FINISHED==========="
