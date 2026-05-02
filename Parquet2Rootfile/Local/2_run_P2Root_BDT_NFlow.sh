#!/bin/bash
echo "==============STARTED=============="

# Run3 NFlow BDT scoring
input="/eos/home-p/pelai/HZa/parquet_DNA_NFlow/"
dataInput="/eos/home-p/pelai/HZa/parquet_DNA/"
target="/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_NFlow/"

years=(2022preEE 2022postEE 2023preBPix 2023postBPix 2024)
systs=("FNUF" "Material" "Electron_scale" "Electron_smear" "Muon_scale" "Muon_smear" "Photon_scale" "Photon_smear")

execute_command() {
    local cmd="$1"
    local max_retries=2
    local attempt=1

    while [ $attempt -le $max_retries ]; do
        echo "Attempt $attempt: $cmd"
        $cmd && break
        echo "Command $cmd failed. Retrying..."
        ((attempt++))
    done

    if [ $attempt -gt $max_retries ]; then
        echo "Error: Maximum retries reached. Command failed: $cmd"
    fi
}

process_sample() {
    local sample="$1"
    local type="$2"
    local input_base="$input"
    if [ -z "$3" ]; then
        local corr="nominal"
    else
        local corr="$3"
    fi

    if [ "$type" = "Data" ]; then
        input_base="$dataInput"
    fi

    for year in "${years[@]}"; do
        command="python /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile/Parque2Root_BDT_NFlow.py "
        command+="-i ${input_base}${type}/${sample}_${year}/merged_${corr}.parquet "
        if [ "$type" = "Data" ]; then
            command+="-o ${target}Data/${year}.root"
        else
            command+="-o ${target}${sample}/${year}.root"
        fi

        if [ "$type" = "Sig_MC" ]; then
            command+=" --split"
        fi

        execute_command "$command" &
        pid_list+=($!)
    done

    for pid in "${pid_list[@]}"; do
        wait $pid
    done

    echo "Sample $sample completed successfully."
}

process_sample_syst() {
    local sample="$1"
    local type="$2"
    local year="$3"
    local uod="$4"

    for syst in "${systs[@]}"; do
        corr="${syst}_${uod}"
        command="python /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile/Parque2Root_BDT_NFlow.py "
        command+="-i ${input}${type}/${sample}_${year}/merged_${corr}.parquet "
        command+="-o ${target}${sample}_${syst}_${uod}/${year}.root"

        if [ "$type" = "Sig_MC" ]; then
            command+=" --split"
        fi

        execute_command "$command" &
        pid_list+=($!)
    done

    for pid in "${pid_list[@]}"; do
        wait $pid
    done

    echo "Sample $sample completed successfully."
}

# Signal nominal samples
samples=(mA_M1 mA_M2 mA_M3 mA_M4 mA_M5 mA_M6 mA_M7 mA_M8 mA_M9 mA_M10 mA_M15 mA_M20 mA_M25 mA_M30)
type="Sig_MC"
for sample in "${samples[@]}"; do
    mkdir -p "$target${sample}/"
    pid_list=()
    process_sample "$sample" "$type"
done

# Signal systematics
samples=(mA_M1 mA_M2 mA_M3 mA_M4 mA_M5 mA_M6 mA_M7 mA_M8 mA_M9 mA_M10 mA_M15 mA_M20 mA_M25 mA_M30)
type="Sig_MC"
for sample in "${samples[@]}"; do
    for sf in "up" "down"; do
        for syst in "${systs[@]}"; do
            mkdir -p "$target${sample}_${syst}_${sf}"
        done
        for year in "${years[@]}"; do
            pid_list=()
            process_sample_syst "$sample" "$type" "$year" "$sf"
        done
    done
done

# Background samples
for i in {1..4}; do
    if [ "$i" = "1" ]; then
        samples=(DYGto2LG_10to50 DYGto2LG_50to100)
        years=(2022preEE 2022postEE)
    elif [ "$i" = "2" ]; then
        samples=(DYGto2LG_10to100)
        years=(2023preBPix 2023postBPix 2024)
    elif [ "$i" = "3" ]; then
        samples=(DYJetsToLL)
        years=(2022preEE 2022postEE 2023preBPix 2023postBPix)
    elif [ "$i" = "4" ]; then
        samples=(DYJetsTo2E DYJetsTo2Mu DYJetsTo2Tau)
        years=(2024)
    fi

    type="Bkg_MC"
    for sample in "${samples[@]}"; do
        mkdir -p "$target$sample"
        pid_list=()
        process_sample "$sample" "$type"
    done
done

# Data samples use nominal, non-NFlow parquet because Data has no corrected branches.
samples=(Data)
years=(2022preEE 2022postEE 2023preBPix 2023postBPix 2024)
type="Data"
for sample in "${samples[@]}"; do
    mkdir -p "$target$sample"
    pid_list=()
    process_sample "$sample" "$type"
done

echo "==============FINISHED==========="
