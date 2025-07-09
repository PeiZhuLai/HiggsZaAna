# # !/bin/bash                                                                                                                                                                       

#########################################################################
# # Prepare annoying DYGto2LG
#########################################################################
# DYGto2LG_path=/eos/home-p/pelai/HZa/Root_Dataset/run3/DYGto2LG
# DYGto2LG_10to50_path=/eos/home-p/pelai/HZa/Root_Dataset/run3/DYGto2LG_10to50
# DYGto2LG_50to100_path=/eos/home-p/pelai/HZa/Root_Dataset/run3/DYGto2LG_50to100
# DYGto2LG_10to100_path=/eos/home-p/pelai/HZa/Root_Dataset/run3/DYGto2LG_10to100

# # hadd total.root file1.root file2.root
# if [ -d "$DYGto2LG_path" ]; then
#     echo "Directory exists: $DYGto2LG_path — removing it."
#     rm -rf "$DYGto2LG_path"
# else
#     echo "Directory does not exist: $DYGto2LG_path — creating it."
# fi
# mkdir -p "$DYGto2LG_path"

# echo "hadd $DYGto2LG_path/2022preEE.root $DYGto2LG_10to50_path/2022preEE.root $DYGto2LG_50to100_path/2022preEE.root"
# hadd $DYGto2LG_path/2022preEE.root $DYGto2LG_10to50_path/2022preEE.root $DYGto2LG_50to100_path/2022preEE.root
# echo "hadd $DYGto2LG_path/2022postEE.root $DYGto2LG_10to50_path/2022postEE.root $DYGto2LG_50to100_path/2022postEE.root"
# hadd $DYGto2LG_path/2022postEE.root $DYGto2LG_10to50_path/2022postEE.root $DYGto2LG_50to100_path/2022postEE.root
# # cp from.root to.root
# echo "cp $DYGto2LG_10to100_path/2023preBPix.root $DYGto2LG_path/2023preBPix.root"
# cp $DYGto2LG_10to100_path/2023preBPix.root $DYGto2LG_path/2023preBPix.root
# echo "cp $DYGto2LG_10to100_path/2023postBPix.root $DYGto2LG_path/2023postBPix.root"
# cp $DYGto2LG_10to100_path/2023postBPix.root $DYGto2LG_path/2023postBPix.root
# # Sum four eras
# years=( 2022preEE 2022postEE 2023preBPix 2023postBPix )

# input_files=""
# for year in "${years[@]}"; do
#     input_files+=" $DYGto2LG_path/${year}.root"
# done

# output_file="$DYGto2LG_path/run3.root"
# echo "hadd $output_file $input_files"
# hadd $output_file $input_files

#########################################################################
# # Prepare DYJetsToLL
#########################################################################
# DYJetsToLL_path=/eos/home-p/pelai/HZa/Root_Dataset/run3/DYJetsToLL
# years=( 2022preEE 2022postEE 2023preBPix 2023postBPix )

# input_files=""
# for year in "${years[@]}"; do
#     input_files+=" $DYJetsToLL_path/${year}.root"
# done

# output_file="$DYJetsToLL_path/run3.root"
# echo "hadd $output_file $input_files"
# hadd $output_file $input_files

#########################################################################
# Add run3.root all background
#########################################################################
# DYJetsToLL_path=/eos/home-p/pelai/HZa/Root_Dataset/run3/DYJetsToLL
# DYGto2LG_path=/eos/home-p/pelai/HZa/Root_Dataset/run3/DYGto2LG
# Bkg_MC_path=/eos/home-p/pelai/HZa/Root_Dataset/run3/All_Bkg
# if [ -d "$Bkg_MC_path" ]; then
#     echo "Directory exists: $Bkg_MC_path — removing it."
#     rm -rf "$Bkg_MC_path"
# else
#     echo "Directory does not exist: $Bkg_MC_path — creating it."
# fi
# mkdir -p "$Bkg_MC_path"

# echo "hadd $Bkg_MC_path/run3.root $DYJetsToLL_path/run3.root $DYGto2LG_path/run3.root"
# hadd $Bkg_MC_path/run3.root $DYJetsToLL_path/run3.root $DYGto2LG_path/run3.root


#########################################################################
# # Prepare Data
#########################################################################
# Data_path=/eos/home-p/pelai/HZa/Root_Dataset/run3/Data
# years=( 2022preEE 2022postEE 2023preBPix 2023postBPix )

# input_files=""
# for year in "${years[@]}"; do
#     input_files+=" $Data_path/${year}.root"
# done

# output_file="$Data_path/run3.root"
# echo "hadd $output_file $input_files"
# hadd $output_file $input_files


#########################################################################
# # Prepare Sig
#########################################################################
base_path=/eos/home-p/pelai/HZa/Root_Dataset/run3
massList=( M5 M15 M30 )
years=( 2022preEE )

# Add years into run3.root
for mass in "${massList[@]}"; do
    input_files=""
    for year in "${years[@]}"; do
        input_files+=" $base_path/ALP_${mass}/${year}.root"
    done

    output_file="$base_path/ALP_${mass}/run3.root"
    
    # Check if output file exists and remove it
    if [ -f "$output_file" ]; then
        echo "File exists: $output_file — removing it."
        rm -f "$output_file"
    fi
    
    echo "hadd $output_file $input_files"
    hadd $output_file $input_files
done

# Add run3.root in all ALP mass points
Sig_MC_path=/eos/home-p/pelai/HZa/Root_Dataset/run3/All_Sig
if [ -d "$Sig_MC_path" ]; then
    echo "Directory exists: $Sig_MC_path — removing it."
    rm -rf "$Sig_MC_path"
else
    echo "Directory does not exist: $Sig_MC_path — creating it."
fi
mkdir -p "$Sig_MC_path"

input_files=""
for mass in "${massList[@]}"; do
    input_files+=" $base_path/ALP_${mass}/run3.root"
done
output_file="$Sig_MC_path/run3.root"

# Check if output file exists and remove it
if [ -f "$output_file" ]; then
    echo "File exists: $output_file — removing it."
    rm -f "$output_file"
fi

echo "hadd $output_file $input_files"
hadd $output_file $input_files