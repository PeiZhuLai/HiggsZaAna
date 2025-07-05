import re
import os

# Define the regex pattern to extract the value
pattern = re.compile(r'- Best category:  -0.1 (\S+) --->')

# Define the list of components to extract the value for
components_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]

# Open the output file to store all extracted values
output_file_name_BDTcut = "BDT_cut_all_ma.txt"
output_file_BDTcut = open(output_file_name_BDTcut, 'w')

output_file_name_significance = "significance_all_ma.txt"
output_file_significance = open(output_file_name_significance, 'w')

# Loop through each component in the list and extract the value
for component in components_list:
    input_file_name = f"/afs/cern.ch/work/p/pelai/HZa/ALP/Plot/optimize_run2UL/categorize_all_M{component}.txt"
    with open(input_file_name, 'r') as file:
        list = open(input_file_name)
        for line in list:
            line.rstrip()
            line.lstrip()
            
            BDT_value = float(line.split()[4])

            if component == 1:
                output_file_BDTcut.write("mvaCuts = {")
                output_file_BDTcut.write(f"{component}:{round(BDT_value, 3)}, ")
            
            if component > 1 and component < 30:
                output_file_BDTcut.write(f"{component}:{round(BDT_value, 3)}, ")

            if component == 30:
                output_file_BDTcut.write(f"{component}:{round(BDT_value, 3)}")
                output_file_BDTcut.write("}")
            
           
            significance = float(line.split()[7])
            output_file_significance.write(f"Significance M{component}: {round(significance, 3)}\n")

# Close the output file
output_file_BDTcut.close()

output_file_significance.close()
