#!/bin/bash

## with pho
echo $(date)
echo 'Searching for MN...'
python3 ../scripts/1.pdb_search_for_ligands_2.6A_choose_metal.py MN > logs/1.log_search.txt. ### Give the metal of interest as their PDB ID
echo $(date)
echo 'Processing...'
python3 ../scripts/2.ligands_process_01.py >> logs/2.log_ligands.txt 
python3 ../scripts/3.ligands_process_02.py >> logs/3.log_lig_process.txt
echo $(date)
echo 'Getting submols...'
python3 ../scripts/4.extract_groups_2.6A_images_opt.py > logs/4.log_submols.txt ### Image generation is optional. If not needed, type an argument after the .py.
#python3 ../scripts/4.extract_groups_2.6A_images_opt.py no > logs/4.log_submols.txt ### No image will be generated.
echo "Finished\n\n"
echo $(date)

echo $(date)
echo 'Searching for MG...'
cd ../mg
python3 ../scripts/1.pdb_search_for_ligands_2.6A_choose_metal.py MG > logs/1.log_search.txt
echo $(date)
echo 'Processing...'
python3 ../scripts/2.ligands_process_01.py >> logs/2.log_ligands.txt
python3 ../scripts/3.ligands_process_02.py >> logs/3.log_lig_process.txt
echo $(date)
echo 'Getting submols...'
python3 ../scripts/4.extract_groups_2.6A_images_opt.py > logs/4.log_submols.txt
echo "Finished\n\n"
echo $(date)

echo $(date)
echo 'Searching for ZN...'
cd ../zn
python3 ../scripts/1.pdb_search_for_ligands_2.6A_choose_metal.py ZN > logs/1.log_search.txt
echo $(date)
echo 'Processing...'
python3 ../scripts/2.ligands_process_01.py >> logs/2.log_ligands.txt
python3 ../scripts/3.ligands_process_02.py >> logs/3.log_lig_process.txt
echo $(date)
echo 'Getting submols...'
python3 ../scripts/4.extract_groups_2.6A_images_opt.py > logs/4.log_submols.txt
echo "Finished\n\n"
echo $(date)
