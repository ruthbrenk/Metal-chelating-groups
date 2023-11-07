# Last edited on 22/09/23
# BSA
# Doesn't cat coord files

import os
from glob import glob
from datetime import datetime

# datetime object containing current date and time
begin = datetime.now()

# dd/mm/YY H:M:S
dt1_string = begin.strftime("%d/%m/%Y %H:%M:%S")
timestamp = begin.strftime('%Y%m%d')
print("Timestamp:", dt1_string, "\n-------------------------------\n")

cwd = os.getcwd()
mf = f'{cwd}/near' # Mother folder

lines_seen = set() # holds lines already seen
outfile = open("logs/2.log_ligands_tmp.txt", "w")
for line in open("logs/2.log_ligands.txt", "r"):
	if line not in lines_seen: # not a duplicate
		outfile.write(line)
		lines_seen.add(line)
outfile.close()

os.system(f"rm logs/2.log_ligands.txt")
os.system(f"sort logs/2.log_ligands_tmp.txt -o logs/2.log_ligands_sorted.txt")
os.system(f"rm logs/2.log_ligands_tmp.txt")

# Count the ligands
text = open("logs/2.log_ligands_sorted.txt", "r")
  
d = {}

print("Each filtered molecule appears in X structures (only for those >1)\n")

for line in text:
    line = line.split(sep=' ')[-1][:-1] # Remove the leading spaces and newline character
    line2 = "%s@" % line

# Create dictionary of ligands
    words1 = line2.split("@") # Split the line into words
    words = list(filter(None, words1))
                        
    for word in words:
        if word in d:
            d[word] = d[word] + 1 # Increment count of word by 1
        else:
            d[word] = 1 # Add the word to dictionary with count 1

i = 2 # initializing range 

ds = dict(sorted(d.items(), key=lambda x: x[1], reverse=True)) # Sort the number of occurrences from the highest to lowest

for key in list(ds.keys()):
        if ds[key] >= i:
            print(key, ":", ds[key]) # Print the contents of ligand dictionary

# Create a file with a list of all ligands ids (compatible to https://www.rcsb.org/downloads/ligands)
with open("logs/ligands_final_list.txt", "w") as newfile:
	for key in sorted(list(d.keys())):
		newfile.write("%s," % (key)) # %s, = as strings separated by comma

newfile.close()

## Remove the last comma from the ligand list file
with open("logs/ligands_final_list.txt", "rb+") as filehandle:
    filehandle.seek(-1, os.SEEK_END)
    filehandle.truncate()
filehandle.close()

# Count ligands in the list file
file = open("logs/ligands_final_list.txt", "r")
read_data = file.read()
search = (",")
word_count = read_data.count(search)
word_count1 = word_count +1
print("\nThere are", word_count1, "ligands\n")

## Download the cifs to get the SMILES 
cifs = []
for line in open("logs/ligands_final_list.txt", "r"):
    ligs = line.split(sep=',')
    os.makedirs(f'{mf}/cifs', exist_ok=True)
    for lig in ligs:
        if os.path.isfile(f'{mf}/cifs/{lig}.cif'):
            continue
        else:
            os.system(f"wget https://files.rcsb.org/ligands/download/{lig}.cif -O {mf}/cifs/{lig}.cif")
            cifs.append(lig)
            #print('Downloaded', str(lig), 'cif file')

print(len(cifs), 'cifs were downloaded')

if len(ligs) - len(cifs) > 0:
    print(len(ligs) - len(cifs), 'were not downloaded:')
    not_downloaded = ligs
    for i in cifs:
        not_downloaded.pop(i)
    print(not_downloaded)

print('\nEND')

end = datetime.now()
dt2_string = end.strftime("%d/%m/%Y %H:%M:%S")
print("\n-------------------------------\nTimestamp:", dt2_string)
