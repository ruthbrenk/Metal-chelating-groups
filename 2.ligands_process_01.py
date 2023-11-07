# Gather ligands from the atom lines (unwanted compounds already removed)
# Last edited on 05/09/23

def search_string_in_file(file_name, string_to_search):
	"""Search for the given string in file and return lines containing that string"""
	line_number = 0
	list_of_results = []
	
	with open(file_name, "r") as read_obj: # Open the file in read only mode
		for line in read_obj: # Read all lines in the file one by one
			line_number += 1 # For each line, check if line contains the string
			if string_to_search in line:
				words = line.split(sep=" ")[2][:-5]
				code = words.rsplit('_')
				het = code[1]
				pdb = code[0]
				if '## Created' in line:
					list_of_results.append((line_number, pdb, het)) # If yes, then add the line number & line as a tuple in the list
				if '## Deleted' in line:
					rm_line = line_number - 2
					list_of_results.remove((rm_line, pdb, het)) # If yes, removes the respective tuple in the list

	return list_of_results # Return list of tuples containing line numbers and lines where string is found

matched_lines = search_string_in_file("logs/1.log_search.txt", ".sdf") # Specify file & string

for elem in matched_lines:
	print("PDB ID = " + elem[1] + " :: " + elem[2])