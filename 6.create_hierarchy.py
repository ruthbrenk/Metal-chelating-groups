import pandas as pd
import sys, copy,string

out_file = open('out.txt', 'w')

df = pd.read_csv('smarts_substructure_match.csv', index_col=0)

 #create dictionary
global family_dict 
family_dict = {}   


def make_subtable(data_frame, row_name,fam_number,head_number_dict):
	col_names = []
	new_df = data_frame.copy(deep=True)

	drop_list = []

	for col in data_frame.columns:
            if col != 'total_x': #we don't need to check for this column
               col_names.append(col)
	for col_name in col_names:
		value =data_frame.loc[row_name][col_name]
		if int(value) ==0:
	    		new_df = new_df.drop(columns = [col_name])
	    		drop_list.append(col_name)
		elif int(value) ==1:
	    		#record familiy membership
                        if fam_number in head_number_dict:
                            head_number_dict[fam_number].append(col_name)
                        else:
                            head_number_dict[fam_number] = [col_name]



                        if col_name in family_dict:
                            family_dict[col_name].append([fam_number,row_name])
                        else:
                            family_dict[col_name] = [[fam_number,row_name]]
	    		
	number = data_frame.loc[row_name]['total_x']

	
	#delete all row for which we have no more a columns
	for name in drop_list:
		#print (name)
		new_df = new_df.drop(name)	
		
	new_df = new_df.drop('total_y')
	
	new_df = new_df.drop(columns = ['total_x'])
	new_df['total_x'] = new_df.sum(axis=1)
	
	new_df.loc['total_y']= new_df.sum()
	
	new_df = new_df.sort_values(by = 'total_x',ascending=False) 
	    		
	print ('*************')
	level = fam_number.count('.')
	level2 = str('{:02d}'.format(level))
	print (fam_number, f':LEVEL{level2}:', row_name, number)
	print ('*************')
	
	return new_df

#------------------------------------------------
def check_if_head(value, row_name,df,group_number):
	is_head = False
	row_name = row.name
	value =df.loc['total_y'][row_name]	   
	if value == 0: #then it is the head and not part of any other group
		is_head=True
	elif  group_number.count('.') > 0: #check if it only occurs in heads a level above, if is the case, it is still head, no need to check in 1. level
		group_number_level = group_number.count('.')
		col_names = []
		for col in df.columns:
			if col != 'total_x': #we don't need to check for this column
				col_names.append(col)
		for col_name in col_names:
			hit =df.loc[row_name][col_name]
			if int(hit) ==1: 
				#check if row_name and col_name are the same molecule, only written diffrentl, e.g. if they have the same number of atoms
				if col_name in head_number_dict:
					#print ('---------> group head <--------------------')
					head_number = head_number_dict[col_name]
					head_number_level = head_number.count('.')
					#print (head_number_level, group_number_level)
					if head_number_level  == group_number_level: #occurs in head one level above
						value = value-1
						print (value, 'value subsructured', row_name, col_name)

		if value == 0:
			is_head = True
	
	return is_head
						
					

#add total row
df.loc['total_x']= df.sum()
    
df['total_y'] = df.sum(axis=1)

#print (df.head)
#print (df)

df =df.transpose()

#print (df.head)

df = df.sort_values(by = 'total_x',ascending=False) 

#print (df[0:3])



#counter = 0
head_dict = {}
head_number_dict = {}
new_df_list = []
old_df_list = [['0',df]] #[[group_number,data_frame]]
while len(old_df_list) > 0: #there are still dfs to process
    for group_number, df in old_df_list:
        counter = 0
        # loop through the rows using iterrows()
        for index, row in df.iterrows():
            if index != "total_y":
                row_name = row.name
                value =df.loc['total_y'][row_name]	
                if check_if_head(value, row_name,df,group_number):   		
                #if value == 0: #then it is the head and not part of any other group
                    #we need to save the head in the head dict
                    counter=counter + 1
                    fam_number  = group_number + '.' + str(counter)
                    head_dict[fam_number] = row_name

                    value =df.loc[row_name]['total_x']	    		
                    if value > 0: #if it 0 we are done		    		
	    
                                #print ('make subtable')
                                sub_frame = make_subtable(df, row_name,fam_number,head_number_dict)
                                if len(sub_frame) > 1: #we need more than just the total row
                                    new_df_list.append([fam_number,sub_frame])

    old_df_list = copy.deepcopy(new_df_list) #need to make a deep copy, otherwise the list are still linked
    new_df_list = []



for key in dict(sorted(head_dict.items())):
    #out_file.write(key + '\t' + head_dict[key] + '\n')
    if key in head_number_dict:
        out_file.write(key + '\t' + head_dict[key] + '\n') #always write the key and the smarts of the familiy if there are familiy members
        out_file.write(str(len(head_number_dict[key])) +'\n') #print out how many members are in this level
        for smarts in head_number_dict[key]:
            out_file.write(smarts + '\t')
        out_file.write('\n')
    elif key.count('.') == 1:
        out_file.write(key + '\t' + head_dict[key] + '\n') #always write if it is a 0 level smarts

print(head_dict)
