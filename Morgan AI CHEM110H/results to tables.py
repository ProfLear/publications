# if I keep working on this, then need to identify those tables that are long, and change this to handle longtable usage for them

def write_surveys_to_tables(file_to_read, file_to_write):
    # first set up the dictionary to hold the responses. 
    with open(file_to_read, "r") as ps:
        for row in ps: # go row by row
            entries = row.split(",")
            print(len(entries))
            data_to_write = {}
            for i in range(1, int(len(entries)/3)+1):
                data_to_write[f"{i}"] = {"question":"", "responses":{}} #append a dictionary for each question that there is
            break
        
    #% now, let us make sure we have all the answers and headers recorded.
    with open(file_to_read, "r") as ps:
        table_numbers = range(0, len(data_to_write))
        i = 0
        for row in ps: # go row by row
            entries = row.split(",")
            for tn in table_numbers:
                print(tn)
                if i == 0: #we are getting headers
                    data_to_write[f"{tn+1}"]["question"] = f"{entries[tn*3].split(":")[1].strip()}"
                else: #we are in the data
                    if f"{entries[tn*3]}".strip() != "": # don't need empty strings
                        data_to_write[f"{tn+1}"]["responses"][f"{entries[tn*3]}"] = f"{entries[tn*3+1]}"
            i = i+1
    
    # now we can write the file, making sure to sort as needed. 
    to_write = list(range(2, 7)) + [8, 9] + list(range(22, 26,)) + list(range(27, len(data_to_write) + 1))
    with open(file_to_write, "w") as tx:   
        for key in data_to_write:
            if int(key) in to_write: # now we have a data table to write
                tx.write(r'''
    \begin{table}[]
        \centering
        \caption{''' + f'''Responses for survey question {key}. The survey question was: ``{data_to_write[key]['question']}"''' + r'''}
        \begin{tabular}{p{0.33\linewidth} r}
            Response &  Counts\\ \hline''' + "\n"       
                        )
                if "3" in data_to_write[key]["responses"]: # this is a likert scale response. 
                    # first, make sure all answer choices are present, even if they don't have answers...
                    for c in range(1, 6):
                        if f"{c}" not in data_to_write[key]["responses"].keys():
                            data_to_write[key]["responses"][f"{c}"] = "0"
                    response_keys = sorted(data_to_write[key]["responses"].keys(), reverse = True)
                else:
                    response_keys = data_to_write[key]["responses"].keys()
                
                running_total = 0
                for r in response_keys: 
                    key_string = r.replace("&", "\&")
                    response_string = data_to_write[key]['responses'][r]
                    response_string.replace("&", r"\&")
                    tx.write("\t\t" + f"{key_string}" + " & " + f"{response_string}" + r"\\" + "\n" )
                    running_total= running_total +  int(data_to_write[key]['responses'][r])
                
                tx.write("\t\t" + r'\hline ' + "\n")
                tx.write("\t\t" + r"\textbf{total counts} & \textbf{" + f"{running_total}" + r"}\\" + "\n")
                
                tx.write("\t" + r'''\end{tabular}
        \label{''' + f'''tab:question{key}''' + r'''}
    \end{table}
            '''
                             )
    
for r, w in zip([
        r"G:\My Drive\PennState\Research\Manuscripts\2025\Morgan JChemEd\CHEM 110H AI Project\Course\CHEM 110H AI Pre-Course_Proccessed.csv",
        "G:/My Drive/PennState/Research/Manuscripts/2025/Morgan JChemEd/CHEM 110H AI Project/Course/CHEM 110H_Post Course Survey_Proccessed.csv",
        ], 
        [
            r"G:\My Drive\PennState\Research\Manuscripts\2025\Morgan JChemEd\CHEM 110H AI Project\presurveyresults.tex",
            "G:/My Drive/PennState/Research/Manuscripts/2025/Morgan JChemEd/CHEM 110H AI Project/postsurveyresults.tex",
            ]):
    print(r)
    write_surveys_to_tables(r, w)
    