# -*- coding: utf-8 -*-
"""
Created on Tue May  6 17:46:10 2025

@author: benle
"""


file = "G:/My Drive/PennState/Research/Manuscripts/2025/Morgan JChemEd/CHEM 110H AI Project/Course/CHEM 110H_Post Course Survey_Proccessed.csv"
with open(file, "r") as ps:
    for row in ps: # go row by row
        entries = row.split(",")
        print(len(entries))
        data_to_write = []
        for i in range(0, int(len(entries)/3)):
            data_to_write.append([])
        break
#
with open(file, "r") as ps:
    table_numbers = range(0, len(data_to_write))
    i = 0
    for row in ps: # go row by row
        entries = row.split(",")
        for tn in table_numbers:
            print(tn)
            if i == 0: #we are getting headers
                data_to_write[tn].append(r'''
\begin{table}[]
    \centering
    \caption{''' + f'''Responses for question {tn+1}. The survey question was: ``{entries[tn*3].split(":")[1].strip()}"''' + r'''}
    \begin{tabular}{p{0.33\linewidth} r}
         Response &  Counts\\ \hline
''' + "\n"
                                        )
            else: #we are in the data
                data_to_write[tn].append("\t\t" + f"{entries[tn*3]}" + " & " + f"{entries[tn*3+1]}" + r"\\" + "\n" )
        i = i+1
        
for i, table in enumerate(data_to_write):
    table.append(r'''\hline
    \end{tabular}
    \label{''' + f'''tab:question{i+1}''' + r'''}
\end{table}
'''
                 )
    
    
with open("G:/My Drive/PennState/Research/Manuscripts/2025/Morgan JChemEd/CHEM 110H AI Project/postsurveyresults.tex", "w") as tx:
    for table in data_to_write[1:9]+data_to_write[21:]: #skip the first table
        for line in table:
            if "\t\t & \\" not in line:
                tx.write(line)
    