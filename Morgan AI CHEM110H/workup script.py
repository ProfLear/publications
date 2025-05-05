# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 12:01:00 2025

@author: benle
"""
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import chi2_contingency # for testing differences in responses.
import scipy.stats


#Go through each csv file (per exam)
# extract data by question
# add to chart


# this will get all the answers from all surveys.  Then we need to know how to handle them

# get files
files = Path("G:/My Drive/PennState/Research/Manuscripts/2025/Morgan JChemEd/CHEM 110H AI Project").rglob("*.csv")
responses = {}
for file in files:
    responses[f"{file.stem}"] = {}
    # start a plot
    fig = make_subplots(rows = 4, cols = 4)
    #open the file
    #extract information
    read_csv = pd.read_csv(file) # actually, have this read the excel files
    print(file)
    #print(read_csv.head())
    for column in read_csv: # returns keys
        if ":" in column and "Unnamed" not in column:
            number, key = column.split(":", 1)
            number = int(number)
            if number not in responses[f"{file.stem}"]:
                responses[f"{file.stem}"][number] = {}
 
                responses[f"{file.stem}"][number]["question"] = key
                options = read_csv[column][read_csv[column].notna()].copy()
                counts = read_csv[f' {number}: counts'][read_csv[f' {number}: counts'].notna()].copy()
                responses[f"{file.stem}"][number]["answers"] = pd.DataFrame({'options':options, 'counts':counts})

#%% Pre-surveys





#%% Plotting pre-exam surveys
for question in range(2, 10):# get the question number
    n_rows = 1
    n_cols = 2
    fig = make_subplots(rows = n_rows, cols = n_cols,)
    i = 0
    for file in responses: # go through each file
        print(i)
        print(f"row: {int((i)/n_cols)+1}")
        print(f"col: {int(i)%n_cols + 1} ")
        row = 1
        if "pre-exam" in file.lower():
            if "Exam 1" in file:
                col = 1
            if "Exam 2" in file:
                col = 2
            fig.add_bar(y = responses[file][question]['answers']["counts"], x = responses[file][question]['answers']["options"], 
                        row = row, col = col, showlegend = False
                        )
            fig.add_annotation(
                text = f"n = {int(np.sum(responses[file][question]['answers']['counts']))}" ,
                showarrow = False,
                #xref = "paper", yref = "paper",
                x = 1, y = np.max(responses[file][question]['answers']["counts"]),
                row = row, col = col
                )
            i = i+1
    max_y = 0
    for data in fig.data:
        if np.max(data.y) > max_y:
            max_y = np.max(data.y) 
    fig.update_xaxes(range = [0.5, 5.5])
    fig.update_yaxes(range = [0, max_y])
    fig.update_layout(
        title = f"Question {question}: {responses[file][question]['question'].replace(' Please select how much you agree or disagree with the following statements: -','')}",
        width = 3.3*300, height = 2*300,
        template = "simple_white"
        )
    fig.show("png") 
        

#%% Now, let us compare the pre-post COURSE
import numpy as np
import scipy.stats as stats


for question in list(range(29, 50)):# get the question number
    question_map = [
        [],
        ]
    pre_responses = []
    post_responses = []
    n_rows = 1
    n_cols = 1
    fig = make_subplots(rows = n_rows, cols = n_cols,)
    i = 0
    for file in responses: # go through each file
        print(i)
        print(f"row: {int((i)/n_cols)+1}")
        print(f"col: {int(i)%n_cols + 1} ")
        row = 1
        if "course" in file.lower():
            if "pre" in file.lower():
                col = 1
                for number, count in zip(responses[file][question]['answers']["options"], responses[file][question]['answers']["counts"]):
                    for i in range(int(count)):
                        pre_responses.append(number)
            if "post" in file.lower():
                col = 1
                for number, count in zip(responses[file][question]['answers']["options"], responses[file][question]['answers']["counts"]):
                    for i in range(int(count)):
                        post_responses.append(number)
            fig.add_bar(x = responses[file][question]['answers']["counts"], y = responses[file][question]['answers']["options"], 
                        row = row, col = col, showlegend = False, orientation = "h"
                        )
            #fig.add_annotation(
            #    text = f"n = {int(np.sum(responses[file][question]['answers']['counts']))}" ,
            #    showarrow = False,
            #    #xref = "paper", yref = "paper",
            #    x = 1, y = np.max(responses[file][question]['answers']["counts"]),
            #    row = row, col = col
            #    )
            i = i+1
            
            #fig.update_layout(
            #    title = f"Question {question}: {responses[file][question]['question'].replace(' Please select how much you agree or disagree with the following statements: -','')}",
            #    width = 3.3*300, height = 2*300,
            #    template = "simple_white"
            #    )
    max_x = 0
    for data in fig.data:
        if np.max(data.x) > max_x:
            max_x = np.max(data.x) 
    fig.update_yaxes(
        range = [0.5, 5.5],
        tickvals = [1, 2, 3, 4, 5],
        ticktext = ["strongly disagree", "agree", "neutral", "agree", "strongly agree"],
        )
    fig.update_xaxes(
        title = "counts",
        range = [0, 35.5]
        )
    
    
    
    

    stat, p_value = stats.mannwhitneyu(pre_responses, post_responses, alternative='two-sided')
    #chi2, p_value, dof, expected = scipy.stats.chi2_contingency([fig.data[0].x, fig.data[1].x])
    
    fig.update_layout(
        title = f"<b>Question {question}:</b> median1 = {np.median(pre_responses):.1f}, median2 = {np.median(post_responses):.1f},  p = {p_value:.3f}",
        #title = f"<b>Q{question}:</b>{responses['CHEM 110H AI Pre-Course_Proccessed'][question]['question'].split('-')[1]}",
        template = "simple_white",
        width = int(3.3*300/2),
        height = int(1.*300),
        margin = dict(l = 50, t = 50, r = 10, b = 50),
        font = dict(size = 10)
        )
    
    fig.show("png") 
    fig.write_image(f"Question {question}.svg")


#%% significant difference in student familiarity with AI

familiar_fig = make_subplots()
familiar_fig.add_bar(x = [15, 33, 2], orientation  = "h", showlegend = False,)
familiar_fig.add_bar(x = [33, 15, 1], orientation  = "h", showlegend = False,)

labels = ["very familiar", "somewhat familiar", "not sure what it is"]

familiar_fig.update_xaxes(title = "counts", range = [0, 33*1.1])
familiar_fig.update_yaxes(range = [2.5, -0.5], tickvals = [0, 1, 2], ticktext = ["very familiar", "somewhat familiar", "not sure what it is"])

chi2, p_value, dof, expected = scipy.stats.chi2_contingency([[15, 33, 2], [33, 15, 1]])


familiar_fig.update_layout(
    title = f"<b>Question 2:  p = {p_value:.3f}",
    template = "simple_white",
    )
familiar_fig.show("png")
#%% no signifcant different in how people percieve the description of AI

describe_fig = make_subplots()
describe_fig.add_bar(x = [40, 10], orientation = "h", showlegend = False)
describe_fig.add_bar(x = [37, 9], orientation = "h", showlegend = False)
"""
labels["
   	options
0	AI is a system that can learn and adapt based on data
	options
1	mimicking certain aspects of human intelligence

       ",
       "
   	options
2	AI is a set of algorithms and technologies designed to perform tasks that typically require human intelligence
	options
3	like understanding language or recognizing images

       "]
"""

chi2, p_value, dof, expected = scipy.stats.chi2_contingency([[40, 10], [37, 9]])


describe_fig.update_layout(
    title = f"<b>Question 2:  p = {p_value:.3f}",
    template = "simple_white",
    )
describe_fig.show("png")
#%% no signficant change in how student percieve the components of AI
components_fig = make_subplots()
components_fig.add_bar(x = [45, 33, 26, 24, 4])
components_fig.add_bar(x = [42, 36, 32, 20, 2])

labels = ["machine learning", "big data", "neural networks", "robotics", "not sure"]

chi2, p_value, dof, expected = scipy.stats.chi2_contingency([[45, 33, 26, 24, 4], [42, 36, 32, 20, 2]])


components_fig.update_layout(
    title = f"<b>Question 2:  p = {p_value:.3f}",
    template = "simple_white",
    )
components_fig.show("png")

#%% How students use...
use_fig = make_subplots()
pre_labels = []
pre_use = []
post_use = []
pre_index = []
post_index = []
for i, option in enumerate(responses["CHEM 110H AI Pre-Course_Proccessed"][25]["answers"]["options"]):
    pre_labels.append(option)
    pre_index.append(i)
    pre_use.append(responses["CHEM 110H AI Pre-Course_Proccessed"][25]["answers"]["counts"][i])
    for j, poption in enumerate(responses["CHEM 110H_Post Course Survey_Proccessed"][25]["answers"]["options"]):
        if option == poption:
            post_use.append(responses["CHEM 110H_Post Course Survey_Proccessed"][25]["answers"]["counts"][j])
            post_index.append(j)
use_fig.add_bar(x = pre_use, orientation = "h", showlegend = False, marker = dict(color = "#008888"),)
use_fig.add_bar(x = post_use, orientation = "h", showlegend = False, marker = dict(color = "#f5653d"),)


use_fig.update_xaxes(title = "counts")
use_fig.update_yaxes(range = [len(pre_use)-0.5, -0.5], tickvals = list(range(len(pre_labels))), ticktext = pre_labels)

chi2, p_value, dof, expected = scipy.stats.chi2_contingency([pre_use, post_use])

use_fig.update_yaxes(
    linecolor = "#a9a9a9",
    tickcolor = "#a9a9a9",
    tickfont = dict(color = "#a9a9a9"),
    )
use_fig.update_xaxes(
    linecolor = "#a9a9a9",
    tickcolor = "#a9a9a9",
    tickfont = dict(color = "#a9a9a9"),
    )
for c in [1, 2]:
    fig.update_xaxes(title = "counts", row = 2, col = c)

use_fig.update_layout(
    template = "simple_white",
    width = int(3.3*300/2),
    height = int(1.*300),
    margin = dict(l = 0, t = 10, r = 10, b = 50),
    font = dict(size = 12, color = "#a9a9a9"),
    )
#use_fig.add_annotation(text = f"n = {np.array(use_fig.data[0].x).sum()}", showarrow = False, font = dict(color = "#008888"), x = 10, y = 8)
#use_fig.add_annotation(text = f"n = {np.array(use_fig.data[1].x).sum()}", showarrow = False, font = dict(color = "#008888"), x = 10, y = 9)
use_fig.show("png")
use_fig.write_image("HowUsed.svg")

#%%

use_comp_plot = make_subplots(cols = 3, rows = 1, horizontal_spacing = 0)
use_comp_plot.add_bar(x = responses["CHEM 110H AI Pre-Course_Proccessed"][25]["answers"]["counts"], 
                orientation = "h", 
                showlegend = False, 
                col = 1, row = 1)
use_comp_plot.add_bar(x = responses["CHEM 110H_Post Course Survey_Proccessed"][25]["answers"]["counts"], 
                orientation = "h", 
                showlegend = False, 
                col = 3, row = 1)

#add lines
for pre_i, post_i in zip(pre_index, post_index):
    use_comp_plot.add_scatter(x = [0, 1], y = [pre_i, post_i],
                              mode="lines",
                              line = dict(color = "grey"),
                              showlegend = False,
                         row = 1, col = 2)

use_comp_plot.update_xaxes(range=[max(responses["CHEM 110H_Post Course Survey_Proccessed"][25]["answers"]["counts"])*1.1, 0],
                           row = 1, col = 1)
use_comp_plot.update_xaxes(showticklabels = False, ticks= "", linecolor = "white", row = 1, col = 2)
use_comp_plot.update_xaxes(range=[0, max(responses["CHEM 110H_Post Course Survey_Proccessed"][25]["answers"]["counts"])*1.1],
                           row = 1, col = 3)

use_comp_plot.update_yaxes(range = [len(responses["CHEM 110H_Post Course Survey_Proccessed"][25]["answers"]["counts"])-1.5, -0.5],
                           ticks = "", linecolor = "white", showticklabels = False)

for i, option in enumerate(responses["CHEM 110H AI Pre-Course_Proccessed"][25]["answers"]["options"]):
    use_comp_plot.add_annotation(text = option,
                                 x = responses["CHEM 110H AI Pre-Course_Proccessed"][25]["answers"]["counts"][i],
                                 y = i,
                                 showarrow = False,
                                 xanchor = "right",
                                 yanchor = "middle")
                                  
use_comp_plot.update_layout(template = "simple_white",
                            margin = dict(l = 200))

use_comp_plot.show("browser")

#%%
use_fig.add_bar(x = [45, 33, 26, 24, 4])
use_fig.add_bar(x = [42, 36, 32, 20, 2])

labels = ["machine learning", "big data", "neural networks", "robotics", "not sure"]

chi2, p_value, dof, expected = scipy.stats.chi2_contingency([[45, 33, 26, 24, 4], [42, 36, 32, 20, 2]])


components_fig.update_layout(
    title = f"<b>Question 2:  p = {p_value:.3f}",
    template = "simple_white",
    )
components_fig.show("png")

#%% Sentiment analysis
pre_comments = ""
for file in responses:
    if "Pre-Course" in file:
        for answer, count in zip(responses[file][53]["answers"]["options"], responses[file][53]["answers"]["counts"]):
            for i in range(int(count)):
                pre_comments = pre_comments + answer + " "

post_CHEM110H_thoughts = ""
for file in responses:
    if "Post Course" in file:
        #print(responses[file][68]["answers"]["counts"])
        for answer, count in zip(responses[file][66]["answers"]["options"], responses[file][66]["answers"]["counts"]):
            for i in range(int(count)):
                post_CHEM110H_thoughts = post_CHEM110H_thoughts + answer + " "


post_comments = ""
for file in responses:
    if "Post Course" in file:
        #print(responses[file][68]["answers"]["counts"])
        for answer, count in zip(responses[file][68]["answers"]["options"], responses[file][68]["answers"]["counts"]):
            for i in range(int(count)):
                post_comments = post_comments + answer + " "

import nltk
from nltk.sentiment import SentimentIntensityAnalyzer

# Download VADER lexicon (only once)
nltk.download('vader_lexicon')

# Initialize Sentiment Analyzer
sia = SentimentIntensityAnalyzer()

# Run Sentiment Analysis on the entire aggregated text
aggregated_pre_sentiment = sia.polarity_scores(pre_comments)

appregated_CHEM110H_sentiment = sia.polarity_scores(post_CHEM110H_thoughts)
aggregated_post_sentiment = sia.polarity_scores(post_comments)

# Print Results
print("Aggregated Sentiment Pre-Scores:", aggregated_pre_sentiment)
print("Aggregated Sentiment CHEM110H:", appregated_CHEM110H_sentiment)
print("Aggregated Sentiment Post-Scores:", aggregated_post_sentiment)

# Determine Overall Sentiment
compound_score = aggregated_pre_sentiment['compound']
overall_sentiment = "Positive" if compound_score > 0.05 else "Negative" if compound_score < -0.05 else "Neutral"

#print(f"Overall Sentiment: {overall_sentiment}")



        
#%% Can be adapted to any number of files

for question in range(2, 12):# get the question number

    fig = make_subplots(rows = n_rows, cols = n_cols,)
    i = 0
    for i, file in enumerate(responses): # go through each file
        print(i)
        print(f"row: {int((i)/n_cols)+1}")
        print(f"col: {int(i)%n_cols + 1} ")

        fig.add_bar(y = responses[file][question]['answers']["counts"], x = responses[file][question]['answers']["options"], 
                    row = int((i)/n_cols)+1, col = int(i-2)%n_cols + 1,
                    )
        fig.add_annotation(
            text = file ,
            showarrow = False,
            xref = "paper", yref = "paper",
            x = 0, y = 0,
            row = int(int(i-2)/n_cols)+1, col = int(i)%n_cols + 1
            )
    
    fig.update_xaxes(range = [0.5, 5.5])
    fig.update_layout(
        width = 9*300, height = 6.5*300,
        template = "simple_white"
        )
    fig.show("png") 
    
    
    
#%% Prepare figures for publication
indepdence = [31,32,33,36]
concerns = [37, 39, 40, 41]
utility = [45, 46, 47, 49]
fig = make_subplots(rows = n_rows, cols = n_cols,
                    shared_xaxes = True, shared_yaxes = True)
for i, question in enumerate(utility):# get the question number
    question_map = [
        [],
        ]
    pre_responses = []
    post_responses = []
    n_rows = 2
    n_cols = 2

    for file in responses: # go through each file
        print(i)
        print(f"row: {int((i)/n_cols)+1}")
        print(f"col: {int(i)%n_cols + 1} ")
        row = int((i)/n_cols)+1
        if "course" in file.lower():
            if "pre" in file.lower():
                col = int(i)%n_cols + 1
                for number, count in zip(responses[file][question]['answers']["options"], responses[file][question]['answers']["counts"]):
                    for j in range(int(count)):
                        pre_responses.append(number)
                fig.add_bar(x = responses[file][question]['answers']["counts"], y = responses[file][question]['answers']["options"], 
                            marker = dict(color = "#008888"),
                            row = row, col = col, showlegend = False, orientation = "h"
                            )
            if "post" in file.lower():
                col = int(i)%n_cols + 1
                for number, count in zip(responses[file][question]['answers']["options"], responses[file][question]['answers']["counts"]):
                    for j in range(int(count)):
                        post_responses.append(number)
                fig.add_bar(x = responses[file][question]['answers']["counts"], y = responses[file][question]['answers']["options"], 
                            marker = dict(color = "#f5653d"),
                            row = row, col = col, showlegend = False, orientation = "h"
                            )
            
            
    max_x = 0
    for data in fig.data:
        if np.max(data.x) > max_x:
            max_x = np.max(data.x) 
    fig.update_yaxes(
        range = [0.5, 5.5],
        tickvals = [1, 2, 3, 4, 5],
        ticktext = ["strongly disagree", "agree", "neutral", "agree", "strongly agree"],
        linecolor = "#a9a9a9",
        tickcolor = "#a9a9a9",
        tickfont = dict(color = "#a9a9a9"),
        )
    fig.update_xaxes(
        range = [0, 35.5],
        linecolor = "#a9a9a9",
        tickcolor = "#a9a9a9",
        tickfont = dict(color = "#a9a9a9"),
        )
    for c in [1, 2]:
        fig.update_xaxes(title = "counts", row = 2, col = c)
    
    fig.update_layout(
        template = "simple_white",
        width = int(3.3*300/2),
        height = int(1.*300),
        margin = dict(l = 0, t = 10, r = 10, b = 50),
        font = dict(size = 12, color = "#a9a9a9"),
        )
    
fig.show("png") 
fig.write_image(f"{utility=}.svg")
