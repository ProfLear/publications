# -*- coding: utf-8 -*-
"""
Created on Fri Apr 25 10:08:09 2025

@author: benle
"""

from plotly.subplots import make_subplots
from collections import Counter
import re

# Raw text input
pre = """
computer Genomics/Oncology
Buzzword, overhyped, undefined
Neural-Network, Limited, Patchwork
information internet 
Impactful, dangerous, potential 
unskilled
Robot, computer, dangerous
Innovative, dangerous, eye-opening
Cheating
cheating, cheating, revolutionary
Robot, computer, neural-network
Technology, advanced, scary
helpful, interesting, powerful
information, replacement, cheating
robot, future, helpful
Tool, research, organization
Information, cheating, unreliable 
Robot, Chat-GPT, futuristic
science, evolving, data
potential, dystopian, encompassing
Processing, Data, Potential
adaptive, helpful, scary
Innovation, Dangerous, Advanced
innovation, practical, futuristic
Challenging cheating 
dangerous, inhumane, scary
robot, science, scary
fast, convenient, information
helpful, fast, evolving
algorithm, data, all-knowing
dangerous, helpful, complementary
evolving, innovative, dangerous
Scary, Fascinating, Mind-blowing
Robot, scary, fast
scary, cheating, helpful
Chat-GPT, Robot, Machine-Learning
fast, unreliable, complex
Automated, computer, Brain
Futuristic, Robot, Matrix
unnatural, helpful, cautious
cheating, dangerous, unknown
helpful tool exciting
Adaptive, manipulative, untrustworthy   
robot, thinking, idea
future, overwhelming, complex
helpful, dangerous, providable
WALLE, chat-gpt, science
helpful, smart, fast
Interesting, evolving, Aid
Data, computer, network
"""

post = '''
Brainstorm Answer Helpful
Chat-GPT, academic integrity, AI detector
Chat-GPT
Robot, generative, informational
helpful tool
Helpful, interesting, evolving 
Summarize, dangerous, Cheap
Chat-GPT, machine-learning, Generative
helpful, limited, exciting 
helpful, powerful, data
helpful, cheating, scary
careful, helpful, time-saver
Assistance generative Tool  
helpful, scary, resource
Data, Analyze, Algorithm
helpful, smart, human-like
dangerous, helpful, evolving
unreliable, robot, helpful
Helpful, valuable, evolving
Robot, Technology, Science
scary, helpful, complex
Helpful, dangerous, unknowing
 Network, innovation, Generation
data,  network, evolving
Robot, tool, nonhuman 
Computer, generative, adaptive
Helpful, intriguing, easy
generative, computer, mind
Learn, adapt, tool
unreliable, evolving, evolving
Computer, assistance, algorithm
dangerous, evolving, innovation
smart, innovation, generative
cheating, time-saver, choice
helpful, helpful, resource
Fast, helpful, research
evolving, innovation, helpful
innovation, evolving, misunderstood
dangerous, helpful, fast
Adaptive, helpful, unreliable
Unreliable, scary, dangerous
Misunderstood, Helpful, Complex
attempt, data, analyze
Chat-GPT, human-like Robot, Machine-Learning
study, learn, data
helpful, unreliable
'''

# Normalize text: lowercase, split words but keep hyphenated and slash-separated words intact
prewords = re.findall(r'\b[\w/-]+\b', pre.lower())
postwords = re.findall(r'\b[\w/-]+\b', post.lower())

# Count word frequency
preword_counts = Counter(prewords)
postword_counts = Counter(postwords)

# Sort counts
# Sort: by descending count, then alphabetically for ties
presorted_counts = sorted(preword_counts.items(), key=lambda x: (-x[1], x[0]))
postsorted_counts = sorted(postword_counts.items(), key=lambda x: (-x[1], x[0]))

# Separate into two lists
prewords_list, precounts = zip(*presorted_counts)
postwords_list, postcounts = zip(*postsorted_counts)


# Create Plotly bar chart
fig = make_subplots(rows = 2, cols = 1)

cuttoff = 8
fig.add_bar(x=precounts[0:cuttoff], y=prewords_list[0:cuttoff], orientation = "h", marker = dict(color = "#008888"), showlegend = False, row = 1, col = 1)
fig.add_bar(x = postcounts[0:cuttoff], y = postwords_list[0:cuttoff], orientation  = "h", marker = dict(color = "#f5653d"), showlegend = False, row = 2, col =1)



fig.update_xaxes(range = [0, max(max(precounts), max(postcounts))])
fig.update_yaxes(autorange="reversed")

# Show figure
fig.update_yaxes(
    linecolor = "#a9a9a9",
    tickcolor = "#a9a9a9",
    tickfont = dict(color = "#a9a9a9"),
    )
fig.update_xaxes(
    linecolor = "#a9a9a9",
    tickcolor = "#a9a9a9",
    tickfont = dict(color = "#a9a9a9"),
    )
for c in [1, 2]:
    fig.update_xaxes(title = "counts", row = 2, col = 1)

fig.update_layout(
    template = "simple_white",
    width = int(3.3*300/2),
    height = int(1.25*300),
    margin = dict(l = 10, t = 10, r = 10, b = 50),
    font = dict(size = 12, color = "#a9a9a9"),
    )
fig.show("png")
fig.write_image("sentiment.svg")
