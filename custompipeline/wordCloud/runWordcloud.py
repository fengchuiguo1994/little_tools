# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 19:06:59 2019

@author: huang
"""

import jieba
from wordcloud import WordCloud
import re
import time
from string import punctuation
import matplotlib.pyplot as plt

add_punc='，。、【 】 “”：；（）《》‘’{}？！⑦()、%^>℃：.”“^-——=&#@￥'
all_punc=punctuation+add_punc
test = "|".join(['\\'+x for x in list(all_punc)])
r = re.compile(test)

text = ''
t = time.time()
with open(r"C:\Users\huang\Desktop\files.txt",'r',encoding='utf-8') as fin:
    for line in fin:
        line = line.strip()
        line = re.sub(r, '', line)
        text = text+line
print(time.time()-t)
# print(text)

text_list = []
t = time.time()
with open(r"C:\Users\huang\Desktop\files.txt",'r',encoding='utf-8') as fin:
    for line in fin:
        line = line.strip()
        line = re.sub(r, '', line)
        text_list.append(line)
text = ''.join(text_list)
print(time.time()-t)
# print(text)

word = {}
black = set(['的'])
for i in jieba.cut(text):
    if i not in word:
        word[i] = 0
    word[i] += 1
for i in black:
    if i in word:
        del(word[i])
print(word)

wordcloud = WordCloud("test",width=600,height=400)
wordcloud.add("",word.keys(),word.values(),word_size_range=[20,100])
wordcloud.render()