import sys
import os
import re
import unicodedata
import numpy
import operator
import math
import json
import networkx
from stemming.porter import stem
from gensim import corpora,models,similarities

def split(paragraph):
    sentenceEnders = re.compile('[, .!?\n]*')
    sentenceList = sentenceEnders.split(paragraph)
    return sentenceList

def sort_nicely( l ): 
  """ Sort the given list in the way that humans expect. 
  """ 
  convert = lambda text: int(text) if text.isdigit() else text 
  alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
  l.sort( key=alphanum_key )

def gram_1(array1):
    array=[]
    for i in range(len(array1)):
        array.append([array1[i],0])

    #print array
    for i in (range(len(array))):
        for j in (range(len(array))):
            if (array[j][0]==array[i][0]) and (array[i][1] is not -1):
                array[i][1]+=1
                if i is not j :
                    array[j][1]=-1
    saturatedarray=[]
    for i in (range(len(array))):
        if array[i][1] is not -1:
            saturatedarray.append(array[i])
    #print saturatedarray
    b=sorted(saturatedarray, key=operator.itemgetter(1))
    return b

def heuristic(node1,node2):
    return distances[node1-1][node2-1]


#readdata

topic=os.listdir('C:/Users/prakash chand gupta/Desktop/AI PROJECT/Articles')
sort_nicely(topic)

'''
from readdata import *
file=open('C:/Users/harshit/Desktop/project_summer/data/a1d5.txt','w')
file.write(read())
file.close()
'''

#avgsentencelength
'''
from avglength import *
a=[]
for i in range(1,6):
    filename='C:/Users/harshit/Desktop/project_summer/data/a1t'+str(i)+'.txt'
    p1=sent_length(filename)
    p2=word_length(filename,p1)
    #print p1
    #print p2
    a.append([1,float("%.2f" % round(p1,2)),float("%.2f" % round(p2,2))])
for i in range(1,6):
    filename='C:/Users/harshit/Desktop/project_summer/data/a1d'+str(i)+'.txt'
    p1=sent_length(filename)
    p2=word_length(filename,p1)
    #print p1
    #print p2
    a.append([2,float("%.2f" % round(p1,2)),float("%.2f" % round(p2,2))])

print a
'''
source_no = 123
target_no = 111
#factor_chosen = 0.3
file=open('C:/Users/prakash chand gupta/Desktop/AI PROJECT/1000_common.txt','r')
p=file.read().replace("'",' ').replace("`",'').replace(";",'').replace("!",'.').lower()
commonwords=split(p)

#print commonwords
'''

file=open('C:/Users/prakash chand gupta/Desktop/AI PROJECT/Articles2/'+str(243)+'.txt','r')
source_p=file.read().replace("'",' ').replace("`",'').replace(";",'').replace("!",'.').replace('"',' ').lower()

print source_p
temp=split(source_p)
print "\nAFTER SPLITTING ARTICLE INTO WORDS\n"
print temp
print len(temp)
for i in range(len(temp)):
	temp[i]=stem(temp[i])
#	wor=list(set(words))
print "\nAFTER STEMMING\n"
print temp
source = [x for x in temp if x not in commonwords]
print "\nAFTER REMOVING STOPWORDS\n"
print source
print len(source)
#source = list(set(temp1))
source = gram_1(source)
print "\nAFTER COUNTING 1 gram frequency\n"
print source
length_source = sum(source[i][1] for i in range(len(source)))
print length_source
'''
'''
file=open('C:/Users/prakash chand gupta/Desktop/AI PROJECT/Articles2/'+str(target_no)+'.txt','r')
target_p=file.read().replace("'",' ').replace("`",'').replace(";",'').replace("!",'.').lower()
temp=split(target_p)

for i in range(len(temp)):
	temp[i]=stem(temp[i])
#	wor=list(set(words))

target = [x for x in temp if x not in commonwords]
target = gram_1(target)
length_target = sum(target[i][1] for i in range(len(target)))
'''

sentences =[]
for i in range(1,476):
    file=open('C:/Users/prakash chand gupta/Desktop/AI PROJECT/Articles2/'+str(i)+'.txt','r')
    p=file.read().replace("'",' ').replace("`",'').replace(";",'').replace(":",'').replace("!",'.').replace('"',"").replace('\xe2',"").replace('\x80',"").replace('\x9d',"").replace('\x91',"").replace('\x92',"").replace('\x93',"").replace('\x94',"").replace('\x95',"").replace('\x96',"").replace('\x97',"").replace('\x98',"").replace('\x99',"").lower()
    p=unicodedata.normalize('NFKD', unicode(p,errors='ignore')).encode('ascii','ignore')
    temp = split(p)
    temp2 = [x for x in temp if x not in commonwords]
    for i in range(len(temp2)):
        temp2[i]=stem(temp2[i])
    temp3 = [x for x in temp2 if x not in commonwords]
    temp3 = gram_1(temp3)
    sentences.append(temp3)

'''
source_distance=[]
target_distance=[]
'''
'''
for i in range(len(sentences)):
    length_temp = sum(sentences[i][j][1] for j in range(len(sentences[i])))
    count=0
    if len(sentences[i]) < len(source):
        for k in range(len(sentences[i])):
                   for j in range(len(source)):
                               if sentences[i][k][0] == source[j][0] and sentences[i][k][0]!='':
                                   count=count+1
    else:
        for k in range(len(source)):
                   for j in range(len(sentences[i])):
                               if source[k][0] == sentences[i][j][0] and source[k][0]!='':
                                   count=count+1
    source_distance.append((float(count)/len(sentences[i]))/len(source))
#print source_distance
'''
'''
distances=[]

for i in range(len(sentences)):
    distances.append([])
    length_tempi = sum(sentences[i][j][1] for j in range(len(sentences[i])))
    count=0
    
    for j in range(len(sentences)):
        count=0
        length_tempj = sum(sentences[j][m][1] for m in range(len(sentences[j])))
        
        if len(sentences[i]) < len(sentences[j]):
            for k in range(len(sentences[i])):
                   for l in range(len(sentences[j])):
                               if sentences[i][k][0] == sentences[j][l][0] and sentences[i][k][0]!='':
                                   count = round(count + round(math.sqrt(float(sentences[i][k][1]*sentences[j][l][1])/(length_tempi * length_tempj)),5),5)
        else:
            for k in range(len(sentences[j])):
                   for l in range(len(sentences[i])):
                               if sentences[j][k][0] == sentences[i][l][0] and sentences[j][k][0]!='':
                                   count = round(count + round(math.sqrt(float(sentences[j][k][1]*sentences[i][l][1])/(length_tempi * length_tempj)),5),5)
        if count == 0.0:
            distances[i].append(0.0)
        else:
            distances[i].append(round(-math.log(count),2))
#print distances
#print target_distance

f=open("distances.txt","w")
json.dump(distances,f)
f.close()
'''


f=open("distances.txt","r")
distances=json.load(f)
f.close()


f=open("corpus.txt","r")
corpus=json.load(f)
f.close()



tfidf = models.TfidfModel(corpus)
index = similarities.SparseMatrixSimilarity(tfidf[corpus], num_features=9343)

new =[]

for i in range(len(corpus)):
    sims = index[tfidf[corpus[i]]]
    new.append(list(enumerate(sims)))

distances2 = []
'''
for i in range(len(new)):
    distances2.append([])
    for j in range(len(new[i])):
        if(new[i][j][1] == 0.0):
            distances2[i].append(100000000)
        else:
            distances2[i].append(1/new[i][j][1])


f=open("distances2.txt","w")
json.dump(distances2,f)
f.close()
'''


f=open("distances2.txt","r")
distances2=json.load(f)
f.close()




'''
for i in range(len(distances)):
    distances[i][target_no]=8
distances[source_no][target_no]=50
'''

a_path=[]
distances[source_no-1][target_no-1]=50
A=numpy.array(distances)
G=networkx.from_numpy_matrix(A,create_using=networkx.DiGraph())
#print (list(networkx.all_simple_paths(G,source_no,target_no,3)))
#print networkx.dijkstra_path_length(G,3,45)
a_path = networkx.dijkstra_path(G,source_no-1,target_no-1)
index_link=[]
#print a_path
index_link.append(a_path[0])

for i in range(len(a_path)-1):
    temp=[]
    distances[a_path[i]][a_path[i+1]]=50
    distances[a_path[i+1]][a_path[i]]=50
    V=numpy.array(distances)
    G2=networkx.from_numpy_matrix(V,create_using=networkx.DiGraph())
    temp = networkx.dijkstra_path(G2,a_path[i],a_path[i+1])
    for j in range(len(temp)):
        if temp[j]!=a_path[i]:
            index_link.append(temp[j])
            l=len(index_link)
            distances[index_link[l-1]][index_link[l-2]]=50
print "\nLINKING ARTICLES INDEX\n"
for i in range(len(index_link)):
    print index_link[i]+1
  
linkart = [];
'''
A=numpy.array(distances2)
G=networkx.from_numpy_matrix(A,create_using=networkx.DiGraph())
#print (list(networkx.all_simple_paths(G,source_no,target_no,3)))
#print networkx.dijkstra_path_length(G,3,45)
a_path = networkx.dijkstra_path(G,source_no-1,target_no-1)

print "\nLINKING ARTICLES INDEX\n"
print a_path

index_link2 = a_path
'''

distances2[source_no-1][target_no-1]=100000000
A=numpy.array(distances2)
G=networkx.from_numpy_matrix(A,create_using=networkx.DiGraph())
#print (list(networkx.all_simple_paths(G,source_no,target_no,3)))
#print networkx.dijkstra_path_length(G,3,45)
a_path = networkx.dijkstra_path(G,source_no-1,target_no-1)
index_link2=[]
#print a_path
index_link2.append(a_path[0])

for i in range(len(a_path)-1):
    temp=[]
    distances2[a_path[i]][a_path[i+1]]=100000000
    distances2[a_path[i+1]][a_path[i]]=100000000
    V=numpy.array(distances2)
    G2=networkx.from_numpy_matrix(V,create_using=networkx.DiGraph())
    temp = networkx.dijkstra_path(G2,a_path[i],a_path[i+1])
    for j in range(len(temp)):
        if temp[j]!=a_path[i]:
            index_link2.append(temp[j])
            l=len(index_link2)
            distances2[index_link2[l-1]][index_link2[l-2]]=100000000
print "\nLINKING ARTICLES INDEX\n"
for i in range(len(index_link2)):
    print index_link2[i]+1


'''
a_path=[]
distances[source_no-1][target_no-1]=50
A=numpy.array(distances)
G=networkx.from_numpy_matrix(A,create_using=networkx.DiGraph())
#print (list(networkx.all_simple_paths(G,source_no,target_no,3)))
#print networkx.dijkstra_path_length(G,3,45)
a_path = networkx.astar_path(G,source_no-1,target_no-1,heuristic)
index_link=[]
#print a_path
index_link.append(a_path[0])

for i in range(len(a_path)-1):
    temp=[]
    distances[a_path[i]][a_path[i+1]]=50
    distances[a_path[i+1]][a_path[i]]=50
    V=numpy.array(distances)
    G2=networkx.from_numpy_matrix(V,create_using=networkx.DiGraph())
    temp = networkx.astar_path(G2,a_path[i],a_path[i+1],heuristic)
    for j in range(len(temp)):
        if temp[j]!=a_path[i]:
            index_link.append(temp[j])
            l=len(index_link)
            distances[index_link[l-1]][index_link[l-2]]=50
print "\nLINKING ARTICLES INDEX\n"
for i in range(len(index_link)):
    print index_link[i]+1
  
linkart = [];
'''
A=numpy.array(distances2)
G=networkx.from_numpy_matrix(A,create_using=networkx.DiGraph())
#print (list(networkx.all_simple_paths(G,source_no,target_no,3)))
#print networkx.dijkstra_path_length(G,3,45)
a_path = networkx.dijkstra_path(G,source_no-1,target_no-1)

print "\nLINKING ARTICLES INDEX\n"
print a_path

index_link2 = a_path
'''

distances2[source_no-1][target_no-1]=100000000
A=numpy.array(distances2)
G=networkx.from_numpy_matrix(A,create_using=networkx.DiGraph())
#print (list(networkx.all_simple_paths(G,source_no,target_no,3)))
#print networkx.dijkstra_path_length(G,3,45)
a_path = networkx.astar_path(G,source_no-1,target_no-1,heuristic2)
index_link2=[]
#print a_path
index_link2.append(a_path[0])

for i in range(len(a_path)-1):
    temp=[]
    distances2[a_path[i]][a_path[i+1]]=100000000
    distances2[a_path[i+1]][a_path[i]]=100000000
    V=numpy.array(distances2)
    G2=networkx.from_numpy_matrix(V,create_using=networkx.DiGraph())
    temp = networkx.astar_path(G2,a_path[i],a_path[i+1],heuristic2)
    for j in range(len(temp)):
        if temp[j]!=a_path[i]:
            index_link2.append(temp[j])
            l=len(index_link2)
            distances2[index_link2[l-1]][index_link2[l-2]]=100000000
print "\nLINKING ARTICLES INDEX\n"
for i in range(len(index_link2)):
    print index_link2[i]+1



'''

'''
for i in range(len(index_link)):
    linkart.append(sentences[index_link[i]])
'''
'''
for i in range(len(sentences)):
    linkart.append(sentences[i])

    
link_arr =[];
for i in range(len(linkart)):
    for j in range(len(linkart[i])):
        counter = 0
        for k in range(len(linkart)):
            for l in range(len(linkart[k])):
                 if linkart[i][j][0] == linkart[k][l][0] and i!=k:
                     counter = counter + linkart[k][l][1]
        link_arr.append([linkart[i][j][0],linkart[i][j][1]+counter])

link_arr2 =[]
#print link_arr

found = [0]*len(link_arr)

for i in range(len(link_arr)):
    for j in range(i,(len(link_arr))):
        if link_arr[i][0]==link_arr[j][0] and i!=j :
            found[j] = 1

#print link_arr
for i in range(len(link_arr)):
    if found[i]==0:
        link_arr2.append(link_arr[i])
'''
#print link_arr2
'''
f = open("words.txt","w")
json.dump(link_arr2,f)
f.close()
'''
f=open("words.txt","r")
link_arr2=json.load(f)
f.close()

corpus =[]
'''
for i in range(len(sentences)):
    corpus.append([])  
    for j in range(len(sentences[i])):
      for k in range(len(link_arr2)):
          if sentences[i][j][0] == link_arr2[k][0]:
              corpus[i].append([k,sentences[i][j][1]])

f = open("corpus.txt","w")
json.dump(corpus,f)
f.close()
'''



print '''--------------------------------------------------COHERENCE 2------------------------------------'''
temp3 = []
temp =1
for i in range(len(index_link)-1):
    temp3.append(new[index_link[i]][index_link[i+1]][1])
    temp = min (temp, new[index_link[i]][index_link[i+1]][1])
    print new[index_link[i]][index_link[i+1]][1]

print temp*100,'% is the coherence'
print temp3
coherence2=temp*100
eg_bhatta = temp3
print "\n"
print "\n"

print '''--------------------------------------------------COHERENCE 1-----------------------------------'''
temp2 = []
for i in range(len(index_link)-1):
    temp = 0
    for j in range(len(sentences[index_link[i]])):
        for k in range(len(sentences[index_link[i+1]])):
            if sentences[index_link[i]][j][0] == sentences[index_link[i+1]][k][0]:
                temp = temp +  1
    temp2.append(temp)    

print min(temp2),'is the coherence'
print temp2
coherence1 = min(temp2)

print "\n"
print "\n"
print '''--------------------------------------------------COHERENCE 2------------------------------------'''
temp3 = []
temp =1
for i in range(len(index_link2)-1):
    temp3.append(new[index_link2[i]][index_link2[i+1]][1])
    temp = min (temp, new[index_link2[i]][index_link2[i+1]][1])
    print new[index_link2[i]][index_link2[i+1]][1]

print temp*100,'% is the coherence'
print temp3
coherence4=temp*100

eg_cosine = temp3
print "\n"
print "\n"

print '''--------------------------------------------------COHERENCE 1------------------------------------'''
temp2 = []
for i in range(len(index_link2)-1):
    temp = 0
    for j in range(len(sentences[index_link2[i]])):
        for k in range(len(sentences[index_link2[i+1]])):
            if sentences[index_link2[i]][j][0] == sentences[index_link2[i+1]][k][0]:
                temp = temp +  1
    temp2.append(temp)    

print min(temp2),'is the coherence'
print temp2
coherence3 = min(temp2)
eg_coherence1 = temp2
print "\n"
print "\n"


        
#counter =counter +1
#print(networkx.astar_path(G,source_no,target_no,heuristic))
'''
st_distance=[]
initial_st_distance=[]
for i in range(len(source_distance)):
    st_distance.append(source_distance[i]*target_distance[i])
    initial_st_distance.append(source_distance[i]*target_distance[i])


st_distance.sort()
final_st_distance = st_distance

#print initial_st_distance
#print final_st_distance

found = [0]*len(source_distance)
index_st_distance = []

for i in range(len(final_st_distance)):
	for j in range(len(initial_st_distance)):
		if final_st_distance[i]==initial_st_distance[j] and found[j]==0:
			index_st_distance.append(j)
			found[j]=1

#print index_st_distance
index_link=[]

max_st_distance = max(st_distance)
for i in range(len(index_st_distance)):
    if st_distance[i]>=factor_chosen*max_st_distance:
        index_link.append(index_st_distance[i])

#print index_link
source_distance_link=[]
initial_source_distance_link=[]
final_source_distance_link=[]

for i in range(len(index_link)):
    source_distance_link.append(source_distance[index_link[i]])
    initial_source_distance_link.append(source_distance[index_link[i]])

#print source_distance_link
source_distance_link.sort(reverse=True)
final_source_distance_link = source_distance_link
found =[0]*len(source_distance_link)
#print source_distance_link

index_link_distance = []
for i in range(len(final_source_distance_link)):
	for j in range(len(initial_source_distance_link)):
		if final_source_distance_link[i]==initial_source_distance_link[j] and found[j]==0:
			index_link_distance.append(j)
			found[j]=1
#print index_link_distance

link = []

for i in range(len(index_link_distance)):
    link.append(index_link[index_link_distance[i]])


'''
'''
sentences2 =[]
for i in range(1,41):
    file=open('C:/Users/prakash chand gupta/Desktop/AI PROJECT/Articles2/'+str(i)+'.txt','r')
    p=file.read().replace("'",' ').replace("`",'').replace(";",'').replace("!",'.').lower()
    temp = split(p)
    temp3 = [x for x in temp if x not in commonwords]
    sentences2.append(list(set(temp3)))

words=[]
for i in range(len(sentences2)):
    words = list(set(list(set(words)) + list(set(sentences2[i]))))

wordgraph=[]
for i in range(len(words)):
    wordgraph.append([])
    for j in range(len(sentences2)):
        for k in range(len(sentences2[j])):
            if words[i]==sentences2[j][k]:
                wordgraph[i].append(j)


#print wordgraph

bp_graph=[]
for i in range(len(wordgraph)):
    bp_graph.append([])
    for j in range(len(sentences2)):
        bp_graph[i].append([0])
        
for i in range(len(wordgraph)):
    for j in range(len(wordgraph[i])):
        bp_graph[i][wordgraph[i][j]]=1
        
print bp_graph
'''
'''
print "\n--------SOURCE-----------\n"
print source_p
print "\n"
for i in range(len(link)):
    print link[i]
    print "\n"
    print link_p[link[i]]
    print "\n"
print "\n--------TARGET-----------\n"    
print target_p
print "\n"
print link
for i in range(len(link)):
    if link[i]<source_no-501 and link[i]>=target_no-501:
        link[i]=link[i]+1
    elif link[i]>=source_no-501 and link[i]<target_no-501:
        link[i]=link[i]+1
    elif link[i]>=source_no-501 and link[i]>=target_no-501:
        link[i]=link[i]+2
    elif link[i]<source_no-501 and link[i]<target_no-501:
        link[i]=link[i]
print "-------------LINKING ARTICLES-----------"
print link        

'''

print "\nLINK TOPICS\n"
for i in range(len(index_link)):
    print topic[index_link[i]]

print "\nLINK TOPICS\n"
for i in range(len(index_link2)):
    print topic[index_link2[i]]
