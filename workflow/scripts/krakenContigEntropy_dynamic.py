#!/usr/bin/env python3

import os
import sys
import argparse
import operator
import math
import numpy as np

parser = argparse.ArgumentParser(description='Process output of Kraken run on contigs.')
parser.add_argument('-c','--contig', default="", help='per-contig output from Kraken')
parser.add_argument('-r','--report', default="", help='report output from Kraken (all taxa must be reported)')
parser.add_argument('-e','--entropy', default=0.0,type=float,help='upper cut-off for entropy, defaults to 0.0')
parser.add_argument('-u','--unknown', action='store_true',help='set -u, if contigs should be reported that were not annotated by Kraken')
parser.add_argument('-o','--outname', default="",help='name for the output can be specified or will be constructed from input file name and entropy cut-off')
parser.add_argument('-s','--silent', action='store_false',help='set -s, to suppress printing the number of contigs not annotated by Kraken')
#parser.add_argument('-d','--taxdir', default=".",type=str,help='directory where the tax files sit')

args = parser.parse_args()
krakFile = args.contig
taxFile = args.report
divthresh = args.entropy

#taxdir = args.taxdir+"/"
if args.outname != "":
    outFile1 = args.outname
elif divthresh>0:
    outFile1 = krakFile + "annoEnt" + str(divthresh) + ".tsv"
else:
    outFile1 = krakFile + "annoUnambig.tsv"


# Definition of the class Node

class Node:
    """Node"""
    def __init__(self):
        self.tax_id = 0       # Number of the tax id.
        self.parent = 0       # Number of the parent of this node
        self.children = []    # List of the children of this node
        self.tip = 0          # Tip=1 if it's a terminal node, 0 if not.
        self.name = ""        # Name of the node: taxa if it's a terminal node, numero if not.       
    def genealogy(self):      # Trace genealogy from root to leaf
        ancestors = []        # Initialise the list of all nodes from root to leaf.
        tax_id = self.tax_id  # Define leaf
        while 1:
            if tax_id in name_object:
                ancestors.append(tax_id)
                tax_id = name_object[tax_id].parent
            else:
                break
            if tax_id == "1":
                # If it is root, we reached the end.
                # Add it to the list and break the loop
                ancestors.append(tax_id)
                break
        return ancestors # Return the list

# Function to find common ancestor between two nodes or more
def common_ancestor(node_list):
    global name_object
    list1 = name_object[node_list[0]].genealogy()  # Define the whole genealogy of the first node
    for node in node_list:
        list2 = name_object[node].genealogy()      # Define the whole genealogy of the second node
        ancestral_list = []                             
        for i in list1:
            if i in list2:                         # Identify common nodes between the two genealogy
                ancestral_list.append(i)                 
        list1 = ancestral_list                     # Reassing ancestral_list to list 1.
    common_ancestor = ancestral_list[0]            # Finally, the first node of the ancestral_list is the common ancestor of all nodes.
    return common_ancestor                         # Return a node


#############################
#                           #
#   Read taxonomy file      #
#                           #
#############################

global name_object
name_object = {}
name_dict = {}          # Initialise dictionary with TAX_ID:NAME
name_dict_reverse = {}  # Initialise dictionary with NAME:TAX_ID #not all names are unique, though
parentList = [None] * 100
rank_dict = {}          # Initialise dictionary with TAX_ID:RANK
rank_dict_reverse = {}  # Initialise dictionary with RANK:[TAX_IDs]

tax_file = open(taxFile, "r")
old_ilevel = 0
while 1:
    line = tax_file.readline()
    if line == "":
        break
    line = line.rstrip()
    tab = line.split("\t")   # 0: total perc, 1:total counts, 2: counts at this level, 3: rank abbrev., 4: tax ID, 5: name (indented)
    if tab[3] == "U":
        next
    rank, tax_id, indentName = tab[3], tab[4], tab[5]      # Assign tax_id and name ...
    name = indentName.lstrip(' ')
    ilevel = int((len(indentName) - len(name))/2)     # get level from indentation
#    print(ilevel)
    parentList[ilevel] = tax_id
    if ilevel > 0:
        tax_id_parent = parentList[ilevel-1]
    if name not in name_dict_reverse:
        name_dict_reverse[name] = tax_id
        name_dict[tax_id] = name                # ... and load them into dictionary
    else:
        if ilevel > 0:
            name_dict_reverse[name + "." + name_dict[tax_id_parent]] = tax_id
            name_dict[tax_id] = name  + "." + name_dict[tax_id_parent]  # ... for kids that have their parent's names or worse
                        
    rank_dict[tax_id] = rank
    if rank not in rank_dict_reverse:
        rank_dict_reverse[rank] = [tax_id]
    else:
        rank_dict_reverse[rank].append(tax_id)
        
        
    if tax_id not in name_object:          # this should always be the case - otherwise we would overwrite in the following lines
        name_object[tax_id] = Node()
    name_object[tax_id].tax_id = tax_id             # Assign tax_id
    if ilevel > 0:
        name_object[tax_id].parent = tax_id_parent  # Assign tax_id parent
    name_object[tax_id].name = name                 # Assign name
    
    if ilevel > 0:
        siblings = name_object[tax_id_parent].children  # Parent is always already in the object
        siblings.append(tax_id)                         # ...we found its children.
        name_object[tax_id_parent].children = siblings  # ... so add them to the parent 
tax_file.close()

########################
#                      #
# reconstruct taxonomy #
#                      #
########################
convention = ["R","D","K","P","C","O","F","G","S"]

ranklist = [rank_dict["1"]]
indx1 = convention.index(ranklist[-1])
all_ranks = list(rank_dict_reverse.keys())
all_main_ranks = [s[0] for s in all_ranks]
#print(" ".join(all_ranks))
#print(" ".join(all_main_ranks))
for main_level in convention[indx1:]:
    curr_ranks = [all_ranks[i] for i, e in enumerate(all_main_ranks) if e == main_level]
    curr_ranks.sort()
    ranklist = [*ranklist, *curr_ranks]
print("All ranks: " + " ".join(ranklist))    

leftover_ranks = np.setdiff1d(all_ranks,ranklist)

print("Ranks in report that were not used: " + " ".join(leftover_ranks))
   

#####################
#                   #
# contig annotation #
#                   #
#####################
#function to calculate diversity of annotations
def shannonDiv(dictionary,sumTax):
    taxDiv = 0.0
    if len(dictionary) > 0 and sumTax > 1:
        for item in dictionary:
            taxDiv += (float(dictionary[item])/sumTax) * math.log(float(dictionary[item])/sumTax,2) / math.log(sumTax,2)
    else:
        taxDiv = 1.0
    return 0.00 - taxDiv

# function to retrieve name and number of annotated bases for a taxon
def orgGenealCount(anc,taxDict,orgCnt,geneaDict,taxSum):
    if anc in geneaDict:
        taxName = geneaDict[anc]
        taxSum += int(orgCnt)
        if taxName not in taxDict:
            taxDict[taxName] = int(orgCnt)
        else:
            taxDict[taxName] += int(orgCnt)
    return taxDict, taxSum

# function to retrieve taxon name
def orgGenealName(anc,geneaDict,taxName):
    if anc in geneaDict:
        taxName = geneaDict[anc]
    return taxName

# function to test if lower taxon is in higher taxon and stop the annotation if necessary
def phyloTester(taxName,testList,retVal,annotationList):
    if taxName != "unknown":
        testList.append(name_dict_reverse[taxName])
        if len(testList) == 2:
            if testList[0] in name_object[testList[1]].genealogy():
                del testList[0]
            else:
                del annotationList[-1]
                retVal = -2.0
                taxName = "unknown"
                    
ucnt = 0

krak_file =  open(krakFile,"r")
out_file1 = open(outFile1, "w")
out_file1.write("contig" + "\t"+ "length" +"\t"+ "entropy" + "\t"+ "annotationLevel" + "\t"+ "\t".join(ranklist) +"\n")
while 1:
    linek = krak_file.readline()
    if linek == "":
        break
    linek = linek.rstrip()
    tabk = linek.split("\t")
    if tabk[0] == "U":
        if args.unknown:
            out_file1.write(tabk[1] + "\t"+ tabk[3] +"\t"+ "NA" + "\t"+ "not annotated" + "\t"+ "\t".join(["NA"] * len(ranklist)) +"\n")
        ucnt += 1
    else:
#        print(tabk[0])
        cotak = tabk[4].split(" ")
        orgs = {}
        for i in cotak:
            orgID = i.split(":")[0]
            orgCount = i.split(":")[1]
            if orgID != "A" and orgID != "0":
#                if orgID not in orgs:
#                    orgs[orgID] = int(orgCount)
#                else:
#                    orgs[orgID] += int(orgCount)
                if orgID in name_dict:
                    for anc in name_object[orgID].genealogy(): #should return taxids
#                        print(rank_dict[anc] + " " + anc)
                        if anc not in orgs:
                            orgs[anc] = int(orgCount)
                        else:
                            orgs[anc] += int(orgCount)
                else:
                    print("Unknown tax ID " + orgID + " !")
        ranksum_dict = {}
        rank_orgs = {}
        for ctax in orgs.keys():
#            print(rank_dict[ctax] + " " + ctax)
#            print(rank_dict[ctax])
            if rank_dict[ctax] not in ranksum_dict:
                ranksum_dict[rank_dict[ctax]] = [orgs[ctax],len(ranklist)-ranklist.index(rank_dict[ctax])]
                rank_orgs[rank_dict[ctax]] = {ctax: orgs[ctax]}
            else:
                ranksum_dict[rank_dict[ctax]][0] += orgs[ctax]
                rank_orgs[rank_dict[ctax]][ctax] = orgs[ctax]
 #       print(ranksum_dict)
        cranks = sorted(ranksum_dict.items(), key=lambda x:(x[1][0],x[1][1]), reverse = True)
        rankwin_dict = {}
        winnerList = []
        first = True
        for rank in cranks:
#            print(rank)
            rank = rank[0]
            cdiv = shannonDiv(rank_orgs[rank],ranksum_dict[rank][0])
            if cdiv <= divthresh and cdiv >= 0:
                cwin = max(rank_orgs[rank].items(), key=operator.itemgetter(1))[0]
                winancs = name_object[cwin].genealogy()
                if first:
                    rankwin_dict[rank] = cwin
                    winnerList.append(cwin)
                    first = False
                elif not set(winnerList).isdisjoint(set(winancs)):
                    rankwin_dict[rank] = cwin
                    winnerList.append(cwin)
            else:
                break
        #out_file1.write(tabk[1] + "\t"+ tabk[3] +"\t"+ "NA" + "\t"+ "not annotated" + "\t"+ "\t".join(["NA"] * len(ranklist)) +"\n")
        ranked_winner_list = []
        went = "NA"
        wrank = "not annotated"
        for crank in ranklist:
            if crank in rankwin_dict:
                ranked_winner_list.append(name_dict[rankwin_dict[crank]])
                went = shannonDiv(rank_orgs[crank],ranksum_dict[crank][0])
                wrank = crank
            else:
                ranked_winner_list.append("NA")
        if went != "NA" or args.unknown:
            out_file1.write(tabk[1] + "\t"+ tabk[3] +"\t"+ str(went) + "\t"+ wrank + "\t"+ "\t".join(ranked_winner_list) +"\n")        
krak_file.close()
out_file1.close()

if args.silent:
    print(ucnt)

