# -*- coding: utf-8 -*-
"""
Created on Oct 3 2023
author: Amanda Erin Wilson
@author: @amandaerinwilson

Code associated with data and figures for Wilson AE, Liberles DA. Expectations of Duplicate Gene Retention Under the Gene Duplicability Hypothesis.
Purpose: 
    1) Generate the 3D Surface Plots t1 vs t2 vs pratio
    2) To create csv files that contain a 3D and 2D version of the arrays that contain the pratios for various t1 and t2 values.
    3) Generate the survival curves over time for each category Alt_func, Dos, and Non.
    
B-F parameters:
F+d (rate at which fully redundant genes get lost from the genome) to d (d is the rate at which non duplicated genes are lost)
B and C describe the dynamics from how you move from instantaneous rate to the asymptotic rate (shape of the curve - approx of how [process behaves)
For genes in the Non category, where both copies can only be retained by chance, the parameter values for the survival curve are b = 0, c = 1, d > 10, f = anything. 
For genes in the Dos category, that are sensitive to dosage balance effects, the parameter values are b < 0, 0 < c < 1, d = -f for λ(t)0.02 < 0.1. 
For genes in the Alt_func category, that have potential to subfunctionalize or neofunctionalize, the parameter values are b > 0, c > 0, d > 0, f > 0.   

Additional Parameters
Alt_func : percent of starting genome in Alt_func category (Alpha_Alt_func)
Dos : percent of starting genome in Dos category (Alpha_Dos)
Non : percent of starting genome in Non category (Alpha_Non)
switch : percent of gene duplicates retained through ALt_func that switch to Non category in t2. This value is 0 under the gene duplicability model and >0 for mutational opportunity model (beta_switch_mo)

"""
import math
import csv
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd

###########################################################################
#initialize parameters
#n_max = 170
n_max = 100

#number of time points in t1 and t2 
time_points = 51
#time_points = 51
#time_points = 81 #may be more clear

file_name_start = 'mutational_opportunity_vers8_' 

#make graph using the pratio OR the log 10 of pratio
p_ratio = 'pratio'
#p_ratio = 'log of pratio'
###############################################################################
#CHOOSE ONE OF THE FOLLOWING
"""
#1)
#INDEPENDENCE HYPOTHESIS
alts = [0]
doses = [0]
nons = [1]
number_of_combos = 1

switch = 0
file_name_end = '_independence'
hypothesis  = "Independence"
"""
"""
#2)
#FOR TESTING ONE GRAPH FOR DUPLICABILITY
alts = [0.3]
doses = [0.45]
nons = [0.25]
number_of_combos = 1

switch = 0
file_name_end = '_duplicability'
hypothesis  = "Duplicability"
"""

"""
#3
#FOR TESTING ONE GRAPH FOR MUTATIONAL OPPORTUNITY
alts = [0.3]
doses = [0.45]
nons = [0.25]
number_of_combos = 1

switch = .75
file_name_end = "_mut_op"
hypothesis = "Mutational Opportunity"
"""

"""
#4
#FOR TESTING MULTIPLE GRAPHS FOR DUPLICABILITY
alts = [0.75, 0.60, 0.45, 0.30, 0.15, 0.00, 
        0.50, 0.40, 0.30, 0.20, 0.10, 0.00, 
        0.25, 0.10, 0.00]

doses = [0.00, 0.15, 0.30, 0.45, 0.60, 0.75, 
         0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 
         0.00, 0.15, 0.25]
nons = [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
        0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 
        0.75, 0.75, 0.75]
number_of_combos = 15

switch = 0
file_name_end = '_duplicability'
hypothesis  = "Duplicability"
"""
#5
#FOR TESTING MULTIPLE GRAPHS FOR MUTATIONAL OPPORTUNITY
alts = [0.75, 0.60, 0.45, 0.30, 0.15, 0.00, 
        0.50, 0.40, 0.30, 0.20, 0.10, 0.00, 
        0.25, 0.10, 0.00]

doses = [0.00, 0.15, 0.30, 0.45, 0.60, 0.75, 
         0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 
         0.00, 0.15, 0.25]
nons = [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
        0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 
        0.75, 0.75, 0.75]
number_of_combos = 15

#switch = 0.75
switch = 0.25
#switch = 0.20
#switch = 0.1
#switch = 0.5
#switch = 0.05
file_name_end = "_mut_op"
hypothesis = "Mutational Opportunity"

###############################################################################

#parameters for Neo-functionalization/Sub-functionalization
b_alt_func = 10.0
c_alt_func = 2.37
d_alt_func = 0.00054
f_alt_func = 5.84

#parameters for Dosage
b_dos = -17.0
c_dos = 0.2573
d_dos = -0.000028
f_dos = 0.000028

#parameters for Non-functionalization
b_non = 0
c_non = 1
d_non = 20
f_non = 5


###########################################################################
#Functions

def calculate_probability_of_survival_of_duplicate_gene_copy_by_time(b, c, d, f, time):  
    summation = 0
    for n in range(0,n_max):
        nfac = math.factorial(n)
        beta = (((-b)**n)*(time**((c*n)+1)))/((c*n*nfac) + nfac)
        summation = summation + beta
        # print("summation: " + str(summation))
    survival_probability = math.exp(-d*time - f*summation)
    return survival_probability
    
def calculate_pratio_2d(st1_alt_func, st1_dos, st1_non, st2_alt_func, st2_dos, st2_non, alt_func_percent, dos_percent, non_percent, alt_switch_percent):    
    #probability of (survival in t2 given survived in t1)/(survival in t2 given lost in t1)
    alt_ret =  (2*alt_func_percent*st1_alt_func)
    alt_ret_ret_switch = alt_ret*st2_non *alt_switch_percent
    alt_ret_ret_noswitch = alt_ret*st2_alt_func * (1-alt_switch_percent)
    dos_ret = (2*dos_percent*st1_dos)
    dos_ret_ret = dos_ret * st2_dos   
    non_ret = (2*non_percent*st1_non)
    non_ret_ret = non_ret * st2_non
    alt_noret = ((1-st1_alt_func)*alt_func_percent)
    alt_noret_ret = alt_noret * st2_alt_func
    dos_noret = ((1-st1_dos)*dos_percent)
    dos_noret_ret = dos_noret * st2_dos    
    non_noret = ((1-st1_non)*non_percent)
    non_noret_ret = non_noret * st2_non
    pratio = ((alt_ret_ret_noswitch + alt_ret_ret_switch + dos_ret_ret + non_ret_ret)/(alt_noret_ret + dos_noret_ret + non_noret_ret)) * ((alt_noret + dos_noret + non_noret)/(alt_ret + dos_ret + non_ret))                                 
    return pratio  

def read_csv_file(file_name):
    file_name_full = (file_name_start + file_name +'.csv')
    df = pd.read_csv(file_name_full,
            header=0,
            usecols=["t1", "t2", "pratio","alt_surv_t1", "dos_surv_t1", "non_surv_t1", "alt_surv_t2", "dos_surv_t2", "non_surv_t2", "log of pratio"])    
    #print(df.head())    
    return df
    
def print_3d_graph(percents, file_name, p_ratio):
    df = read_csv_file(file_name)
    minimum_pratio = df['pratio'].min()
    print("Minimum Pratio: " + str(minimum_pratio))  
    #plot 3D scatter
    fig2 = plt.figure()
    ax2 = plt.axes(projection='3d')
    surf2 = ax2.scatter3D(df['t1'], df['t2'], df[p_ratio], c = df[p_ratio], cmap=cm.cividis)
    fig2.colorbar(surf2)
    ax2.set_title('Pratio for t1 and t2: ' + percents, fontsize=14)
    ax2.set_xlabel('$t1$', fontsize=12)
    ax2.set_ylabel('$t2$', fontsize=12)
    ax2.set_zlabel(r'probability ratio', fontsize=11)
    ax2.view_init(15, 45)
    plt.draw()
    plt.pause(.001)
    
def plot_survival_curves(file_name):
    df = read_csv_file(file_name)
    #Plot survival curves
    fig, t1_plot = plt.subplots()
    t1_plot.plot(df['t1'], df['alt_surv_t1'], color = 'red', label = 'Alt_func')
    t1_plot.plot(df['t1'], df['dos_surv_t1'], color = 'blue', label = 'Dos')
    t1_plot.plot(df['t1'], df['non_surv_t1'], color = 'yellow', label = 'Non')
    t1_plot.legend(loc = 'upper right', shadow = True, fontsize = '12')
    t1_plot.set_title('Survival over Time (t1)', fontsize=14)
    t1_plot.set_xlabel('Time Since Duplication Event', fontsize=12)
    t1_plot.set_ylabel('Proportion Gene Duplicate Copies Surviving', fontsize=12)
    plt.show()   
    
def calculate_and_make_csv(file_name, alt_func_percent, dos_percent, non_percent, alt_switch_percent):
    file = open(file_name_start+ file_name +'.csv', 'w', newline='')
    writer = csv.writer(file, delimiter=',')
    write_header_row = ["t1", "t2", "pratio", "alt_surv_t1", "dos_surv_t1", "non_surv_t1", "alt_surv_t2", "dos_surv_t2", "non_surv_t2", "log of pratio"]
    writer.writerow(write_header_row)
    each_t1 = 0.01
    for i in range(0, time_points):  
        alt_func_survival_t1 = calculate_probability_of_survival_of_duplicate_gene_copy_by_time(b_alt_func, c_alt_func, d_alt_func, f_alt_func, each_t1)
        dos_survival_t1 = calculate_probability_of_survival_of_duplicate_gene_copy_by_time(b_dos, c_dos, d_dos, f_dos, each_t1)
        non_survival_t1 = calculate_probability_of_survival_of_duplicate_gene_copy_by_time(b_non, c_non, d_non, f_non, each_t1)
        each_t2 = 0.01
        for j in range(0, time_points):
            #print("t1: " + str(each_t1))
            #print("t2: " + str(each_t2))        
            alt_func_survival_t2 = calculate_probability_of_survival_of_duplicate_gene_copy_by_time(b_alt_func, c_alt_func, d_alt_func, f_alt_func, each_t2)
            dos_survival_t2 = calculate_probability_of_survival_of_duplicate_gene_copy_by_time(b_dos, c_dos, d_dos, f_dos, each_t2)
            non_survival_t2 = calculate_probability_of_survival_of_duplicate_gene_copy_by_time(b_non, c_non, d_non, f_non, each_t2)       
            probability_ratio_2d = calculate_pratio_2d(alt_func_survival_t1, dos_survival_t1, non_survival_t1, alt_func_survival_t2, dos_survival_t2, non_survival_t2, alt_func_percent, dos_percent, non_percent, alt_switch_percent)
            #print("Probability ratio: " + str(probability_ratio_2d))
            log_pratio = math.log10(probability_ratio_2d)
            row_to_write = [each_t1, each_t2, probability_ratio_2d, alt_func_survival_t1, dos_survival_t1, non_survival_t1, alt_func_survival_t2, dos_survival_t2, non_survival_t2, log_pratio]
            writer.writerow(row_to_write)
            each_t2 = each_t2+0.01
        each_t1 = each_t1+0.01                     
    file.close()
#############################################################################

def main(percent_alt_func, percent_dos, percent_non, percent_alt_switch):    
    percentages = str(100*switch) +"% switch,  \n" +str(100*percent_alt_func) +"% Alt_func, "+str(100*percent_dos)+"% Dos, "+str(100*percent_non)+"% Non, \n (" + hypothesis + " Hypothesis)"
    percentages_file_name = str(100*percent_alt_func)+'_'+str(100*percent_dos)+'_'+str(100*percent_non)+file_name_end
    print(percentages)
    calculate_and_make_csv(percentages_file_name, percent_alt_func, percent_dos, percent_non, percent_alt_switch)
    print_3d_graph(percentages, percentages_file_name, p_ratio)
    

#############################################################################

for i in range (0, number_of_combos):
    alt = alts[i]
    dos = doses[i]
    non = nons[i]
    main(alt, dos, non, switch)

file_name_for_surv_curv = str(100*alts[0])+'_'+str(100*doses[0])+'_'+str(100*nons[0])+file_name_end
plot_survival_curves(file_name_for_surv_curv)

########################################
