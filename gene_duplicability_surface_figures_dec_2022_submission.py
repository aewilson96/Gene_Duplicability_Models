# -*- coding: utf-8 -*-
"""
Published on Github Dec 2 2022
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

"""

import math
#import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import numpy as np
import csv

###########################################################################
#initialize parameters
# n_max = 170
n_max = 100

# q = time points
q = 51 #number of time points in t1 and t2 
#q = 81 #may be more clear

#number_percent_combos = 1
number_percent_combos = 16
percentages = [' 100% Alt_func, 0% Dos, 0% Non \n (Independence Hypothesis)', ' 75% Alt_func, 0% Dos, 25% Non \n (Duplicability Hypothesis)', ' 60% Alt_func, 15% Dos, 25% Non \n (Duplicability Hypothesis)', ' 45% Alt_func, 30% Dos, 25% Non \n (Duplicability Hypothesis)', ' 30% Alt_func, 45% Dos, 25% Non \n (Duplicability Hypothesis)', ' 15% Alt_func, 60% Dos, 25% Non \n (Duplicability Hypothesis)', ' 0% Alt_func, 75% Dos, 25% Non \n (Duplicability Hypothesis)', ' 50% Alt_func, 0% Dos, 50% Non \n (Duplicability Hypothesis)', ' 40% Alt_func, 10% Dos, 50% Non \n (Duplicability Hypothesis)', ' 30% Alt_func, 20% Dos, 50% Non \n (Duplicability Hypothesis)', ' 20% Alt_func, 30% Dos, 50% Non \n (Duplicability Hypothesis', ' 10% Alt_func, 40% Dos, 50% Non \n (Duplicability Hypothesis)', ' 0% Alt_func, 50% Dos, 50% Non \n (Duplicability Hypothesis)', ' 25% Alt_func, 0% Dos, 75% Non \n (Duplicability Hypothesis)', ' 10% Alt_func, 15% Dos, 75% Non \n (Duplicability Hypothesis)', ' 0% Alt_func, 25% Dos, 75% Non \n (Duplicability Hypothesis)']
percentages_file_name = ['100_Alt_func_0_Dos_0_Non', '75_Alt_func_0_Dos_25_Non', '60_Alt_func_15_Dos_25_Non)', '45_Alt_func_30_Dos_25_Non', '30_Alt_func_45_Dos_25_Non', '15_Alt_func_60_Dos_25_Non', '0_Alt_func_75_Dos_25_Non', '50_Alt_func_0_Dos_50_Non', '40_Alt_func_10_Dos_50_Non', '30_Alt_func_20_Dos_50_Non', '20_Alt_func_30_Dos_50_Non', '10_Alt_func_40_Dos_50_Non', '0_Alt_func_50_Dos_50_Non', '25_Alt_func_0_Dos_75_Non', '10_Alt_func_15_Dos_75_Non', '0_Alt_func_25_Dos_75_Non']
alt_func_percentages = [1.0, 0.75, 0.6, 0.45, 0.30, 0.15, 0.0, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0, 0.25, 0.1, 0.0]
dos_percentages = [0.0, 0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.0, 0.15, 0.25]
non_percentages = [0.0, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75]

survival_immediately_post_wgd = 0.9999999999999 #needs to not be 1 for calculation, and can make sense because perhaps can assume two wgd events can't happen at exactly the same time, so SOMETHING had to be lost


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

def calculate_probability_of_survival_of_duplicate_gene_copy_by_time(b, c, d, f, time_range):  
    list_of_probabilities_of_survival = []
    for time in time_range:
        summation = 0
        for n in range(0,n_max):
            nfac = math.factorial(n)
            beta = (((-b)**n)*(time**((c*n)+1)))/((c*n*nfac) + nfac)
            summation = summation + beta
#        print(summation)
        survival_probability = math.exp(-d*time - f*summation)
        list_of_probabilities_of_survival.append(survival_probability)
    return list_of_probabilities_of_survival

def calculate_st2_ret1 (percent_of_genome, st1_ret1):
    st2_ret1 = 2*percent_of_genome*st1_ret1
    return st2_ret1

def calculate_st2_notret1 (percent_of_genome, st1_ret1):
    st2_notret1 = (1-st1_ret1)*percent_of_genome
    return st2_notret1
        
def calculate_pratio(st1_alt_func, st1_dos, st1_non, st2_alt_func, st2_dos, st2_non):
    #Calculate s(t2) given retained/not in t1
    prob_surv = np.zeros(shape=(q,q))
    
    for i in range(0,q):
        #calculating st2 given retained in t1
        st2_ret1_alt_func = calculate_st2_ret1(alt_func_percent, st1_alt_func[i])
        st2_ret1_dos = calculate_st2_ret1(dos_percent, st1_dos[i])
        st2_ret1_non = calculate_st2_ret1(non_percent, st1_non[i])
        
        #calculating st2 given not retained in t1
        st2_notret1_alt_func = calculate_st2_notret1(alt_func_percent, st1_alt_func[i])  
        st2_notret1_dos = calculate_st2_notret1(dos_percent, st1_dos[i])
        st2_notret1_non = calculate_st2_notret1(non_percent, st1_non[i])
        
        for j in range(0,q):
            #probability of (survival in t2 given survived in t1)/(survival in t2 given lost in t1)
            psurv = (((st2_ret1_alt_func * st2_alt_func[j])+(st2_ret1_dos * st2_dos[j])+(st2_ret1_non * st2_non[j]))/((st2_notret1_alt_func * st2_alt_func[j])+(st2_notret1_dos * st2_dos[j])+(st2_notret1_non * st2_non[j]))) * (((st2_notret1_alt_func)+(st2_notret1_dos)+(st2_notret1_non))/((st2_ret1_alt_func)+(st2_ret1_dos)+(st2_ret1_non)))
            prob_surv[i][j] = psurv
            
#    print("this is prob_surv")
#    print(prob_surv)
    return prob_surv
    
def calculate_pratio_2d(t1, t2, st1_alt_func, st1_dos, st1_non, st2_alt_func, st2_dos, st2_non):
    #Calculate s(t2) given retained/not in t1
#    prob_surv_2d_array = np.array(['t1', 't2', 'prob_ratio'])
    prob_surv_2d_array = np.empty([1, 3])
    prob_surv_2d_array[0, 0] = 0
    prob_surv_2d_array[0, 1] = 1
    prob_surv_2d_array[0, 2] = 2

    for i in range(0,q):
        #calculating st2 given retained in t1
        st2_ret1_alt_func = calculate_st2_ret1(alt_func_percent, st1_alt_func[i])
        st2_ret1_dos = calculate_st2_ret1(dos_percent, st1_dos[i])
        st2_ret1_non = calculate_st2_ret1(non_percent, st1_non[i])
        
        #calculating st2 given not retained in t1
        st2_notret1_alt_func = calculate_st2_notret1(alt_func_percent, st1_alt_func[i])  
        st2_notret1_dos = calculate_st2_notret1(dos_percent, st1_dos[i])
        st2_notret1_non = calculate_st2_notret1(non_percent, st1_non[i])
                
        for j in range(0,q):
            #probability of (survival in t2 given survived in t1)/(survival in t2 given lost in t1)
            psurv = (((st2_ret1_alt_func * st2_alt_func[j])+(st2_ret1_dos * st2_dos[j])+(st2_ret1_non * st2_non[j]))/((st2_notret1_alt_func * st2_alt_func[j])+(st2_notret1_dos * st2_dos[j])+(st2_notret1_non * st2_non[j]))) * (((st2_notret1_alt_func)+(st2_notret1_dos)+(st2_notret1_non))/((st2_ret1_alt_func)+(st2_ret1_dos)+(st2_ret1_non)))
            p_survs_array_part = np.empty([1, 3])
            t1_value = t1[i]
            t2_value = t2[j]
            p_survs_array_part[0, 0] = t1_value
            p_survs_array_part[0, 1] = t2_value
            p_survs_array_part[0, 2] = psurv                           
            prob_surv_2d_array = np.concatenate((prob_surv_2d_array, p_survs_array_part), axis = 0)                               
           
#    print("this is prob_surv")
#    print(prob_surv)
    return prob_surv_2d_array    

###########################################################################
def main_calculations(alt_func_percent, dos_percent, non_percent):
    t1_range = []
    for i in range(0,q):
        t1_range.append(i/100)
#    print(t1_range)
    
    #calculate s(t1) for alt_functionalization
    list_of_probabilities_of_survival_by_alt_functionalization = calculate_probability_of_survival_of_duplicate_gene_copy_by_time(b_alt_func, c_alt_func, d_alt_func, f_alt_func, t1_range)
#    print ("List of prob of surv by alt_func in t1: ")
#    print (list_of_probabilities_of_survival_by_alt_functionalization)
    list_of_probabilities_of_survival_by_alt_functionalization_in_t1 = []
    list_of_probabilities_of_survival_by_alt_functionalization_in_t1.append(survival_immediately_post_wgd)
    for each_time_point in range(1,q):
        list_of_probabilities_of_survival_by_alt_functionalization_in_t1.append(list_of_probabilities_of_survival_by_alt_functionalization[each_time_point])
#    print ("List of prob of surv by alt_func in t1 adjusted to avoid division by zero: ")
#    print(list_of_probabilities_of_survival_by_alt_functionalization_in_t1)

        
    #calculate s(t1) for dosage
    list_of_probabilities_of_survival_by_dosage = calculate_probability_of_survival_of_duplicate_gene_copy_by_time(b_dos, c_dos, d_dos, f_dos, t1_range)
#    print ("List of prob of surv by dos in t1: ")
#    print(list_of_probabilities_of_survival_by_dosage)
    list_of_probabilities_of_survival_by_dosage_in_t1 = []
    list_of_probabilities_of_survival_by_dosage_in_t1.append(survival_immediately_post_wgd)
    for each_time_point in range(1,q):
        list_of_probabilities_of_survival_by_dosage_in_t1.append(list_of_probabilities_of_survival_by_dosage[each_time_point])    
#    print ("List of prob of surv by dos in t1 adjusted to avoid division by zero: ")
#    print(list_of_probabilities_of_survival_by_dosage_in_t1)

    #calculate s(t1) for nonfunctionalization
    list_of_probabilities_of_survival_of_genes_that_can_only_be_nonfunctionalized = calculate_probability_of_survival_of_duplicate_gene_copy_by_time(b_non, c_non, d_non, f_non, t1_range)
#    print ("List of prob of surv in non category in t1: ")
#    print(list_of_probabilities_of_survival_of_genes_that_can_only_be_nonfunctionalized)    
    list_of_probabilities_of_survival_of_genes_that_can_only_be_nonfunctionalized_in_t1 = []
    list_of_probabilities_of_survival_of_genes_that_can_only_be_nonfunctionalized_in_t1.append(survival_immediately_post_wgd)
    for each_time_point in range(1,q):
        list_of_probabilities_of_survival_of_genes_that_can_only_be_nonfunctionalized_in_t1.append(list_of_probabilities_of_survival_of_genes_that_can_only_be_nonfunctionalized[each_time_point])
#    print ("List of prob of surv in non category in t1 adjusted to avoid division by zero: ")
#    print (list_of_probabilities_of_survival_of_genes_that_can_only_be_nonfunctionalized_in_t1)

    #Establish t2_range
    t2_range = []
    i = 0
    for i in range(q):
        t2_range.append(i/100)
#    print(t2_range)

    #calculate sn(t2) for alt_functionalization
    list_of_probabilities_of_survival_by_alt_functionalization_in_t2 = list_of_probabilities_of_survival_by_alt_functionalization_in_t1
#    print ("List of prob of surv by alt_func in t2: ")
#    print (list_of_probabilities_of_survival_by_alt_functionalization_in_t2)
    #calculate sd(t2) for dosage compensation
    list_of_probabilities_of_survival_by_dosage_in_t2 = list_of_probabilities_of_survival_by_dosage_in_t1
#    print ("List of prob of surv by dos in t2: " )
#    print(list_of_probabilities_of_survival_by_dosage_in_t2)    
    #calculate so(t2) for non-functionalization
    list_of_probabilities_of_survival_of_genes_that_can_only_be_nonfunctionalized_in_t2 = list_of_probabilities_of_survival_of_genes_that_can_only_be_nonfunctionalized_in_t1
#    print ("List of prob of surv in non category in t2: ")
#    print(list_of_probabilities_of_survival_of_genes_that_can_only_be_nonfunctionalized_in_t2) 


    probability_ratio = calculate_pratio(list_of_probabilities_of_survival_by_alt_functionalization_in_t1, list_of_probabilities_of_survival_by_dosage_in_t1, list_of_probabilities_of_survival_of_genes_that_can_only_be_nonfunctionalized_in_t1, list_of_probabilities_of_survival_by_alt_functionalization_in_t2, list_of_probabilities_of_survival_by_dosage_in_t2, list_of_probabilities_of_survival_of_genes_that_can_only_be_nonfunctionalized_in_t2)

    probability_ratio_2d = calculate_pratio_2d(t1_range, t2_range, list_of_probabilities_of_survival_by_alt_functionalization_in_t1, list_of_probabilities_of_survival_by_dosage_in_t1, list_of_probabilities_of_survival_of_genes_that_can_only_be_nonfunctionalized_in_t1, list_of_probabilities_of_survival_by_alt_functionalization_in_t2, list_of_probabilities_of_survival_by_dosage_in_t2, list_of_probabilities_of_survival_of_genes_that_can_only_be_nonfunctionalized_in_t2)

    return t1_range, t2_range, probability_ratio, list_of_probabilities_of_survival_by_alt_functionalization_in_t1, list_of_probabilities_of_survival_by_dosage_in_t1, list_of_probabilities_of_survival_of_genes_that_can_only_be_nonfunctionalized_in_t1, list_of_probabilities_of_survival_by_alt_functionalization_in_t2, list_of_probabilities_of_survival_by_dosage_in_t2, list_of_probabilities_of_survival_of_genes_that_can_only_be_nonfunctionalized_in_t2, probability_ratio_2d
    
#############################################################################
#Main
for each_percentage_combo in range (0, number_percent_combos):
    #Run calculations and extract info from it
    alt_func_percent = alt_func_percentages[each_percentage_combo]
    dos_percent = dos_percentages[each_percentage_combo]
    non_percent = non_percentages[each_percentage_combo]
    
    calc = main_calculations(alt_func_percent, dos_percent, non_percent)
    t1_range = calc[0]
    t2_range = calc[1]
    probability_ratio = calc[2]
    print (probability_ratio)
    list_of_probabilities_of_survival_by_alt_functionalization_in_t1 = calc[3]
    list_of_probabilities_of_survival_by_dosage_in_t1 = calc[4]
    list_of_probabilities_of_survival_of_genes_that_can_only_be_nonfunctionalized_in_t1 = calc[5]
    list_of_probabilities_of_survival_by_alt_functionalization_in_t2 = calc[6]
    list_of_probabilities_of_survival_by_dosage_in_t2 = calc[7]
    list_of_probabilities_of_survival_of_genes_that_can_only_be_nonfunctionalized_in_t2 = calc[8]
    probability_ratio_2d = calc[9]
    print (probability_ratio_2d)

    #plot 3D surface
    # fig = plt.figure()
    # ax = plt.axes(projection='3d')
    # x, y = np.meshgrid(t1_range, t2_range, indexing='ij') 
    # surf = ax.plot_surface(x, y, probability_ratio, cmap=cm.cividis)
    # fig.colorbar(surf)
    # ax.set_title('Pratio for t1 and t2:' + percentages[each_percentage_combo], fontsize=14)
    # ax.set_xlabel('$t1$', fontsize=12)
    # ax.set_ylabel('$t2$', fontsize=12)
    # ax.set_zlabel(r'probability ratio', fontsize=11)
    # ax.view_init(15, 45)
    # plt.draw()
    # plt.pause(.001)
    
    # file = open('array_practice_file_3D_' + percentages_file_name[each_percentage_combo] +'.csv', 'a', newline='')
    # writer = csv.writer(file)
    # i = 0
    # for i in range (0, q):
    #     writer.writerow(probability_ratio[i])
    # file.close()
 

    #plot 3D scatter
    fig2 = plt.figure()
    ax2 = plt.axes(projection='3d')
    surf2 = ax2.scatter3D(probability_ratio_2d[1:,0], probability_ratio_2d[1:,1], probability_ratio_2d[1:,2], c = probability_ratio_2d[1:,2], cmap=cm.cividis)
    fig2.colorbar(surf2)
    ax2.set_title('Pratio for t1 and t2:' + percentages[each_percentage_combo], fontsize=14)
    ax2.set_xlabel('$t1$', fontsize=12)
    ax2.set_ylabel('$t2$', fontsize=12)
    ax2.set_zlabel(r'probability ratio', fontsize=11)
    ax2.view_init(15, 45)
    plt.draw()
    plt.pause(.001)
    
    file2 = open('pratio_array_2D_' + percentages_file_name[each_percentage_combo] +'.csv', 'a', newline='')
    writer2 = csv.writer(file2)
    writer2.writerow('abp')
    i = 1
    for i in range (1, q*q):
        writer2.writerow(probability_ratio_2d[i])
    file2.close()    


#Plot survival curves
fig, t1_plot = plt.subplots()
t1_plot.plot(t1_range, list_of_probabilities_of_survival_by_alt_functionalization_in_t1, color = 'red', label = 'Alt_func')
t1_plot.plot(t1_range, list_of_probabilities_of_survival_by_dosage_in_t1, color = 'blue', label = 'Dos')
t1_plot.plot(t1_range, list_of_probabilities_of_survival_of_genes_that_can_only_be_nonfunctionalized_in_t1, color = 'yellow', label = 'Non')
legend = t1_plot.legend(loc = 'upper right', shadow = True, fontsize = '12')
t1_plot.set_title('Survival over Time', fontsize=14)
t1_plot.set_xlabel('Time Since Duplication Event', fontsize=12)
t1_plot.set_ylabel('Proportion Gene Duplicate Copies Surviving', fontsize=12)
plt.show()