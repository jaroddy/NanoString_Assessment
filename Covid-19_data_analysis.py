# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 12:50:08 2021

@author: jarod
"""

#import Packages

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA



def excel_loader(excel_name):
    #Split Excel Sheet into 2 dataframes
    
    annotations = pd.read_excel(excel_name, sheet_name ="Annotations")

    normalized_counts = pd.read_excel(excel_name, sheet_name = "NormalizedCounts")

    return annotations, normalized_counts

def ROI_per_patient(annotations):
    #Count Instances of ROI per patient
    
    patient_list = annotations.patient.unique()
    
    ROIs = pd.DataFrame(index = range(1), columns = patient_list)
    
    counter = 0
    
    for patients in range(len(patient_list)):
        
        temp_patient = patient_list[patients]
        
        for data in range(len(annotations.patient)):
            
            if temp_patient == annotations.patient[data]:
                
                counter = counter + 1
                
        ROIs[temp_patient][0] = counter
        
        
    return ROIs, patient_list
        

def seqsat_graph(annotations, patient_list):
    
    fig = plt.figure(figsize = (14,10), dpi = 300)
    
    ax = fig.add_subplot(111, projection='3d')
    
    annotations_num = annotations
    
    #Changing 3 variables to numerical form for graphing
    
    num_pl = range(len(patient_list))
    
    segtype_list = annotations['Segment Type'].unique()
    
    num_sg = range(len(segtype_list))
    
    tissue_list = annotations['Tissue Substructure'].unique()
    
    num_ts = range(len(tissue_list))
    
        
    for x in range(len(patient_list)):
        
        annotations_num.patient = annotations_num.patient.replace(patient_list[x],num_pl[x])
        
    for y in range(len(segtype_list)):
        
        annotations_num['Segment Type'] = annotations_num['Segment Type'].replace(segtype_list[y],num_sg[y])
        
    for z in range(len(tissue_list)):
        
        annotations_num['Tissue Substructure'] = annotations_num['Tissue Substructure'].replace(tissue_list[z],num_ts[z])

    #Plot Creation
        
    img = ax.scatter(annotations_num.patient, 
                  annotations_num['SequencingSaturation'], 
                  annotations_num['Segment Type'], 
                  c=annotations_num['Tissue Substructure'], 
                  cmap=plt.hot())

    #Setting Axis labels, tickmarks and title
    
    ax.set_title("Sequencing Saturation by Patient ID, Segment Type and Tissue Substructure", fontsize=18)
    
    plt.xticks(num_pl,patient_list, fontsize=10)
    plt.xlabel("Patient ID")
    
    plt.yticks(range(0,100,10), fontsize=10)
    plt.ylabel("Sequencing Saturation")
    
    ax.set_zticks(num_sg)
    ax.set_zticklabels(segtype_list, fontsize=10, rotation=30)
    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel("Segment Type", rotation=90, labelpad=15)
    
    cbar = fig.colorbar(img, orientation='vertical', ticks=num_ts)
    cbar.ax.set_yticklabels(tissue_list)
    cbar.ax.set_ylabel("Tissue Substructure")
    
    ax.view_init(30,300)
    
    #Saving the figure
    
    plt.savefig("Sequencing Saturation by Patient ID, Segment Type and Tissue Substructure.png")
    
def dim_reduction(normalized_counts):
    
    #Preparing the counts by log2 transformation
    
    counts_log2 = normalized_counts
    
    counts_log2 = normalized_counts.drop(['TargetName'],axis=1)
    
    counts_log2 = counts_log2.set_index(normalized_counts.TargetName)
    
    counts_log2 = np.log2(counts_log2)
    
    counts_log2 = counts_log2.transpose()
    
    sample_name = counts_log2.index.tolist()
    
    #Actual PCA Steps
    
    PCA_95 = PCA(n_components=30)
    
    PCA_95.fit(counts_log2)
    
    PCA_counts = PCA_95.transform(counts_log2)
    
    variance = PCA_95.explained_variance_ratio_
    
    return PCA_counts, variance, sample_name
    
def pca_graph(PCA_counts, annotations, patient_list, normalized_counts):
    
    #Creating 4 graphs based on 3 most relevant principal components
    
    graph_comp = pd.DataFrame(index=range(len(annotations)),columns=['SARS-CoV-2 RNA ISH', 'Patient ID', 'Tissue Substructure', 'Expression ORF1ab'])
    
    fig = plt.figure(figsize = (10,10), dpi = 300)
    
    ax1 = fig.add_subplot(111, projection='3d')
    
    annotations_num = annotations
    
    #Converting RNA ISH values to numerical form
    
    rnaish_list = annotations['SARS-CoV-2 RNA ISH'].unique()
    
    num_rnaish = range(len(rnaish_list))
    
    for x in range(len(rnaish_list)):
        
        annotations_num['SARS-CoV-2 RNA ISH'] = annotations_num['SARS-CoV-2 RNA ISH'].replace(rnaish_list,num_rnaish[x])
    
    #Grabbing the ORF1ab Expression information
    
    normalized_counts = normalized_counts.transpose()
    
    normalized_counts.columns = normalized_counts.iloc[0]
    
    normalized_counts = normalized_counts.drop(['TargetName'])
    
    index_nc = normalized_counts.index
    
    orf1ab_exp = normalized_counts["ORF1ab"]
    
    #Filling out variables for easy graphing
     
    for y in range(len(graph_comp)):
        
        for z in range(len(normalized_counts)):
            
            if index_nc[z] == annotations['Sample Name'][y]:
        
                graph_comp['SARS-CoV-2 RNA ISH'][y] = annotations.patient[y]
                
                graph_comp['Patient ID'][y] = annotations.patient[y]
                
                graph_comp['Tissue Substructure'][y] = annotations['Tissue Substructure'][y]
                
                graph_comp['Expression ORF1ab'][y] = orf1ab_exp[y]
   
    #Setting Axis labels, tickmarks and title for 4 different graphs
    #4 full size graphs were made with the colorbar changing to visualize expression by different metrics
    
    img = ax1.scatter(PCA_counts[:,0],
                  PCA_counts[:,1] , 
                  PCA_counts[:,2] , 
                  c=graph_comp['Expression ORF1ab'], 
                  cmap=plt.hot())
    
    ax1.set_title("Gene Expression by 3 Principal Components and ORF1ab", fontsize=16)

    plt.xlabel("Principal Component 1")
    
    plt.ylabel("Principal Component 2")
    
    ax1.zaxis.set_rotate_label(False)
    ax1.set_zlabel("Principal Component 3", rotation=90, labelpad=15)
    
    cbar = fig.colorbar(img, orientation='vertical')
    
    cbar.ax.set_ylabel("ORF1ab Expression")
    
    ax1.view_init(30,300)
    
    plt.savefig("Gene Expression by 3 Principal Components and ORF1ab.png")
    
    img = ax1.scatter(PCA_counts[:,0],
                  PCA_counts[:,1] , 
                  PCA_counts[:,2] , 
                  c=graph_comp['Patient ID'], 
                  cmap=plt.hot())
    
    cbar.ax.set_ylabel("Patient ID")
    ax1.set_title("Gene Expression by 3 Principal Components and Patient ID", fontsize=16)
    plt.savefig("Gene Expression by 3 Principal Components and Patient ID.png")   
    
    img = ax1.scatter(PCA_counts[:,0],
                  PCA_counts[:,1] , 
                  PCA_counts[:,2] , 
                  c=graph_comp['SARS-CoV-2 RNA ISH'], 
                  cmap=plt.hot())
    
    cbar.ax.set_ylabel("SARS-CoV-2 RNA ISH")
    ax1.set_title("Gene Expression by 3 Principal Components and RNA ISH", fontsize=16)
    plt.savefig("Gene Expression by 3 Principal Components and RNA ISH.png")  
    
    img = ax1.scatter(PCA_counts[:,0],
                  PCA_counts[:,1] , 
                  PCA_counts[:,2] , 
                  c=graph_comp['Tissue Substructure'], 
                  cmap=plt.hot())
    
    cbar.ax.set_ylabel("Tissue Substructure")
    ax1.set_title("Gene Expression by 3 Principal Components and Tissue Substructure", fontsize=16)
    plt.savefig("Gene Expression by 3 Principal Components and Tissue Substructure.png")  

    return normalized_counts, graph_comp
    

#Defining Global Variables

excel_name = "COVID-19_data.xlsx"


#Begin script

annotations, normalized_counts = excel_loader(excel_name)

ROIs, patient_list = ROI_per_patient(annotations)

seqsat_graph(annotations, patient_list)

PCA_counts, variance, sample_name = dim_reduction(normalized_counts)

pca_graph(PCA_counts, annotations, patient_list, normalized_counts)

ROIs.to_csv('ROIs per Patient.csv', index=False)









