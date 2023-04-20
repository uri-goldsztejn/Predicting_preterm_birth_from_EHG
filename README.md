# Predicting preterm births from electrohysterogram recordings via deep learning

This is the code developed for our paper "[Predicting preterm births from electrohysterogram recordings via deep learning](https://www.medrxiv.org/content/10.1101/2022.12.25.22283937v1)"  

Uri Goldsztejn, Arye Nehorai  
Washington University in St. Louis, 2023

*If you find this code useful, please consider citing.*

## Content
* Overview
* Files
* Requirements
* Citation

## Overview

About one in ten babies is born preterm, i.e., before completing 37 weeks of gestation, which can result in permanent neurologic deficit and is a leading cause of child mortality. However, identifying mothers at high risk of preterm labor can prompt beneficial treatments to reduce the incidence of preterm births and improve their outcomes.

Here, we introduce a deep learning approach to predict preterm births directly from electrohysterogram (EHG) measurements of pregnant mothers recorded at around 31 weeks of gestation. Our model, which includes a recurrent neural network, predicts preterm births using short-time Fourier transforms of EHG recordings and clinical information from two public datasets. 

Our findings suggest that preterm births can be predicted from EHGs and clinical information around the 31st week of gestation with comparable accuracy as current clinical practices which predict labor within one week. Our method holds promise for reducing newborn morbidity and mortality, especially in populations with limited access to healthcare, who suffer more from preterm birth.



## Files

* *metadata_classification.m* - Implements the classification model based on clinical information.
* *metadata_regression.m* - Implements the regression model based on clinical information.
* *EHG_classification.m* - Implements the classification model based on EHG measurements.
* *EHG_regression.m* - Implements the regression model based on clinical measurements.
* *combined_classification.m* - Implements the classification model based on EHG measurements and clinical information.
* *combined_regression.m* - Implements the regression model based on clinical measurements and clinical information.
* *results/* - Scripts for recreating the figures in the publication.
* *occlusions/* - Scripts for recreating the ablation study described in the publication.
* *utils/* - Auxiliary functions.
* *dataset/* - EHG and clinical information from the Physionet repository.

## Data

As described in our paper, we combined data from the TPEHG DB and the TPEHGT DS.
The data from the TPEHG DB is available at https://physionet.org/content/tpehgdb/1.0.1/, doi: https://doi.org/10.13026/C2FW2V.
The data from the TPEHGT DS is available at https://physionet.org/content/tpehgt/1.0.0/, doi: https://doi.org/10.13026/C2166R.

## Software requirements

MATLAB>=R2021a 

## Citation

Goldsztejn, Uri, and Arye Nehorai. "Predicting preterm births from electrohysterogram recordings via deep learning." medRxiv (2022): 2022-12.
