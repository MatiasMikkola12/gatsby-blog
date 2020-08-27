---
title: "Combinatorial drug dose-response data analysis"
date: "2020-08-22T22:12:03.284Z"
description: "Step-by-step tutorial on how to analyze combinatorial drug dose-response data and derive measures of combinatorial drug efficacy."
---


Even though there are many anticancer drugs available in the market, cancer is still one of the deadliest diseases in the world. One of the reasons for this is that cancer is a complex disease, and it affects several redundant pathways or mechanisms in a cell to achieve the same result (over-production of cells). However, most anticancer drugs are designed to target only one or a few cancer-causing gene mutations. Because of this, even if a drug is designed to target a highly relevant gene(s), mutations in other (redundant) genes can reactivate the growth of cancer cells (this causes what sometimes is referred to as a ‘drug-resistant’ cancer). 

One way to address this is to provide drug combinations, with each drug targeting a different but redundant cancer-causing mechanism, for example, different genes involved in the same biological mechanism. So, how do we know if a drug combination works better than each drug individually? Because we’re talking about therapies given to patients, before being used as therapies, each drug combination would need to be tested through clinical trials. But clinical trials are expensive and can take several years, so there’s no way we can have clinical trials for all the hundreds of potential anticancer drug combinations. To prioritize the most promising drugs and drug combinations, researchers test drug combinations’ effects in cell cultures with dose-response studies.

In single drug dose-response studies, we have cancerous cell cultures to which a drug is administered in increasing concentration. For each of drug concentrations, the inhibition of cell growth is measured. The higher the inhibition of cell measured, the more effective the drug is. Combinatorial drug dose-response studies are similar to the single drug studies. The only difference is that 2 drugs (or N drugs) are administered in varying concentrations, instead of just one drug. 

We can observe the difference in output of experiments for single and combinatorial studies in Figure 1. Here, we have drugs Navitoclax and AT-406. In single drug experiments, we have 2D plots for each drug of the cell inhibition vs the drug concentration. In combinatorial experiments, we have a heatmap showing the cell inhibition vs concentrations of both drugs. From the plots, several measures of single and combinatorial drug efficacy can be obtained. For individual drugs, probably the most widely used measure of drug efficacy is IC50, which is the concentration of the drug that is required to achieve a 50% inhibition of cell growth. The lower the IC50, the more potent the drug is. For drug combinations, the efficacy is measured in comparison with the individual drugs using a measure of synergism, which tests whether 2 drugs in combination produce a greater effect than each one individually. 


![alt text](images/figure_1.jpg "Figure 1")

<br>
EXPLAIN HOW SYNERGY IS MEASURED IN MOST POPULAR MODELS (LOEWE, HSA, BLISS, ZIP; WITH EQUATIONS)
STEP-BY-STEP DATA ANALYSIS IN R