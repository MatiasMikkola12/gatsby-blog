---
title: "Combinatorial drug dose-response data analysis"
date: "2020-08-22T22:12:03.284Z"
description: "Step-by-step tutorial on how to analyze combinatorial drug dose-response data and derive measures of combinatorial drug efficacy."
---

This is part II of analysis of drug dose-response data, where we look at *combinatorial* drug dose-response experiments. For *single* drug dose-response data analysis and introduction of basic concepts on measures of drug efficacy, go to my post on [drug dose-response data analysis](../post-single-dose-response-data-analysis).
 
Even though there are many anticancer drugs available in the market, cancer is still one of the deadliest diseases in the world. One of the reasons for this is that cancer is a complex disease, affecting several redundant pathways or mechanisms in a cell to achieve the same result (over-production of cells). However, most anticancer drugs are designed to target only one or a few cancer-causing gene mutations. Because of this, even if a drug is designed to target a highly relevant gene(s), mutations in other (redundant) genes can reactivate the growth of cancer cells (this leads to what sometimes is referred to as a ‘drug-resistant’ cancer). 

One way to address this is to provide drug combinations, with each drug targeting a different but redundant cancer-causing mechanism, for example, different genes involved in the same biological mechanism. So, how do we know if a drug combination works better than each drug individually? Because we’re talking about drugs given to patients, before being used as therapies, each drug combination would need to be tested through clinical trials. But clinical trials are expensive and can take several years, so it's not possible to have clinical trials for all the hundreds of potential anticancer drug combinations. To prioritize the most promising drugs and drug combinations, researchers test drug combinations’ effects in cell cultures with dose-response studies.

In *single* drug dose-response studies, we have cancerous cell cultures to which a drug is administered in increasing concentrations. For each of the drug concentrations, the inhibition of cell growth is measured. The higher the inhibition of cell growth measured, the more effective the drug is. Combinatorial drug dose-response studies are similar to the single drug studies. The only difference is that 2 drugs (or N drugs) are administered in varying concentrations, instead of just one drug. 

We can observe the difference in output of experiments for single and combinatorial studies in the next figure:
![alt text](images/figure_1.jpg "Figure 1")

Here, we have drugs Navitoclax and AT-406. Single drug experiments produce a measure of cell growth inhibition for each concentration of the drug that is administered. Combinatorial experiments produce a measure of cell growth inhibition for each combination of concentrations of the drugs being administered. 

Efficacy of a drug combination is measured by comparing the measured effect of the drug combination vs the expected effect of the drugs if they did not interact with each other. Let's assume *y1* is the effect of one drug at dose *x1*, *y2* is the effect of the other drug at dose *x2*, and *yc* is the effect produced by combining them. A drug combination can be classified as synergistic, antagonistic or non-interactive, depending on whether yc is greater, less or equal to the expected effect under the assumption of non-interaction, respectively. To calculate the expected effect, we can use a variety of reference models, which use different assumptions to define the expected effect under the hypothesis of non-interaction. There are four frequently-used reference models: 

- In the **HSA Model** [1] the expected combination effect equals to the higher effect of individual drugs (*y<sub>HSA</sub> = max(y1, y2)* ). Any additional effect over the higher single drug will be considered as an HSA synergy
- The **Loewe Model** [2] states that the expected effect must satisfy
<p align="center">
<img width="180" src="https://render.githubusercontent.com/render/math?math=\frac{x_1}{X_{LOEWE}^1} %2B \frac{x_2}{X_{LOEWE}^2} = 1"/>
</p>

where *X <sup>1</sup><sub>LOEWE</sub>*  and *X <sup>2</sup><sub>LOEWE</sub>* are the doses of drug 1 and 2 alone that produce the expected effect *y<sub>LOEWE</sub>*. If we suppose that the observed combination effect *yc* *>* *y<sub>LOEWE</sub>* then we have *X <sup>1</sup><sub>c</sub>* *>* *X <sup>1</sup><sub>LOEWE</sub>* and *X<sup>2</sup><sub>c</sub>* *>* *X<sup>2</sup><sub>LOEWE</sub>*, where *X<sup>1</sup><sub>c</sub>* and *X<sup>2</sup><sub>c</sub>* are defined similarly as the doses of drug 1 and 2 alone that produce *yc*. This leads to:

<p align="center">
<img width="120" src="https://render.githubusercontent.com/render/math?math=\frac{x_1}{X_{c}^1} %2B \frac{x_2}{X_{c}^2} < 1"/>
</p>


The left part of this equation is also called Combination Index (CI), with CI < 1 corresponding to *yc* *>* *y<sub>LOEWE</sub>* and thus a Loewe synergy.
                 
                 
- The **Bliss Model** [3] assumes a stochastic process in which two drugs elicit their effects independently, and the expected combination effect can be calculated based on the probability of independent events as *y<sub>BLISS</sub>* *=* *y1+y2-y1y2*

- The **Zero Interaction Potency Model (ZIP)** [4] combines the advantages of the BLISS and Loewe model. It assumes that two non-interacting drugs are expected to experience minimal changes in their dose-response curves. It considers the dose-response curves of individual drugs and of their combination to derive a score to quantify their interactions. Check [4] for complete mathematical derivation of the ZIP synergy score. 

<br>

Now that we know how to measure efficacy of drug combinations, let's go to the data analysis in R step by step.

### Step-by-step data analysis

We will analyze drug combination data taken from [DrugComb](https://drugcomb.fimm.fi) which you can download <a href="toy_drug_combination_data.csv" download>here</a>.

First, we load libraries we will use and load the data:
```
library(drc)
library(ggplot2)
library(magrittr)
library(synergyfinder)
library(Cairo)

data = data.frame(read.csv('toy_drug_combination_data.csv'))
head(data)
```
![alt text](images/head_data_.png "Figure 2")
<br>

The data matrix contains measures of cell growth inhibition at varying concentrations of drug 5-FU (`conc_r`) and drug BEX-235 (`conc_c`). We can fit dose-response curves to individual drugs by subselecting the rows in which `conc_r` or `conc_c` are equal to zero and using the `drm` model from the `drc` package (for more details, check my post on [drug dose-response data analysis](../post-single-dose-response-data-analysis)):

```
##Drug 5-FU
data_row = data[data$conc_c == 0,]
fitted_curve <- drm(formula = inhibition ~ conc_r,
             data = data_row,
             fct = LL.4())
plot(fitted_curve, 
     log='x', 
     xlab = 'Drug concentration (microM)', 
     ylab= 'Inhibition')
```

![alt text](images/image_fitted_dose_response_5fu.png "DR fitted 5-FU" )

```
##Drug BEZ-235
data_col = data[data$conc_r == 0,]
fitted_curve <- drm(formula = inhibition ~ conc_c,
                    data = data_col,
                    fct = LL.4())
plot(fitted_curve, 
     log='x', 
     xlab = 'Drug concentration (microM)', 
     ylab= 'Inhibition')
```

![alt text](images/image_fitted_dose_response_bez.png "DR fitted 5-FU" )

Then we measure ZIP drug synergism with the package `synergyfinder` (other synergy reference models are available). The overall synergy score is the average score of drug synergy across all concentrations measured:
```
names(data)[names(data) == "inhibition"] <- "response"

data <- ReshapeData(
  data,
  impute = TRUE,
  noise = TRUE, correction = "non", data.type = "inhibition"
)


ZIP <- CalculateSynergy(data, method = "ZIP", adjusted = TRUE)
synergy_score <- mean(ZIP$scores$'10'[ZIP$scores$'10' != 0])
```
![lastimage](images/synergy_score_.png "Figure synergy score")

<br>
And our analysis is done!





<br>
<br>
<br>
<br>
<br>

##### References
[1] M.C. Berenbaum. What is synergy? Pharmacol Rev, 41 (1989), pp. 93-141

[2] S. Loewe. The problem of synergism and antagonism of combined drugs. Arzneimittelforschung, 3 (1953), pp. 285-290

[3] C.I. Bliss. The toxicity of poisons applied jointly. Ann Appl Biol, 26 (1939), pp. 585-615, 10.1111/j.1744-7348.1939.tb06990.x

[4] Yadav, B., Wennerberg, K., Aittokallio, T., & Tang, J. Searching for Drug Synergy in Complex Dose-Response Landscapes Using an Interaction Potency Model. Computational and structural biotechnology journal, 13 (2015), pp. 504-513. 10.1016/j.csbj.2015.09.001






