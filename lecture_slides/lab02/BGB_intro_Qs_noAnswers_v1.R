


Q. What is Psychotria?

An interesting psychadelic drug

The name of a horror movie

A large genus of shrubs (1500 species) in the coffee family (Rubiaceae), many of them found in the Pacific.

A popular baby name in 2020




Q. How many "species" (at least they are OTUs: Operational Taxonomic Units) are in the Psychotria phylogeny we are using?

18
19
20
21



Q. About how old is the clade we are studying (Hawaiian Psychotria), according to the phylogeny we are using (taken from Ree & Smith 2008)?

5.2 years old

5,200 years old

5,200,000 years old

5,200,000,000 years old




Q. In what order did Hawaii's "high islands" erupt from the ocean?

Kauai, Oahu, Maui Nui, Hawaii Big Island
Oahu, Kauai, Maui Nui, Hawaii Big Island
Hawaii Big Island, Maui Nui, Oahu, Kauai
Kauai, Hawaii Big Island, Maui Nui, Oahu 



Q. Look at the "tipranges" object, which stores the geographic ranges of the tip species. How many OTUs inhabit Kauai?

1
4
5
8


Q. Look at the "tipranges" object, which stores the geographic ranges of the tip species. How many OTUs inhabit Hawaii?

1
4
5
8




Q. Use the "numstates_from_numareas" function. If you have a 2-area problem, and if the maximum number of areas per range is 2, and if the null range is an allowed state, how many states (=possible geographic ranges) do you have in the state space?

2
4
8
16
1024
1,048,576



Q. Use the "numstates_from_numareas" function. If you have a 3-area problem, and if the maximum number of areas per range is 3, and if the null range is an allowed state, how many states (=possible geographic ranges) do you have in the state space?

2
4
8
16
1024
1,048,576


Q. Use the "numstates_from_numareas" function. If you have a 4-area problem, and if the maximum number of areas per range is 4, and if the null range is an allowed state, how many states (=possible geographic ranges) do you have in the state space?

2
4
8
16
1024
1,048,576



Q. Use the "numstates_from_numareas" function. If you have a 4-area problem, and if the maximum number of areas per range is 2, and if the null range is an allowed state, how many states (=possible geographic ranges) do you have in the state space?

2
4
10
11
20
1024



Q. Use the "numstates_from_numareas" function. If you have a 10-area problem, and if the maximum number of areas per range is 10, and if the null range is an allowed state, how many states (=possible geographic ranges) do you have in the state space?

2
4
8
16
21
211
1024
1,048,576


Q. Use the "numstates_from_numareas" function. If you have a 20-area problem, and if the maximum number of areas per range is 20, and if the null range is an allowed state, how many states (=possible geographic ranges) do you have in the state space?

2
4
8
16
21
211
1024
1,048,576



Q. Use the "numstates_from_numareas" function. If you have a 20-area problem, and if the maximum number of areas per range is 2, and if the null range is an allowed state, how many states (=possible geographic ranges) do you have in the state space?

2
4
8
16
21
211
1024
1,048,576


Q. Use the "numstates_from_numareas" function. If you have a 20-area problem, and if the maximum number of areas per range is 1, and if the null range is an allowed state, how many states (=possible geographic ranges) do you have in the state space?

2
4
8
16
21
211
1024
1,048,576


Q. What is a likelihood, in statistics?

The probability of a model, conditional on the data

The probability of the data, conditional on a model

An informal term that just means probability, just like usage in the general public

The prior probability of the model

The posterior probability of the model




Q. What is log-likelihood?


The log10 of the likelihood, used to help represent very small likelihoods that are hard for humans to read and hard for computers to store in memory.

The natural logarithm of the likelihood, used to help represent very small likelihoods that are hard for humans to read and hard for computers to store in memory.

A confusing math thing we shouldn't bother understanding, because math is hard.

A logbook of stored likelihoods







Q. Run the DEC and DEC+J models.


Approximately what log-likelihood does the DEC model settle on during the Maximum Likelihood search?

-34.5
-20.9
20.9
34.5



Q. Approximately what log-likelihood does the DEC+J model settle on during the Maximum Likelihood search?

-34.5
-20.9
20.9
34.5





Q. What is the most-probable ancestral range (at the root node) of Hawaiian Psychotria, under the DEC model?
K
O
KO
KOM
KOMH


Q. What is the approximate probability of the most-probable ancestral range (at the root node) of Hawaiian Psychotria, under the DEC model?
25%
49%
55%
100%


Q. What is the most-probable ancestral range (at the root node) of Hawaiian Psychotria, under the DEC+J model?
K
O
KO
KOM
KOMH

Q. What is the approximate probability of the most-probable ancestral range (at the root node) of Hawaiian Psychotria, under the DEC+J model?
25%
49%
55%
100%


Q. What is the displayed most-probable ancestral range for the common ancestor of P_hawaiiensis_Makaopuhi and P_mariniana_MauiNui, under the DEC model?

M
H
MH
OMH

Q. What is the approximate probability of the most-probable ancestral range for the common ancestor of P_hawaiiensis_Makaopuhi and P_mariniana_MauiNui, under the DEC model?

About 10%
About 50%
About 90%
About 100%


Q. What is the displayed most-probable ancestral range for the common ancestor of P_hawaiiensis_Makaopuhi and P_mariniana_MauiNui, under the DEC+J model?

M
H
MH
OMH

Q. What is the approximate probability of the most-probable ancestral range for the common ancestor of P_hawaiiensis_Makaopuhi and P_mariniana_MauiNui, under the DEC+J model?

About 10%
About 50%
About 90%
About 100%









Q. What lnL is higher, -20.9 or -34.5?

-34.5
-20.9


Q. The opposite of the natural logarithm is the exp() function.  It takes Euler's number (2.718282, a fundamental mathematical constant) to the power of the input.  This converts log-likelihoods to regular likelihoods (conditional probabilities).

Calculate exp(-20.9). This will give you the maximized regular likelihood under the DEC+J model.
Calculate exp(-34.5). This will give you the maximized regular likelihood under the DEC model.


What is the ratio of these two probabilities?  I.e., how much more probable is the data under the DEC+J model, compared to the DEC model?

About 2.0
About 8.38
About 13.6
About 806,130





Q. After running all six models, some summary tables are produced comparing the models. Look at:

conditional_format_table(restable_AICc_rellike)



Q. What is the ML estimate of the "d" parameter (range-expansion dispersal rate) under the DEC model?

About 0.0
About 0.01
About 0.028
About 0.034
About 0.11

Q. What is the ML estimate of the "e" parameter (range-contraction rate, i.e. extirpation or single-area local extinction) under the DEC model?

About 0.0
About 0.01
About 0.028
About 0.034
About 0.11


Q. What is the ML estimate of the "d" parameter (range-expansion dispersal rate) under the DEC+J model?

About 0.0
About 0.01
About 0.028
About 0.034
About 0.11

Q. What is the ML estimate of the "e" parameter (range-contraction rate, i.e. extirpation or single-area local extinction) under the DEC+J model?

About 0.0
About 0.01
About 0.028
About 0.034
About 0.11

Q. What is the ML estimate of the "j" parameter (the per-event relative weight of jump dispersal at cladogenesis) under the DEC+J model?

About 0.0
About 0.01
About 0.028
About 0.034
About 0.11




Q. Akaike Information Criterion (AIC) is a form of "penalized likelihood", where the log-likelihood is penalized based on model complexity (meaning, the number of free parameters estimated from the data. Models with more free parameters are more penalized). AICc is sample-size corrected AIC, which increases the penalty on model complexity for small datasets.

AIC or AICc values can be transformed to model weights, which assign all compared models a fractional weight. The weights add up to 1.0, and can be thought of as percentages of support.

Looking at conditional_format_table(restable_AICc_rellike), which model is the best fit, having the highest AICc model weight?

DEC
DEC+J
DIVALIKE
DIVALIKE+J
BAYAREALIKE
BAYAREALIKE+J



Q. Which other models are also highly credible (in the top 95% of model weight)?

DEC & DIVALIKE
DIVALIKE & DIVALIKE+J
BAYAREALIKE & BAYAREALIKE+J
DIVALIKE+J & BAYAREALIKE+J


Q. Adding up the model weight of the 3 "+J" models, how muh of the total model weight do they garner?

About 0%
About 42%
About 99%
About 100%


Q. Adding up the model weight of the 3 non-"+J" models, how muh of the total model weight do they garner?

About 0%
About 42%
About 99%
About 100%







