#######################################################
# Level the tips of the example Psychotria tree
#
# The Psychotria tree has slight variations in tip
# heights that mean they don't all come exactly to
# 0.0.
#
# This script fixes that.
#######################################################
library(BioGeoBEARS)


# Example Psychotria tree
trstr_orig = "((((((((P_hawaiiensis_WaikamoiL1:0.9656850499,P_mauiensis_Eke:0.9656850499):0.7086257935,(P_fauriei2:1.230218511,P_hathewayi_1:1.230218511):0.4440923324):0.1767115552,(P_kaduana_PuuKukuiAS:1.851022399,P_mauiensis_PepeAS:1.851022399):0.0008897862802):0.3347375986,P_kaduana_HawaiiLoa:2.185759997):0.302349378,(P_greenwelliae07:1.131363255,P_greenwelliae907:1.131363255):1.35674612):1.689170274,((((P_mariniana_MauiNui:1.994011054,P_hawaiiensis_Makaopuhi:1.994011054):0.7328279804,P_mariniana_Oahu:2.726839034):0.2574151709,P_mariniana_Kokee2:2.984254205):0.4601084855,P_wawraeDL7428:3.444362691):0.732916959):0.7345185743,(P_grandiflora_Kal2:2.479300491,P_hobdyi_Kuia:2.479300491):2.432497733):0.2873119899,((P_hexandra_K1:2.363984189,P_hexandra_M:2.363984189):0.4630447802,P_hexandra_Oahu:2.826939991):2.372081244);"

tr = read.tree(text=trstr_orig)

tr2 = level_tree_tips(tr, method="highest", printflag=FALSE, fossils_older_than=0.01)

trstr = write.tree(tr2, file="")
trstr

"((((((((P_hawaiiensis_WaikamoiL1:0.9665748366,P_mauiensis_Eke:0.9665748366):0.7086257935,(P_fauriei2:1.231108298,P_hathewayi_1:1.231108298):0.4440923324):0.1767115552,(P_kaduana_PuuKukuiAS:1.851022399,P_mauiensis_PepeAS:1.851022399):0.0008897862802):0.3347375986,P_kaduana_HawaiiLoa:2.186649784):0.302349378,(P_greenwelliae07:1.132253042,P_greenwelliae907:1.132253042):1.35674612):1.689170274,((((P_mariniana_MauiNui:1.99490084,P_hawaiiensis_Makaopuhi:1.99490084):0.7328279804,P_mariniana_Oahu:2.72772882):0.2574151709,P_mariniana_Kokee2:2.985143991):0.4601084855,P_wawraeDL7428:3.445252477):0.732916959):0.7345185743,(P_grandiflora_Kal2:2.480190277,P_hobdyi_Kuia:2.480190277):2.432497733):0.2873119899,((P_hexandra_K1:2.364873976,P_hexandra_M:2.364873976):0.4630447802,P_hexandra_Oahu:2.827918756):2.372081244);"
