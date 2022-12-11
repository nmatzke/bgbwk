# How to query BSM event tables

# How many tables / BSMs do you have?
length(clado_events_tables)
length(ana_events_tables)

# The structure of the events tables is the same as the prt(tr) table structure
trtable = prt(tr)
# double brackets [[1]] gives you the first list item
clado_events_table = clado_events_tables[[1]]

# Check dimensions
dim(trtable)
dim(clado_events_table)

# they both have 37 rows because
# the Psychotria tree has 19 tip nodes and 18 internal nodes
# 19 + 18 = 37

# To see what information is stored in the table,
# use names to get column names
names(trtable)
names(clado_events_table)


# To see what node numbers programs use:
# http://phylo.wikidot.com/example-biogeobears-scripts#lagrange_DIVA_node_numbers

###################################################
# Plot APE/BioGeoBEARS node numbers
###################################################
ntips = length(tr$tip.label)
Rnodenums = (ntips+1):(ntips+tr$Nnode)
tipnums = 1:ntips
plot(tr, label.offset=0.25, cex=1.25)
axisPhylo()
tiplabels(cex=1.5)
nodelabels(text=Rnodenums, node=Rnodenums, cex=1.5)
title("APE/BioGeoBEARS node numbers")

# How to get the list of states:
# Get your list of ranges/states, from your list of areas
# http://phylo.wikidot.com/example-biogeobears-scripts#list_of_ranges
# Get your states list (assuming, say, 4-area analysis, with max. rangesize=4)
max_range_size = 4
#areas = getareas_from_tipranges_object(tipranges)
areas = c("K", "O", "M", "H")

# This is the list of states/ranges, where each state/range
# is a list of areas, counting from 0
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)

# Make the list of ranges
ranges_list = NULL
for (i in 1:length(states_list_0based))
{    
  if ( (length(states_list_0based[[i]]) == 1) && (is.na(states_list_0based[[i]])) )
  {
    tmprange = "_"
  } else {
    tmprange = paste(areas[states_list_0based[[i]]+1], collapse="")
  }
  ranges_list = c(ranges_list, tmprange)
}

# Look at the ranges list
ranges_list

# So e.g. the 2nd state is the 2nd range in 
# ranges_list:
ranges_list[2]
# "K"



# To see what information is stored in the table,
# use names to get column names
names(trtable)
names(clado_events_table)

# Look at the BSM information for some of the last nodes
tail(clado_events_table[,13:21])

# For a particular node / row in the trtable:
# 
# sampled_states_AT_nodes - BSM sampled state AT the node
# sampled_states_AT_brbots - BSM state at branch bottom below the node (i.e.the corner below)
# left_desc_nodes - node number of the left descendant node
# right_desc_nodes - node number of the right descendant node
# samp_LEFT_dcorner - BSM sampled state for that left desc. corner
# samp_RIGHT_dcorner - BSM sampled state for that right desc. corner
# event type, txt=description, dispersal from/to:
#clado_event_type
#clado_event_txt
#clado_dispersal_to
#anagenetic_events_txt_below_node - along the branch, IF there are anagenetic events, they are listed as text here

# Example query - dates and destinations of founder events
founder_event_destination = NULL
founder_event_date = NULL

for (i in 1:50)
{
  clado_events_table = clado_events_tables[[i]]
  TF = clado_events_table$clado_event_type == "founder (j)"
  row_nums = (1:nrow(clado_events_table))[TF]
  for (j in 1:length(row_nums))
  {
    tmprow = row_nums[j]
    tmp_dest = clado_events_table$clado_dispersal_to[tmprow]
    tmp_time = clado_events_table$time_bp[tmprow]
    founder_event_destination = c(founder_event_destination, tmp_dest)
    founder_event_date = c(founder_event_date, tmp_time)
      }
}

# Histogram for each jump dispersal destination
TF = founder_event_destination == "K"
hist(founder_event_date[TF], main="times of jumps to K", xlim=c(5.2,0))
mean(founder_event_date[TF])
1.96*sd(founder_event_date[TF])
# Mean and 95% CI on date of jumps to K: 3.14 +/- 1.63 

TF = founder_event_destination == "O"
hist(founder_event_date[TF], main="times of jumps to O", xlim=c(5.2,0))
mean(founder_event_date[TF])
1.96*sd(founder_event_date[TF])

TF = founder_event_destination == "M"
hist(founder_event_date[TF], main="times of jumps to M", xlim=c(5.2,0))
mean(founder_event_date[TF])
1.96*sd(founder_event_date[TF])

TF = founder_event_destination == "H"
hist(founder_event_date[TF], main="times of jumps to H", xlim=c(5.2,0))
mean(founder_event_date[TF])
1.96*sd(founder_event_date[TF])


# To load a saved clado_events_tables:

# Loads to: RES_clado_events_tables
load(file="RES_clado_events_tables.Rdata")

# How to convert multi-page PDFs to animations
# http://phylo.wikidot.com/biogeographical-stochastic-mapping-example-script#PDFs_to_animations
 