

# Won't work:
remove.packages(“BioGeoBEARS”)
remove.packages(“BioGeoBears”)

# Yes will work:
remove.packages("BioGeoBEARS")

# Notice the difference between:
""
“”


# Try this:

remove.packages("rlang")
install.packages("rlang")
library(rlang)

library(devtools)
devtools::install_github(repo="nmatzke/BioGeoBEARS", INSTALL_opts="--byte-compile")