
source("functions/2_spillover_compensation.R")

options(nwarnings = 10000)

# run workflow step
step2_spillover_compensation()

message('\nDone!\n')
summary(warnings())
