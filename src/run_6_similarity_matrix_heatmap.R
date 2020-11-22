
source("functions/6_similarity_matrix_heatmap.R")

options(nwarnings = 10000)      # shows all warnings (default is last 50)

# run workflow step
step6_similarity_matrix_heatmap()

message('\nDone!\n')
summary(warnings())
