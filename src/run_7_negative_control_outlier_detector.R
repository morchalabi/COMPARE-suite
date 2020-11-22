
source("functions/7_negative_control_outlier_detector")

options(nwarnings = 10000)      # shows all warnings (default is last 50)

# run workflow step
step7_negative_control_outlier_detector()

message('\nDone!\n')
summary(warnings())
