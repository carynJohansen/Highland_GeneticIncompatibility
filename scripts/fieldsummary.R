# Understanding the field report

# https://docs.google.com/spreadsheets/d/1OPwonbK8CAdtxbd8sERE02LdGr6NAtsETawSIJ4UvOg/edit#gid=0

# The Runcie panel

male <- c("SA_16-Low","SA_16-High","Check-Low","Check-High","SA_1-Low","SA_1-High","SA_8-High",
          "SA_8-Low","SA_21-Low","SA_21-High","SA_32-High","SA_32-Low","Mex_8-High","Mex_8-Low",
          "Mex_29-High","Mex_29-Low","Mex_18-High","Mex_18-Low","Mex_3-Low","Mex_3-High",
          "Check-Low","Check-High","Mex_17-Low","Mex_17-High","Mex_10-High","Mex_10-Low",
          "Mex_19-High","Mex_19-Low","Mex_22-Low","Mex_22-High","SA_26-Low","SA_26-High",
          "SA_22-High","SA_22-Low","SA_25-Low","SA_25-High","Mex_16-Low","Mex_16-High",
          "SA_33-High","SA_33-Low","Mex_7-Low","Mex_7-High","SA_27-High","SA_27-Low",
          "Mex_24-High","Mex_24-Low","Mex_15-Low","Mex_15-High","SA_24-Low","SA_24-High",
          "SA_29-Low","SA_29-High","Check-High","Check-Low","SA_2-High","SA_2-Low","Mex_20-High",
          "Mex_20-Low","Mex_25-High","Mex_25-Low","Check-High","Check-Low","Mex_9-Low",
          "Mex_9-High","SA_14-High","SA_14-Low","SA_15-Low","SA_15-High","Mex_14-Low",
          "Mex_14-High","Mex_5-Low","Mex_5-High","Mex_4-High","Mex_4-Low","SA_4-Low",
          "SA_4-High","Mex_1-High","Mex_1-Low","SA_6-Low","SA_6-High","Mex_23-Low","Mex_23-High",
          "Mex_28-High","Mex_28-Low","Check-Low","Check-High","Mex_11-High","Mex_11-Low",
          "SA_7-High","SA_7-Low","Mex_6-High","Mex_6-Low","SA_20-High","SA_20-Low","SA_23-Low",
          "SA_23-High","SA_19-Low","SA_19-High","SA_28-Low","SA_28-High","SA_30-High",
          "SA_30-Low","SA_3-Low","SA_3-High","Mex_12-Low","Mex_12-High","SA_9-High","SA_9-Low",
          "Mex_30-Low","Mex_30-High","Check-Low","Check-High","SA_17-Low","SA_17-High",
          "Mex_13-Low","Mex_13-High","Mex_31-High","Mex_31-Low","Mex_26-High","Mex_26-Low",
          "SA_7-Low","SA_7-High","Mex_6-High","Mex_6-Low","SA_1-Low","SA_1-High","Mex_24-High",
          "Mex_24-Low","Check-Low","Check-High","Mex_12-Low","Mex_12-High","SA_3-Low",
          "SA_3-High","Mex_7-Low","Mex_7-High","Mex_8-High","Mex_8-Low","SA_23-Low",
          "SA_23-High","SA_16-Low","SA_16-High","SA_28-Low","SA_28-High","Mex_18-High",
          "Mex_18-Low","Check-High","Check-Low","Mex_19-High","Mex_19-Low","Mex_16-Low",
          "Mex_16-High","SA_2-High","SA_2-Low","SA_4-High","SA_4-Low","Mex_13-High",
          "Mex_13-Low","Mex_11-Low","Mex_11-High","SA_32-Low","SA_32-High","Mex_26-Low",
          "Mex_26-High","Mex_9-Low","Mex_9-High","Mex_14-High","Mex_14-Low","Mex_1-Low",
          "Mex_1-High","SA_24-Low","SA_24-High","SA_21-High","SA_21-Low","Check-Low",
          "Check-High","SA_17-High","SA_17-Low","Mex_3-Low","Mex_3-High","Mex_30-High",
          "Mex_30-Low","Mex_29-High","Mex_29-Low","SA_6-High","SA_6-Low","SA_8-Low",
          "SA_8-High","Check-High","Check-Low","SA_14-Low","SA_14-High","Mex_28-Low",
          "Mex_28-High","Mex_23-Low","Mex_23-High","SA_19-High","SA_19-Low","SA_15-High",
          "SA_15-Low","Mex_10-Low","Mex_10-High","SA_30-High","SA_30-Low","SA_20-Low",
          "SA_20-High","SA_27-High","SA_27-Low","Check-High","Check-Low","Mex_4-High",
          "Mex_4-Low","Mex_25-Low","Mex_25-High","SA_9-High","SA_9-Low","Mex_15-Low",
          "Mex_15-High","SA_29-Low","SA_29-High","Mex_22-Low","Mex_22-High","SA_26-Low",
          "SA_26-High","Mex_31-Low","Mex_31-High","SA_25-High","SA_25-Low","Mex_5-Low",
          "Mex_5-High","Mex_20-Low","Mex_20-High","SA_33-High","SA_33-Low","Check-High",
          "Check-Low","Mex_17-High","Mex_17-Low","SA_22-High","SA_22-Low")
table(male)

male2 <- gsub("-Low", "", male)
male2 <- gsub("-High", "", male2)

table(male2)
#4 plants from each Male parent, 2 each for high and low.

# What does the high and low mean here? Is this the latitudinal pairing?

#How many Mexican and how many South American?

male_mex <- male2[grep("Mex", male2)]
male_SA <- male2[grep("SA", male2)]

length(table(male_mex)) #28
length(table(male_SA)) #26


