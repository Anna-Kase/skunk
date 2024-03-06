

# Striped skunk dataset prep 

# can be sourced, but shouldn't need to be after writing new file

# read in old and new datasets
og <- read.csv("./data/full_capture_history.csv")

updated <- read.csv("../../../Downloads/full_capture_history.csv")

# remove two columns that the new dataset had that the old one did not
updated <- updated %>% dplyr::select(-c(Start, End))

# full join both datasets
both <- dplyr::left_join(updated, og)

# remove sites that were too close to other sites
both <- both %>% dplyr::filter(!Site %in% c("C05-EGG1", "C05-BMP2", "D02-MOP1",
                                            "D05-BRP1", "D07-NOC1", "R06-SLT2",
                                            "R06-SLT3", "S07-SAG2", "S07-STJ1",
                                            "S07-WGL0"))

# copy the J column so one can be the original and the new can
# have zeros converted to NAs
both <- both %>%
  mutate(naJ = J)

# convert zeros in new J column to be NAs
both["naJ"][both["naJ"] == 0] <- NA

# calculate the mean number of NAs in the new J column -
# these will later become NAs in the detection data (Y) as
# well - but they will be estimated by the model
perc_na <- both %>%
  group_by(Season) %>%
  summarize(
    perc_na <- mean(is.na(naJ))
  )

# remove a few early seasons of data to create a complete seasonal
# detection history
# JU16 also was removed because it was not behaving

complete <- both %>% 
  dplyr::filter(!Season %in% c("JA11", "JA13", "JU13", "JA14", "JU16"))

# save this dataset as a .csv because it needs to be used in the 
# scripts that create the covariate document so the sites are correct
# and everything matches
write.csv(complete, "./data/complete_data.csv")


