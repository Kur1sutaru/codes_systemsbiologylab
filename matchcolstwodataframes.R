# Match two columns of different dataframes and save in a new dataframe

setwd("C:/Users/GCVillalba/Downloads")

# Upper case all the info
`Drug_sensitivity_(PRISM_Repurposing_Primary_Screen)_19Q4` = toupper(`Drug_sensitivity_(PRISM_Repurposing_Primary_Screen)_19Q4`)


# To macth the collumns
deepmapdrugsensib <- merge(deepmapcellinfo, `Drug_sensitivity_(PRISM_Repurposing_Primary_Screen)_19Q4`,
             by.x = "depmap_id", by.y = "X" )

write.csv(deepmapdrugsensib, "deepmapdrugsensib.csv")



