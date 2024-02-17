##### CURATE OTUs WITH LULU #####
# These OTUs were made from the ASVs above, this is a second round of lulu

otu_table = read.delim("FINAL_FOR_LULU_otutab.txt", quote = "", row.names = 1)
names(otu_table) <- gsub("\\.(?=.*_)","-",names(otu_table), perl = TRUE) # At some point in the LULU pipeline, all of the hyphens got replaced with periods, so we're just going to go back in and replace the periods with hyphens again to maintain compatibility.
otu_table = as.data.frame(otu_table)

matchlist = read.delim("match_list.txt", stringsAsFactors = F)

curated_result = lulu(otu_table, matchlist, minimum_ratio_type = "min", minimum_ratio = 1, minimum_match = 84, minimum_relative_cooccurence = 0.95)