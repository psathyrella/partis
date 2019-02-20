#code from Nick Matzke's code on BioGeoBEARS wiki
.events_txt_list_into_events_table<-function(events_txt_list, trtable=NULL, recalc_abs_ages=TRUE)
	{
	
	if (is.null(events_txt_list))
		{
		errortxt = paste("\nWARNING in events_txt_list_into_events_table(): your events_txt_list has NO events!\n\nThis means your tree has NO d/e/a events across the whole tree.\nThis is *expected* e.g. if you inferred d=e=0 under DEC+J. Input a list of '' or NA to avoid this error.\n\n", sep="")
		cat(errortxt)
		errortxt2 = paste("events_txt_list_into_events_table() is returning NULL which will might cause issues later.\n\n", sep="")
		cat(errortxt2)
		return(NULL)
		}
	
	# Convert NAs to "none"
	events_txt_list[is.na(events_txt_list)] = "none"
	
	# Remove lines with no events or NA:
	noneTF = events_txt_list == "none"
	keepTF = (noneTF == FALSE)
	events_txt_list = events_txt_list[keepTF]
	
	
	# If no anagenetic events, return NULL
	if (length(events_txt_list) == 0)
		{
		events_table = NULL
		return(events_table)
		}


	# Include the trtable, if that is input
	if (length(trtable) > 0)
		{
		trtable_subset = NULL
		}

	
	# Convert the events text back into a table:
	tmptable = NULL
	for (i in 1:length(events_txt_list))
		{
		#print(events_txt_list)
		tmptable_rows = .events_txt_into_events_table(events_txt_list[i])
		rownums_in_trtable = as.numeric(tmptable_rows$nodenum_at_top_of_branch)
		#print(tmptable_rows)
		num_newrows = nrow(tmptable_rows)
		tmptable = rbind(tmptable, tmptable_rows)
		} # END for (i in 1:length(events_txt_list))
	events_table = .dfnums_to_numeric(.adf2(tmptable))
	names(events_table) = c("nodenum_at_top_of_branch", "trynum", "brlen", "current_rangenum_1based", "new_rangenum_1based", "current_rangetxt", "new_rangetxt", "abs_event_time", "event_time", "event_type", "event_txt", "new_area_num_1based", "lost_area_num_1based", "dispersal_to", "extirpation_from")	
	
	return(events_table)
	}
	
### adf2 from CRAN version of BioGeoBEARS

.adf2<-function (x) 
{
    rownames = 1:nrow(x)
    return(as.data.frame(x, row.names = rownames, stringsAsFactors = FALSE))
}

### dfnums_to_numeric from CRAN version of BioGeoBEARS
.dfnums_to_numeric<-	function (dtf, max_NAs = 0.5, printout = FALSE, roundval = NULL) 
{
    dtf_classes = .cls.df(dtf, printout = FALSE)
    dtf_names = names(dtf)
    numcols = ncol(dtf)
    cls_col_list = c()
    for (i in 1:numcols) {
        cls_col = NA
        cmdstr = paste("cls_col = class(dtf$'", dtf_names[i], 
            "')", sep = "")
        eval(parse(text = cmdstr))
        cls_col_list[i] = cls_col
    }
    for (i in 1:numcols) {
        if (cls_col_list[i] == "list") {
            (next)()
        }
        if (cls_col_list[i] != "numeric") {
            newcol = NA
            cmdstr = paste("newcol = as.numeric(as.character(dtf$'", 
                dtf_names[i], "'))", sep = "")
            suppressWarnings(eval(parse(text = cmdstr)))
            if (sum(is.na(newcol)) < (max_NAs * length(newcol))) {
                cmdstr = paste("dtf$'", dtf_names[i], "' = newcol", 
                  sep = "")
                suppressWarnings(eval(parse(text = cmdstr)))
                if (!is.null(roundval)) {
                  cmdstr = paste("dtf$'", dtf_names[i], "' = round(dtf$'", 
                    dtf_names[i], "', digits=roundval)", sep = "")
                  suppressWarnings(eval(parse(text = cmdstr)))
                }
            }
        }
    }
    tmp_classes = .cls.df(dtf)
    dtf_classes$newclasses = tmp_classes[, ncol(tmp_classes)]
    if (printout) {
        cat("\n")
        cat("dfnums_to_numeric(dtf, max_NAs=", max_NAs, ") reports: dataframe 'dtf_classes' has ", 
            nrow(dtf_classes), " rows, ", ncol(dtf_classes), 
            " columns.\n", sep = "")
        cat("...names() and classes() of each column below...\n", 
            sep = "")
        cat("\n")
        print(dtf_classes)
    }
    return(dtf)
}


###cls.df from CRAN version of BioGeoBEARS
.cls.df<-function (dtf, printout = FALSE) 
{
    if (class(dtf) == "matrix") {
        dtf = as.data.frame(dtf, stringsAsFactors = FALSE)
    }
    dtf_names = names(dtf)
    numcols = ncol(dtf)
    cls_col_list = c()
    for (i in 1:numcols) {
        cls_col = NA
        cmdstr = paste("cls_col = class(dtf$'", dtf_names[i], 
            "')", sep = "")
        eval(parse(text = cmdstr))
        cls_col_list[i] = cls_col
    }
    colnum = 1:numcols
    dtf_classes = cbind(colnum, dtf_names, cls_col_list)
    dtf_classes = data.frame(dtf_classes, row.names = colnum)
    if (printout) {
        cat("\n")
        cat("cls.df(dtf) reports: dataframe 'dtf' has ", nrow(dtf), 
            " rows, ", numcols, " columns.\n", sep = "")
        cat("...names() and classes() of each column below...\n", 
            sep = "")
        cat("\n")
        print(dtf_classes)
        cat("\n")
    }
    return(dtf_classes)
}

#events_txt_into_events_table from Nick Matzke's code on BioGeoBEARS wiki
.events_txt_into_events_table<-function(branch_events_txt)
	{
	words = strsplit(branch_events_txt, split=";")[[1]]
	
	events_table_for_branch = t(sapply(X=words, FUN=.event_txt_into_events_row))
	row.names(events_table_for_branch) = NULL
	events_table_for_branch
	
	events_table_for_branch = .adf2(events_table_for_branch)
	events_table_for_branch
	names(events_table_for_branch) = c("nodenum_at_top_of_branch", "trynum", "brlen", "current_rangenum_1based", "new_rangenum_1based", "current_rangetxt", "new_rangetxt", "abs_event_time", "event_time", "event_type", "event_txt", "new_area_num_1based", "lost_area_num_1based", "dispersal_to", "extirpation_from")
	
	return(events_table_for_branch)
	}

#code from Nick Matzke's code on BioGeoBEARS wiki
.event_txt_into_events_row<-function(word)
	{

	split_key_item<-function(word2)
		{
		output_pair = c("", "")
		words3 = strsplit(word2, split=":")[[1]]
		numwords = length(words3)
		
		output_pair[1:numwords] = words3[1:numwords]
		
		return(output_pair)
		}


	words2 = strsplit(word, split=",")[[1]]
	output = sapply(X=words2, FUN=split_key_item)
	tmprow = matrix(data=output[2,], nrow=1)
	return(tmprow)
	}
