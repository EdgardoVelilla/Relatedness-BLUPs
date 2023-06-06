justCheck <- function(pedigree) {
     library(data.table, quietly = TRUE)

     if (is.data.frame(pedigree)) {
          pedigree <- as.data.table(pedigree)
     } else if (is.matrix(pedigree)) {
          pedigree <- as.data.table(pedigree)
     } else if (!is.data.table(pedigree)) {
          stop("nothing to do...")
     }

     ped3 <- pedigree[, c(1:3)]
     setnames(ped3, c("TreeID", "mum", "dad"))

     # Replace 0s in mother and father columns with NAs
     ped3[mum == 0, mum := NA]
     ped3[dad == 0, dad := NA]

     # Check that there are no missing mothers and fathers for all individuals
     if (all(is.na(ped3[, mum])) & all(is.na(ped3[, dad]))) {
          stop("All mothers and fathers are missing")
     }

     # Check that all mothers and fathers are present in the ID column
     if (sum((na.omit(ped3[, mum]) %in% ped3[, TreeID]) == FALSE) > 0 &
          any(is.na(ped3[, mum]) == FALSE)) {
          warning("individuals appearing as mother but not in pedigree")
     }

     if (sum((na.omit(ped3[, dad]) %in% ped3[, TreeID]) == FALSE) > 0 &
          any(is.na(ped3[, dad]) == FALSE)) {
          warning("individuals appearing as father but not in pedigree")
     }

     # Check for duplicated individuals in the TreeID column
     if (sum(duplicated(ped3[, TreeID])) > 0) {
          stop("some individuals appear more than once in the pedigree")
     }

     # Remove rows with missing or zero IDs
     if (any(ped3[, TreeID] == 0 | is.na(ped3[, TreeID]))) {
          warning("Missing value in the TreeID column - row will be discarded")

          # ped3 <- ped3[-which(ped3[, TreeID] == 0 | is.na(ped3[, TreeID])), ]
     }

     # # Identify individuals with missing mother
     # unique.mother <- unique(ped3[!is.na(mum), mum])
     # miss.mother <- unique.mother[which(is.na(match(unique.mother, ped3[, TreeID])))]

     # # Identify individuals with missing father
     # unique.father <- unique(ped3[!is.na(dad), dad])
     # miss.father <- unique.father[which(is.na(match(unique.father, ped3[, TreeID])))]

     # # Add missing mothers and fathers as new rows to pedigree
     # if(length(miss.mother) == 0 & length(miss.father) == 0) {
     #    ped.adjusted <- ped3
     # } else {
     #    addPed <- data.table(c(miss.mother, miss.father),
     #                         rep(NA, length(miss.mother) + length(miss.father)),
     #                         rep(NA, length(miss.mother) + length(miss.father)),
     #                         matrix(NA, nrow = (length(miss.mother) + length(miss.father)),
     #                         ncol = ncol(ped3) - 3))
     # 	setnames(addPed, c('TreeID','mum', 'dad'))
     #    ped.adjusted <- rbind(addPed, ped3)
     #     cat('pedigree was fixed...! \n')
     # }

     # # Replace NAs with 0
     # ped.adjusted[, mum:= ifelse(is.na(mum),0,mum)]
     # ped.adjusted[, dad:= ifelse(is.na(dad),0,dad)]

     # # Orders the pedigree so that offspring follow parents
     # ord <- pedigree::orderPed(ped.adjusted)
     # ped.adjusted <- ped.adjusted[order(ord),]

     #  ped.adjusted[]
}
