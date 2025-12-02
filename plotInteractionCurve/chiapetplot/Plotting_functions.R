## =================================================================== ##
## Set layout for plotting
## =================================================================== ##
plot_layout <- function(nrow, ncol, heights = NULL, widths = NULL) {
    # Function to test if sth is a interger
    is.wholenumber <-
        function(x, tol = .Machine$double.eps ^ 0.5)  abs(x - round(x)) < tol
    
    # Check if nrow or ncol are "null unit"
    if (is.null(heights)) {
        heights = unit(rep_len(1, nrow), "null")
    } 
    if (is.null(widths)) {
        widths = unit(rep_len(1, ncol), "null")
    }
    
    # Set layput for later plotting
    grid.newpage()
    layout <- grid.layout(nrow = nrow, ncol = ncol, 
                          heights = heights,
                          widths = widths)
    vp_overall <- viewport(name = "vp_overall", width = unit(0.95, "npc"), height = unit(0.95, "npc"), layout = layout)
    pushViewport(vp_overall)
}

## =================================================================== ##
## function to draw intra chromosomal interactions
## =================================================================== ##
plot_intra_chr_interaction <- function(interaction, cyto, peaks, return = FALSE) {
    # Remove interactions on chrM
    interaction <- interaction[as.character(interaction[[1]]) != "chrM" & as.character(interaction[[4]]) != "chrM", ]
    peaks <- peaks[peaks[[1]] != "chrM", ]
    cyto <- cyto[cyto[[1]] != "chrM", ]
    
    # Set color scheme for each "stain" in cytoband data
    colnames(cyto) <- c("chrom", "start", "end", "name", "stain")
    cyto_color <- data.frame(stain = c("gpos100", "gpos75", "gpos66", "gpos50", "gpos33", "gpos25", "gneg", "acen", "gvar", "stalk"),
                             col = c("#4D4D4D", "#969696", "#AAAAAA", "#C3C3C3", "#C2C2C2", "#E6E6E6", "white", "white", "#F3E6C4", "#5A81AA"))
    cyto <- merge(cyto, cyto_color, by = "stain", all.x = T)
    
    # Order chromosome from chr1 to chrY
    idx <- grepl("[0-9]+", cyto$chrom)
    chrom_no <- paste("chr", unique(c(sort(as.numeric(gsub("chr", "", cyto$chrom[idx]))), gsub("chr", "", sort(cyto$chrom[!idx])))), sep = "")
    
    # Whether return color columns for RCiscos plot
    # interaction <- interaction[, 1:7]  
    interaction$PlotColor <- apply(interaction, 1, function(x){x[7] <- as.numeric(x[7]); if(x[7] == 2){a = "#EF667A"}
                                                               else if(x[7] == 7){a = "#FF0000"}   #E92F4A
                                                               else if(x[7] == 8 | x[7] == 9){a = "#AE2F4A"}
                                                               else if(x[7] >=10 & x[7] <= 15){a = "#7C2F4A"}   
                                                               else if(x[7] >= 16 & x[7] <= 18){a = "#7C2F4A"}
                                                               else(a = "#5A2F63")})
    if (return == TRUE) {
        return(interaction)
    } else {
      #  interaction <- interaction[, c(1:6, 8)]
        # Order cytoband data
        cyto$chrom_no_f <- factor(cyto$chrom, levels = chrom_no)
        cyto <- cyto[with(cyto, order(chrom_no_f, start)), ]
        # Seperate centromere data from cytoband data
        centro <- cyto[cyto$stain == "acen", ]
        centro <- centro[, c("chrom", "start", "end")]
        cyto <- cyto[cyto$stain != "acen", ]
        cyto <- cyto[, c("chrom", "start", "end", "col")]
        
        # Convert variables to correct type
        cyto[, c(1, 4)] <- sapply(cyto[, c(1, 4)], as.character)
        cyto[, 2:3] <- sapply(cyto[, 2:3], as.numeric)
        centro[, 1] <- as.character(centro[, 1])
        centro[, 2:3] <- sapply(centro[, 2:3], as.numeric)
        peaks[, 1] <- as.character(peaks[, 1])
        peaks[, 2:4] <- sapply(peaks[, 2:4], as.numeric)
        
        # Calculate max coordinate for each chromosome and the max peak signal
        cyto_max <- max(cyto[["end"]])
        max_score <- max(peaks[[4]])
        
        # Draw ideogram, interaction curves and peaks signals for each chromosome
        for (i in seq.int(length(chrom_no))) { 
            # Get chromosome name
            if (grepl("X", chrom_no[i])) {
                chrom <- paste("chr", "X", sep = "")
            } else if (grepl("Y", chrom_no[i])) {
                chrom <- paste("chr", "Y", sep = "") 
            }  else { chrom <- paste("chr", i, sep = "") }
            
            # Filter data for ideogram
            chrom_cyto <- cyto[which(cyto[[1]] == chrom), ]
            chrom_cyto <- chrom_cyto[order(chrom_cyto[[2]]), ]
            chrom_max <- max(chrom_cyto[["end"]])
            chrom_centro <- centro[which(centro[["chrom"]] == chrom), ]
            
            # Filter data for peaks
            chrom_peaks <- peaks[which(peaks[[1]] == chrom), ]
            chrom_peaks <- chrom_peaks[order(chrom_peaks[[2]]), ]
            chrom_peaks[[4]] <- chrom_peaks[[4]] / max_score * 2
            chrom_peaks[[4]][chrom_peaks[[4]] > 1] <- 1
            
            # Filter data for inetraction curves
            chrom_interaction <- interaction[which(interaction[[1]] == chrom), ]
            chrom_interaction <- chrom_interaction[, c(2, 3, 5, 6, 8)]
            
            # Plot chromosome name
            vp_chr = viewport(layout.pos.row = 3 * i - 1, layout.pos.col = 1)
            grid.text(chrom, x = unit(0.5, "npc"), y = unit(0, "npc"),
                      just = c("right", "bottom"), gp = gpar(col = "#323745", fontsize = 15) , vp = vp_chr) 
            
            # Plot chromosome ideogram
            vp = viewport(layout.pos.row = 3 * i - 1, layout.pos.col = 2, name = "vp", xscale = c(0, cyto_max))
            grid.rect(x = unit(chrom_cyto[["start"]], "native"), y = unit(1, "npc"), width = unit(chrom_cyto[["end"]] - chrom_cyto[["start"]], "native"), height = unit(1, "npc"),
                                just = c("left", "top"), gp = gpar(fill = chrom_cyto[["col"]], col = NA, col = "#9F9F9F"), vp = vp)    

            if (nrow(chrom_centro) == 0) {
                grid.rect(x = unit(0, "npc"), y = unit(1, "npc"), width = unit(chrom_max, "native"), height = unit(1, "npc"),
                                            just = c("left", "top"), gp = gpar(fill = NA, col = "#9F9F9F"), vp = vp)   

            } else {
                # Left part
                chrom_left <- chrom_centro[1, ]
                chrom_left$start <- 0
                grid.rect(x = unit(0, "npc"), y = unit(1, "npc"), width = unit(chrom_left[["end"]] - chrom_left[["start"]], "native"), height = unit(1, "npc"),
                                                just = c("left", "top"), gp = gpar(fill = NA, col = "#9F9F9F"), vp = vp)   

                # Right part
                chrom_right <- chrom_centro[2, ]
                chrom_right$end <- chrom_max
                grid.rect(x = unit(chrom_right[["start"]], "native"), y = unit(1, "npc"), width = unit(chrom_right[["end"]] - chrom_right[["start"]], "native"), height = unit(1, "npc"),
                                                just = c("left", "top"), gp = gpar(fill = NA, col = "#9F9F9F"), vp = vp)
            }

            draw_Bezier_Curve <- function(vector) {
              x_start = (as.numeric(vector[1]) + as.numeric(vector[2])) / 2
              x_end = (as.numeric(vector[3]) + as.numeric(vector[4])) / 2
              x1 <- x_start; y1 <- -0.2;
              x4 <- x_end; y4 <- -0.2;
              if (abs(x4 - x1) <= 500000) {
                x2 <- (x_start + x_end)/2 + (x_start - x_end)/2 * 30
                y2 <- -3/3
                x3 <- (x_start + x_end)/2 - (x_start - x_end)/2 * 30
                y3 <- -3/3
              } else if(abs(x4 - x1) <= 1000000) {
                x2 <- (x_start + x_end)/2 + (x_start - x_end)/2 * 6
                y2 <- -5/3
                x3 <- (x_start + x_end)/2 - (x_start - x_end)/2 * 6
                y3 <- -5/3
              } else if(abs(x4 - x1) <= 10000000) {
                x2 <- (x_start + x_end)/2 + (x_start - x_end)/2 * 0.2
                y2 <- -7/3
                x3 <- (x_start + x_end)/2 - (x_start - x_end)/2 * 0.2
                y3 <- -7/3
              } else {
                x2 <- (x_start + x_end)/2 + (x_start - x_end)/2 * 0.2
                y2 <- -9/3
                x3 <- (x_start + x_end)/2 - (x_start - x_end)/2 * 0.2
                y3 <- -9/3
              }
              # print(vector)
              bezierGrob(x = unit(c(x1, x2, x3, x4), "native"), y = c(y1, y2, y3, y4),
                gp = gpar(col = as.character(vector[5]), fill = as.character(vector[5]), alpha = 0.5, lwd = 0.5))
            }
            
            curve <- apply(chrom_interaction, 1, draw_Bezier_Curve)
            curves <- gList()
            for (k in seq.int(curve)) {
                curves[[k]] <- gList(curve[[k]])
            }
            vp_interactions <- viewport(layout.pos.row = 3 * i - 1, layout.pos.col = 2, name = "vp_interactions", xscale = c(0, cyto_max), just = c("left", "top"))
            pushViewport(vp_interactions)
            grid.draw(curves)
            upViewport()
            
            #  Plot peaks
            vp_peaks <- viewport(layout.pos.row = 3 * i - 2, layout.pos.col = 2, name = "vp_peaks", xscale = c(0, cyto_max), yscale = c(0, max_score), just = c("left", "bottom"))
            if (nrow(chrom_peaks) == 0) {
                next
            } else {
                grid.rect(x = unit(chrom_peaks[[2]], "native"), y = unit(0.09, "npc"), width = unit(chrom_peaks[[3]] - chrom_peaks[[2]], "native"), height = unit(chrom_peaks[[4]]/2, "npc"),
                          just = c("left", "bottom"), gp = gpar(fill = "#77C4D2", col = "#77C4D2"), vp = vp_peaks)
            }


        }
    }
}