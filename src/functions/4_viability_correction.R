# This module removes intra- and inter-plate cell viability drift bias. Running it with no request for correction reveals possible cell viability bias.
# To remove cell viability bias, it needs to know the direction along which the bias has occurred like the order by which wells have been read.
# Each well is represented by the number of live cells.
# Input arguments are:
#   CORRECT (boolean): should the bias be corrected? like T/TRUE/true or F/FALSE/false
#   drctn_ (string): direction of bias like column or row
#   FITPLOT (boolean): should the regression plots be output? like T/TRUE/true or F/FALSE/false
#   HEATPLOT (boolean): should the plate heatmaps be output? like T/TRUE/true or F/FALSE/false
#   inURL (string): address to input data files like ../data
#   outURL (string): address to output result like ../out

# Algorithm designed and implemented by Mori C.H., mor.chalabi@gmail.com

require(flowCore)
require(pheatmap)
require(gridExtra)
require(ggplot2)

step4_viability_correction = function(CORRECT, drctn_, FITPLOT, HEATPLOT, inURL = '../data/', outURL = '../out/')
{
  # STEP 1: Computing cell viabilities ####
  
  # reading in annotation file
  
  annot_ = read.table(file = paste0(inURL,'/Annotations.txt'),header = T, sep = '\t',as.is = T, check.names = F, stringsAsFactors = F)
  plates_ = sort(unique(annot_$plate))
  
  # reading fcs files of each plate
  
  viability_mats = list()      # a list to store plate matrices storing number of live cells per well for all plates
  for(plate_ in plates_)
  {
    # reading fcs files
    message('Reading fcs files of plate #',plate_,':')
    
    # reading in annotation of files of current plate
    
    annot_tmp = annot_[annot_$plate %in% plate_,]     # annotations of wells (files) on current plate
    rows_ = sort(unique(annot_tmp$row))
    cols_ = sort(unique(annot_tmp$column))
    m_ = length(rows_)                                # number of rows in current plate
    n_ = length(cols_)                                # number of columns in current plate
    
    # reading in fcs file of each well
    
    viability_mat = matrix(data = 0, nrow = length(rows_), ncol = length(cols_))      # plate matrix to store number of live cells per well
    rownames(viability_mat) = rows_
    colnames(viability_mat) = cols_
    fcs_flNms = list()                                                                # a list to store well file names
    for(rw_ in 1:nrow(annot_tmp))
    {
      fcs_dt = read.FCS(filename = paste0(inURL,annot_tmp$file[rw_],'.fcs'), transformation = F, truncate_max_range = F)
  
      # computing offset from beginning of plate matrix
      
      i_ = which(rows_ == annot_tmp$row[rw_])                                   # row index
      j_ = which(cols_ == annot_tmp$column[rw_])                                # column index
      offset_ = if(drctn_ == 'column'){ m_*(j_-1)+i_ }else{ n_*(i_-1)+j_ }      # column-wise (left to right) or row-wise (left to right) offset
      viability_mat[i_,j_] = nrow(fcs_dt@exprs)                                 # computing viabilities
      fcs_flNms[[offset_]] = annot_tmp$file[rw_]                                # fcs file name in this offset of plate
    }
  
    # STEP 2: Intra-plate correcting fcs file ####
    
    if(CORRECT)
    {
      message('intra-plate correction of plate #',plate_,':')
  
      # step 2.1: computing intra-plate correction coefficients of current plate!
  
      y_ = if(drctn_ == 'column'){ matrix(data = viability_mat, ncol = 1, nrow = length(viability_mat), byrow = F)[,,drop = T] }else      # reading matrix column-wise
                                 { matrix(data = viability_mat, ncol = 1, nrow = length(viability_mat), byrow = T)[,,drop = T] }          # reading matrix row-wise
      x_ = 1:length(y_)                                                                                                                   # x_ is then offset of y_
  
      # removing outliers from plate viabilities before regression using inter-quantile range
  
      IQR_ = IQR(y_)
      quartiles_ = quantile(y_, probs = c(.25, .75))
      lowWhisker_ = max(min(y_), quartiles_[1] - IQR_*1.5)
      upWhisker_ = min(max(y_), quartiles_[2] + IQR_*1.5)
      y_tmp = y_[lowWhisker_ < y_ & y_ < upWhisker_]      # valid viabilities of current plate
  
      # fitting regression line
  
      fit_ = lm(data = data.frame(offset = 1:length(y_tmp), live_cells = y_tmp), formula = live_cells~offset)     # fitting simple linear regression
      a_ = fit_$coefficients["offset"]                                                                            # slope of the regression line
      if(a_ <= 0)
      {
        warning(paste0('Slope for plate #',plate_,' was not positive, no intra-plate correction is necessary!'))      # if slope was not positive, no correction
        viability_mats[[plate_]] = viability_mat
        next()                                                                                                        # goes to next plate
      }
      b_ = fit_$coefficients["(Intercept)"]                                                                       # intercept of the regression line
      alpha_ = b_/(a_*x_ + b_)                                                                                    # correction coefficients per well: y*/y = b/ax+b -> y* = y(b/ax+b), y* is translated y
  
      # step 2.2: correcting viabilities of current fcs
  
      for(offset_ in 1:length(fcs_flNms))      # fcs files in fcs_flNms are already ordered by their offset in viability matrix (viability_mat)
      {
        if(is.null(fcs_flNms[[offset_]])) { next() }                    # some wells could be empty on the plate! It goes to next well
        
        if(drctn_ == 'column')                                          # extracting (i,j) from column-wise offset
        {
          j_ = round(ceiling(offset_/m_))     # round is added because arithmetic on integer and double is not exact in R
          i_ = round(offset_ - m_*(j_-1))
        }else
        {
          i_ = round(ceiling(offset_/n_))
          j_ = round(offset_ - n_*(i_-1))
        }
        viability_mat[i_,j_] = viability_mat[i_,j_]*alpha_[offset_]     # y*/y = b/ax+b -> y* = y(b/ax+b)
      }
    }
    
    viability_mats[[plate_]] = viability_mat
  }
  
  if(FITPLOT)
  {
    # plotting regression plots for each plate per page
    
    if(CORRECT) { pdf(file = paste0(outURL,'/viability_fit_intra-plate_corrected.pdf')) }else{ pdf(file = paste0(outURL,'/viability_fit_intra-plate_no_correction.pdf')) }
    
    for(plate_ in plates_)
    {
      y_ = if(drctn_ == 'column'){ matrix(data = viability_mats[[plate_]], ncol = 1, nrow = length(viability_mats[[plate_]]), byrow = F)[,,drop = T]/1e4 }else      # reading matrix column-wise
                                 { matrix(data = viability_mats[[plate_]], ncol = 1, nrow = length(viability_mats[[plate_]]), byrow = T)[,,drop = T]/1e4 }          # reading matrix row-wise
      x_ = 1:length(y_)                                                                                                                                             # x_ is then offset of y_
      
      # removing viability outliers before regression analysis
      
      IQR_ = IQR(y_)
      quartiles_ = quantile(y_, probs = c(.25, .75))
      lowWhisker_ = max(min(y_), quartiles_[1] - IQR_*1.5)
      upWhisker_ = min(max(y_), quartiles_[2] + IQR_*1.5)
      y_tmp = y_[lowWhisker_ < y_ & y_ < upWhisker_]      # valid viabilities
      
      # fitting regression line
      
      fit_ = lm(data = data.frame(offset = 1:length(y_tmp), live_cells = y_tmp), formula = live_cells~offset)     # simple linear regression
      a_ = if(CORRECT){ 0 }else{ fit_$coefficients["offset"] }                                                    # slope of learned line
      b_ = fit_$coefficients["(Intercept)"]                                                                       # intercept of the fit line
      
      p_ = ggplot(data = data.frame(x = x_, y = y_, Plate = plate_), aes(x = x_, y = y_, color = as.factor(Plate)))+
              theme(axis.line = element_line(color = 'black'), panel.background = element_blank(), legend.key = element_rect(fill = NA),
                    text = element_text(face = 'bold',size = 20), plot.title = element_text(vjust = 0.5, hjust = 0.5), aspect.ratio = 1)+
              guides(color = guide_legend(override.aes = list(size = 5)))+
              labs(x = 'Pseudo-time', y= expression('Live cells ('%*%10^-4*')'), color = 'Plate')+
              geom_point(pch = 20, size = 2)+
              scale_color_manual(values = 'red')+
              geom_abline(intercept = b_, slope = a_, color = 'darkblue')+
              scale_x_continuous(expand = c(0.01, 0.01)) + scale_y_continuous(expand = c(0.01, 0.01))
      plot(p_)
    }
    graphics.off()
  }
  
  # STEP 3: Correcting for inter-plate bias ####
  
  if(CORRECT)
  {
    message('inter-plate correction')
    
    # computing inter-plate correction coefficient
    
    old_intcpts = double()                # a numeric vector to store baseline (median) of each plate
    for(plate_ in plates_) { old_intcpts[plate_] = median(viability_mats[[plate_]]) }
    new_intcpts = median(old_intcpts)     # new baseline across all plates
    
    # inter-plate correction
    
    for(plate_ in plates_)
    {
      # reading fcs files
      
      annot_tmp = annot_[annot_$plate %in% plate_,]     # annotations of well files of current plate
      rows_ = sort(unique(annot_tmp$row))               # number of rows in current plate
      cols_ = sort(unique(annot_tmp$column))            # number of rows in current plate
      for(rw_ in 1:nrow(annot_tmp))                     # for each well on current plate
      {
        fcs_dt = read.FCS(filename = paste0(inURL,annot_tmp$file[rw_],'.fcs'), transformation = F, truncate_max_range = F)     # event data of current well
        
        # computing offset from beginning of plate matrix
        
        i_ = which(rows_ == annot_tmp$row[rw_])         # row index of current well
        j_ = which(cols_ == annot_tmp$column[rw_])      # row index of current well
        
        # computing correction coefficient
        
        alpha_ = new_intcpts/old_intcpts[plate_]      # correction coefficient: b*/b = y*/y -> y* = y(b*/b)
  
        # correcting viabilities of current well
        
        viability_mats[[plate_]][i_,j_] = round(viability_mats[[plate_]][i_,j_]*alpha_)     # b*/b = y*/y -> y* = y(b*/b)
  
        # adding viability keyword to the current well's fcs file
        
        keyword(fcs_dt)[['viability']] = viability_mats[[plate_]][i_,j_]
        
        # writing corrected fcs files
        
        keyword(fcs_dt)[['$FIL']] = paste0(annot_tmp$file[rw_],'_compensated_corrected_corrected')      # updating $FIL keyword
        write.FCS(x = fcs_dt, filename = paste0(inURL,annot_tmp$file[rw_],'.fcs'))                      # matrix of events in this well
      }
      
      # writing viability matrices
      write.table(x = viability_mats[[plate_]], file = paste0(outURL,'/viability_mat_P',plate_,'.txt'), quote = F, sep = '\t', row.names = T, col.names = T)
    }
  }
  
  # plotting fit lines
  
  if(FITPLOT)
  {
    dt_ = list()                # a matrix with these columns: viability and plate
    for(plate_ in plates_)      # for each plate
    {
      tmp = if(drctn_ == 'column'){ matrix(data = viability_mats[[plate_]], ncol = 1, nrow = length(viability_mats[[plate_]]), byrow = F)[,,drop = T]/1e4 }else     # reading viability matrix column-wise
                                  { matrix(data = viability_mats[[plate_]], ncol = 1, nrow = length(viability_mats[[plate_]]), byrow = T)[,,drop = T]/1e4 }         # reading viability matrix row-wise
      dt_[[plate_]] = data.frame(viability = tmp, plate = as.character(plate_), stringsAsFactors = T)
    }
    dt_ = do.call(what = rbind, args = dt_)
    
    # plotting all dot plots on one PDF page
    
    if(CORRECT) { pdf(file = paste0(outURL,'/viability_fit_inter-plate_corrected.pdf')) }else{ pdf(file = paste0(outURL,'/viability_fit_inter-plate_no_correction.pdf')) }
    p_ = ggplot(data = dt_, aes(x = 1:nrow(dt_), y = viability, color = plate))+
            theme(axis.line = element_line(color = 'black'), panel.background = element_blank(), legend.key = element_rect(fill = NA),
                  text = element_text(face = 'bold',size = 20), plot.title = element_text(vjust = 0.5, hjust = 0.5), aspect.ratio = 1)+
            guides(color = guide_legend(override.aes = list(size = 5)))+
            labs(x = 'Pseudo-time', y= expression('Live cells ('%*%10^-4*')'), color = 'Plate')+
            geom_point(pch = 20, size = 2)+
            geom_hline(yintercept = median(dt_$viability), color = 'darkred')+
            scale_x_continuous(expand = c(0.01, 0.01)) + scale_y_continuous(expand = c(0.01, 0.01))
    plot(p_)
    graphics.off()
  }
  
  # STEP 4: Plotting plate heatmaps ####
  
  if(HEATPLOT)
  {
    numOfPlates = length(plates_)                             # number of plates in this assay
    nrow_ = floor(sqrt(numOfPlates))                          # number of rows per page
    ncol_ = ceiling(numOfPlates/floor(sqrt(numOfPlates)))     # number of columns per page
    width_ = height_ = max(nrow_, ncol_)*8.45                 # every heatmap takes up this width/height in inches
    
    if(!CORRECT)                                              # if no correction was requested
    {
      fl_ = paste0(outURL,'/plates_heatmap_no_correction_cyto.pdf')
      pdf(file = fl_, width = width_, height = height_)
    
      plot_list = list()
      
      # finding max in all plates
      
      max_ = 0
      for(plate_ in 1:numOfPlates)
      {
        max_tmp = max(viability_mats[[plate_]]/1e4)
        if(max_ < max_tmp) { max_ = max_tmp }
      }
      
      for(plate_ in 1:numOfPlates)      # plotting heatmaps of plates all in one PDF page
      {
        dt_mat = viability_mats[[plate_]]/1e4
        dt_mat = rbind(dt_mat, MAX = c(max_, rep(0,ncol(dt_mat)-1)))      # last row of each plate heatmap contain max value across all plates
        
        # normalizing heatmap's color palette
        
        dt_ = unique(as.numeric(dt_mat))
        IQR_ = IQR(dt_)
        quartiles_ = quantile(dt_, probs = c(.25, .75))
        lowWhisker_ = max(min(dt_), quartiles_[1] - IQR_*1.5)
        upWhisker_ = min(max(dt_), quartiles_[2] + IQR_*1.5)
        
        dt_ = sort(dt_, decreasing = F)
        
        cols1_ = u1_ = NULL
        inds_ = which( dt_ <= lowWhisker_)
        len_ = length(inds_)
        if(0 < len_)
        {
          cols1_ = colorRampPalette(colors = c('grey80','red','red4'))(len_)
          u1_ = dt_[inds_[length(inds_)]]
          if(u1_ == 0) { u1_ = NULL }     # to avoid redundant legend key on heatmap
        }
        
        cols2_ = u2_ = NULL
        inds_ = which(lowWhisker_ < dt_ & dt_ < quartiles_[1])
        len_ = length(inds_)
        if(0 < length(inds_))
        {
          cols2_ = colorRampPalette(colors = c('orange','orange4'))(len_)
          u2_ = dt_[inds_[length(inds_)]]
        }
        
        cols3_ = u3_ = NULL
        inds_ = which(quartiles_[1] <= dt_ & dt_ <= quartiles_[2])
        len_ = length(inds_)
        if(0 < length(inds_))
        {
          cols3_ = colorRampPalette(colors = c("yellow",'yellow4'))(len_)
          u3_ = dt_[inds_[length(inds_)]]
        }
        
        cols4_ = u4_ = NULL
        inds_ = which(quartiles_[2] < dt_ & dt_ <= upWhisker_)
        len_ = length(inds_)
        if(0 < length(inds_))
        {
          cols4_ = colorRampPalette(colors = c('green','green4'))(len_)
          u4_ = dt_[inds_[length(inds_)]]
        }
        
        cols5_ = NULL
        u5_ = if(max(dt_) < round(max(dt_),2)) { round(max(dt_)-0.01,2) }else{ round(max(dt_),2) }
        inds_ = which(upWhisker_ < dt_)
        len_ = length(inds_)
        if(0 < length(inds_))
        {
          cols5_ = colorRampPalette(colors = c('skyblue','blue4'))(len_)
        }
        
        cols_ = c(cols1_, cols2_, cols3_, cols4_, cols5_)
        col_step = 2*(diff(range(dt_))/length(cols_))
        
        p_ = pheatmap(mat = dt_mat,
                      cluster_cols = F,
                      cluster_rows = F,
                      show_rownames = T, cellheight = 20,
                      show_colnames = T, cellwidth = 20,
                      main = paste0('Plate: ',plate_),
                      fontsize = 7,
                      border_color = 'grey90', color = cols_,
                      breaks = c(min(dt_)-col_step,dt_),     # in pheatmap color intervals (showing lower and uper bounds) are open on the left and closed on the right
                      legend = T,
                      legend_breaks = round(c(0,u1_,u2_, u3_, u4_, u5_),2),
                      annotation_col = NA, silent = T)
        plot_list[[plate_]] = p_[[4]]
      }
      grid.arrange(arrangeGrob(grobs= plot_list, nrow = nrow_, ncol = ncol_))
    
      graphics.off()
      
    }else     # if correction was requested
    {
      fl_ = paste0(outURL,'/plates_heatmap_correction_cyto.pdf')
      pdf(file = fl_, width = width_, height = height_)
    
      plot_list = list()
  
      # finding max in all plates
      
      max_ = 0
      dt_ = NULL
      for(plate_ in 1:numOfPlates)
      {
        dt_ = c(dt_, as.numeric(viability_mats[[plate_]]/1e4))
      }
      max_ = max(dt_)
      
      # normalizing heatmap's color palette
      
      dt_ = unique(dt_)
      IQR_ = IQR(dt_)
      quartiles_ = quantile(dt_, probs = c(.25, .75))
      lowWhisker_ = max(min(dt_), quartiles_[1] - IQR_*1.5)
      upWhisker_ = min(max(dt_), quartiles_[2] + IQR_*1.5)
      
      dt_ = sort(dt_, decreasing = F)
      
      cols1_ = u1_ = NULL
      inds_ = which( dt_ <= lowWhisker_)
      len_ = length(inds_)
      if(0 < len_)
      {
        cols1_ = colorRampPalette(colors = c('grey80','red','red4'))(len_)
        u1_ = dt_[inds_[length(inds_)]]
        if(u1_ == 0) { u1_ = NULL }     # to avoid redundant legend key on heatmap
      }
      
      cols2_ = u2_ = NULL
      inds_ = which(lowWhisker_ < dt_ & dt_ < quartiles_[1])
      len_ = length(inds_)
      if(0 < length(inds_))
      {
        cols2_ = colorRampPalette(colors = c('orange','orange4'))(len_)
        u2_ = dt_[inds_[length(inds_)]]
      }
      
      cols3_ = u3_ = NULL
      inds_ = which(quartiles_[1] <= dt_ & dt_ <= quartiles_[2])
      len_ = length(inds_)
      if(0 < length(inds_))
      {
        cols3_ = colorRampPalette(colors = c("yellow",'yellow4'))(len_)
        u3_ = dt_[inds_[length(inds_)]]
      }
      
      cols4_ = u4_ = NULL
      inds_ = which(quartiles_[2] < dt_ & dt_ <= upWhisker_)
      len_ = length(inds_)
      if(0 < length(inds_))
      {
        cols4_ = colorRampPalette(colors = c('green','green4'))(len_)
        u4_ = dt_[inds_[length(inds_)]]
      }
      
      cols5_ = NULL
      u5_ = if(max(dt_) < round(max(dt_),2)) { round(max(dt_)-0.01,2) }else{ round(max(dt_),2) }
      inds_ = which(upWhisker_ < dt_)
      len_ = length(inds_)
      if(0 < length(inds_))
      {
        cols5_ = colorRampPalette(colors = c('skyblue','blue4'))(len_)
      }
      
      cols_ = c(cols1_, cols2_, cols3_, cols4_, cols5_)
      col_step = 2*(diff(range(dt_))/length(cols_))
  
      for(plate_ in 1:numOfPlates)      # plotting heatmaps of plates all in one PDF page
      {
        LEGEND_SHOW = F
        if(plate_ == 1) { LEGEND_SHOW = T }
        
        dt_mat = viability_mats[[plate_]]/1e4
        dt_mat = rbind(dt_mat, MAX = c(max_, rep(0,ncol(dt_mat)-1)))
        
        p_ = pheatmap(mat = dt_mat,
                      cluster_cols = F,
                      cluster_rows = F,
                      show_rownames = T, cellheight = 20,
                      show_colnames = T, cellwidth = 20,
                      main = paste0('Plate: ',plate_),
                      fontsize = 7,
                      border_color = 'grey90', color = cols_,
                      breaks = c(min(dt_)-col_step,dt_),     # in pheatmap color intervals (showing lower and uper bounds) are open on the left and closed on the right
                      legend = LEGEND_SHOW,
                      legend_breaks = round(c(0,u1_,u2_, u3_, u4_, u5_),2),
                      annotation_col = NA, silent = T)
        plot_list[[plate_]] = p_[[4]]
      }
      grid.arrange(arrangeGrob(grobs= plot_list, nrow = nrow_, ncol = ncol_))
    
      graphics.off()
    }
  }

  return(NULL)
}
