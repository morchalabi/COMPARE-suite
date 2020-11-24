# This module removes intra- and inter-plate signal drift bias. Running it with no request for correction reveals possible sources of bias.
# Possible sources of bias are like edge effect, signal drift, cell viability, auto-fluorescence and carry-over effect.
# To remove signal drift bias, it needs to know the direction along which the bias has occurred like the order by which wells have been read.
# Each well is represented by the median of each marker.
# Input arguments are:
#   chnls_ (quoted string): channel (not marker) names like 'chnl1,chn2,chnl3'
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

step3_signal_drift_correction = function(chnls_, CORRECT, drctn_, FITPLOT, HEATPLOT, inURL = '../data/', outURL = '../out/')
{
  # STEP 1: Computing MFIs ####
  
  # reading in annotation file
  
  annot_ = read.table(file = paste0(inURL,'/Annotations.txt'), header = T, sep = '\t',as.is = T, check.names = F, stringsAsFactors = F)
  plates_ = sort(unique(annot_$plate))      # how many plates are there in this assay?
  
  # reading fcs files of each plate
  
  MFI_mats = list()     # a list containing MFI matrices of each plate per channel
  for(plate_ in plates_)
  {
    message('Reading in fcs files of plate #',plate_,':')
  
    # reading in annotation of files of current plate
    
    annot_tmp = annot_[annot_$plate %in% plate_,]     # part of annotation table containing info for current plate
    rows_ = sort(unique(annot_tmp$row))               # rows in current plate
    cols_ = sort(unique(annot_tmp$column))
    m_ = length(rows_)                                # number of rows in current plate
    n_ = length(cols_)                                # number of columns in current plate
    
    # creating empty MFI matrices for current plate for each channel
    
    MFI_mat = list()      # a list containing MFI matrix of current plate per channel
    for(chnl_ in chnls_)
    {
      mat_ = matrix(data = 0, nrow = length(rows_), ncol = length(cols_))
      rownames(mat_) = rows_
      colnames(mat_) = cols_
      MFI_mat[[chnl_]] = mat_
    }
    
    # reading in fcs files (wells) of current plate and computing their MFIs, one per channel
    
    fcs_dt_ls = list()          # a list containing all fcs files (wells) of current plate
    fcs_flNms = character()     # all fcs file names
    for(rw_ in 1:nrow(annot_tmp))
    {
      fcs_dt = read.FCS(filename = paste0(inURL,annot_tmp$file[rw_],'.fcs'), transformation = F, truncate_max_range = F)
  
      # computing offset from the beginning of the MFI matrix
      
      i_ = which(rows_ == annot_tmp$row[rw_])                                   # row index
      j_ = which(cols_ == annot_tmp$column[rw_])                                # column index
      offset_ = if(drctn_ == 'column'){ m_*(j_-1)+i_ }else{ n_*(i_-1)+j_ }      # column-wise (left to right) or row-wise (left to right) offset
      fcs_dt_ls[[offset_]] = fcs_dt                                             # fcs data at location (i,j), offset, of current plate
      fcs_flNms[offset_] = annot_tmp$file[rw_]                                  # fcs file name at location (i,j), offset, of current plate
  
      # computing one MFI per channel for current fcs file at position (i,j)
      # outlier event values (like cell marker expressions) are removed first; Because here we are not assigning negatives to 0, it is fine to use IQR for removing outliers, if any
      
      for(chnl_ in chnls_)
      {
        dt_tmp = fcs_dt_ls[[offset_]]@exprs[,chnl_, drop = T]                                     # event values of the fcs file in current offset
        dt_tmp = dt_tmp[which(0 <= dt_tmp)]                                                       # N.B.: which() ignores NA and NaNs. Non-positives are always non-positives even in the presence of drift
        IQR_ = IQR(dt_tmp)                                                                        # inter-quantile range
        quartiles_ = quantile(dt_tmp, probs = c(.25, .75))                                        # 25th and 75th percentiles
        lowWhisker_ = max(min(dt_tmp), quartiles_[1] - IQR_*1.5)                                  # lower whisker
        upWhisker_ = min(max(dt_tmp), quartiles_[2] + IQR_*1.5)
        MFI_mat[[chnl_]][i_,j_] = median(dt_tmp[lowWhisker_ < dt_tmp & dt_tmp < upWhisker_])      # median is taken over valid event values which are within lower-whisker and upper-whiskers
      }
    }
  
    # STEP 2: Intra-plate correcting of fcs files (wells) in current plate per channel ####
    
    if(CORRECT)
    {
      message('intra-plate correction of plate #',plate_,':')
  
      # computing intra-plate correction coefficients (alpha) for each channel of current plate
  
      for(chnl_ in chnls_)
      {
        # converting MFI of current plate of current channel into a vector needed for regression analysis
        
        y_ = if(drctn_ == 'column'){ matrix(data = MFI_mat[[chnl_]], ncol = 1, nrow = length(MFI_mat[[chnl_]]), byrow = F)[,,drop = T] }else      # reading matrix column-wise
                                   { matrix(data = MFI_mat[[chnl_]], ncol = 1, nrow = length(MFI_mat[[chnl_]]), byrow = T)[,,drop = T] }          # reading matrix row-wise
        x_ = 1:length(y_)                                                                                                                         # x_ is then offset of y_
  
        # removing outliers from plate MFIs before regression using inter-quantile range
  
        IQR_ = IQR(y_)
        quartiles_ = quantile(y_, probs = c(.25, .75))
        lowWhisker_ = max(min(y_), quartiles_[1] - IQR_*1.5)
        upWhisker_ = min(max(y_), quartiles_[2] + IQR_*1.5)
        y_tmp = y_[lowWhisker_ < y_ & y_ < upWhisker_]      # valid MFIs of current plate of current channel
  
        # fitting regression line
  
        fit_ = lm(data = data.frame(offset = 1:length(y_tmp), chnl = y_tmp), formula = chnl~offset)     # fitting simple linear regression
        a_ = fit_$coefficients["offset"]                                                                # slope of the regression line
        if(a_ <= 0)
        {
          warning(paste0('Slope for plate #',plate_,' channel #',chnl_,' was not positive, no intra-plate correction is necessary!'))     # if slope was not positive, no correction
          next()                                                                                                                          # goes to next channel
        }
        b_ = fit_$coefficients["(Intercept)"]                                                           # intercept of regression line
        alpha_ = b_/(a_*x_ + b_)                                                                        # correction coefficients for each well per channel: y*/y = b/ax+b -> y* = y(b/ax+b), y* is translated y
  
        # intra-plate correction of all files of current plate of current channel
  
        for(offset_ in 1:length(fcs_dt_ls))     # fcs files in fcs_dt_ls are already ordered according to their offset in MFI matrix (MFI_mat)
        {
          if(is.null(fcs_dt_ls[[offset_]])) { next() }      # some wells could be empty on the plate! It goes to next well
  
          # correcting MFI of current fcs for current channel (to avoid recomputing IQR)
          
          if(drctn_ == 'column')      # extracting (i,j) from column-wise offset
          {
            j_ = round(ceiling(offset_/m_))     # round is added because arithmetic on integer and double is not exact in R
            i_ = round(offset_ - m_*(j_-1))
          }else
          {
            i_ = round(ceiling(offset_/n_))
            j_ = round(offset_ - n_*(i_-1))
          }
          MFI_mat[[chnl_]][i_,j_] = MFI_mat[[chnl_]][i_,j_]*alpha_[offset_]     # y*/y = b/ax+b -> y* = y(b/ax+b)
          
          # correcting fcs file
          
          fcs_dt_ls[[offset_]]@exprs[,chnl_] = fcs_dt_ls[[offset_]]@exprs[,chnl_]*alpha_[offset_]     # Basically we should have used y* = y + (y*m - ym) for translation of cell/event values along with their MFI for current channel where ym is y-coordinate of the MFI and y* is translated y.
                                                                                                      # However, since this formula would translate some values to x+y- quadrant, we used the same way we translated MFI values, i.e. translation by proportion not addition.
        }
      }
  
      # writing intra-plate corrected fcs files temporarily in "data" directory
  
      for(offset_ in 1:length(fcs_dt_ls))
      {
        if(is.null(fcs_dt_ls[[offset_]])) { next() }
        write.FCS(x = fcs_dt_ls[[offset_]], filename = paste0(inURL,fcs_flNms[offset_],'.fcs'))     # matrix of events in this well
      }
    }
    rm(fcs_dt_ls)
    MFI_mats[[plate_]] = MFI_mat
  }
  
  if(FITPLOT)
  {
    # plotting regression plots for each channel per page
    
    if(CORRECT) { pdf(file = paste0(outURL,'/MFI_fit_intra-plate_corrected.pdf')) }else{ pdf(file = paste0(outURL,'/MFI_fit_intra-plate_no_correction.pdf')) }
    
    for(plate_ in plates_)
    {
      for (chnl_ in chnls_)
      {
        # converting MFI of current plate of current channel into a vector needed for regression analysis
        
        y_ = if(drctn_ == 'column'){ matrix(data = MFI_mats[[plate_]][[chnl_]], ncol = 1, nrow = length(MFI_mats[[plate_]][[chnl_]]), byrow = F)[,,drop = T]/1e4 }else      # reading matrix column-wise
                                   { matrix(data = MFI_mats[[plate_]][[chnl_]], ncol = 1, nrow = length(MFI_mats[[plate_]][[chnl_]]), byrow = T)[,,drop = T]/1e4 }          # reading matrix row-wise
        x_ = 1:length(y_)                                                                                                                                                   # x_ then is offset of y_
        
        # removing outliers from plate MFIs before regression using inter-quantile range
        
        IQR_ = IQR(y_)
        quartiles_ = quantile(y_, probs = c(.25, .75))
        lowWhisker_ = max(min(y_), quartiles_[1] - IQR_*1.5)
        upWhisker_ = min(max(y_), quartiles_[2] + IQR_*1.5)
        y_tmp = y_[lowWhisker_ < y_ & y_ < upWhisker_]      # valid MFIs of current plate of current channel needed for regression analysis
        
        # fitting regression line
        
        fit_ = lm(data = data.frame(offset = 1:length(y_tmp), chnl = y_tmp), formula = chnl~offset)     # fitting simple linear regression
        a_ = if(CORRECT){ 0 }else{ fit_$coefficients["offset"] }                                        # slope of the regression line
        b_ = fit_$coefficients["(Intercept)"]                                                           # intercept of the regression line
  
        p_ = ggplot(data = data.frame(x = x_, y = y_, Plate = plate_), aes(x = x_, y = y_, color = as.factor(Plate)))+
                theme(axis.line = element_line(color = 'black'), panel.background = element_blank(), legend.key = element_rect(fill = NA),
                      text = element_text(face = 'bold',size = 20), plot.title = element_text(vjust = 0.5, hjust = 0.5), aspect.ratio = 1)+
                guides(color = guide_legend(override.aes = list(size = 5)))+
                labs(x = 'Pseudo-time', y = expression('MFI ('%*%10^-4*')'), title = chnl_, color = 'Plate')+
                geom_point(pch = 20, size = 2)+
                scale_color_manual(values = 'red')+
                geom_abline(intercept = b_, slope = a_, color = 'darkblue')+
                scale_x_continuous(expand = c(0.01, 0.01)) + scale_y_continuous(expand = c(0.01, 0.01))
        plot(p_)
      }
    }
    graphics.off()
  }
  
  # STEP 3: Correcting for inter-plate signal drift ####
  
  if(CORRECT)
  {
    message('inter-plate correction')
    
    # computing inter-plate correction coefficient
    
    old_intcpts = matrix(data = NA, nrow = length(plates_), ncol = length(chnls_))      # a matrix to store current intercepts per plate and per channel
    rownames(old_intcpts) = plates_
    colnames(old_intcpts) = chnls_
    new_intcpts = list()                                                                # a list to store new intercepts (baselines) after correction per channel
    for(chnl_ in chnls_)                                                                # for each channel
    {
      dt_ = numeric()
      for(plate_ in plates_)      # for each plate
      {
        dt_[plate_] = median(MFI_mats[[plate_]][[chnl_]])     # current plate's median for current channel
        old_intcpts[plate_,chnl_] = dt_[plate_]
      }
      new_intcpts[[chnl_]] = median(dt_)                      # median of plate medians for current channel
    }
    
    # inter-plate correction
    
    for(plate_ in plates_)      # for each plate
    {
      # reading fcs files
      
      annot_tmp = annot_[annot_$plate %in% plate_,]     # annotations related to wells of current plate
      rows_ = sort(unique(annot_tmp$row))               # number of rows in current plate
      cols_ = sort(unique(annot_tmp$column))            # number of columns in current plate
      for(rw_ in 1:nrow(annot_tmp))                     # for each well on the current plate
      {
        fcs_dt = read.FCS(filename = paste0(inURL,annot_tmp$file[rw_],'.fcs'), transformation = F, truncate_max_range = F)      # event data of current well
        
        # extracting row and column indices of current well
        
        i_ = which(rows_ == annot_tmp$row[rw_])         # row index of current well
        j_ = which(cols_ == annot_tmp$column[rw_])      # column index of current well
        
        # inter-plate correction per channel
        
        for(chnl_ in chnls_)      # for each channel
        {
          # computing correction coefficient
          
          alpha_ = new_intcpts[[chnl_]]/old_intcpts[plate_,chnl_]     # correction coefficient: b*/b = y*/y -> y* = y(b*/b)
  
          # correcting MFI of current fcs (to avoid recomputing IQR)
          
          MFI_mats[[plate_]][[chnl_]][i_,j_] = MFI_mats[[plate_]][[chnl_]][i_,j_]*alpha_     # b*/b = y*/y -> y* = y(b*/b)
          
          # correcting fcs file of current well for current channel
          
          fcs_dt@exprs[,chnl_] = fcs_dt@exprs[,chnl_]*alpha_      # Basically we should have used y* = y + (y*m - ym) for translation of cell/event values along with their MFI for current channel where ym is y-coordinate of the MFI and y* is translated y.
                                                                  # However, since this formula would translate some values to x+y- quadrant, we used the same way we translated MFI values, i.e. translation by proportion not addition.
        }
        
        # writing corrected fcs files
        
        keyword(fcs_dt)[['$FIL']] = paste0(annot_tmp$file[rw_],'_compensated_corrected')      # updating $FIL keyword
        write.FCS(x = fcs_dt, filename = paste0(inURL,annot_tmp$file[rw_],'.fcs'))            # matrix of cells/events in this well
      }
      
      # writing MFI matrices
      
      for(chnl_ in chnls_)
      {
        write.table(x = MFI_mats[[plate_]][[chnl_]], file = paste0(outURL,'/MFI_mat_P',plate_,chnl_,'.txt'), quote = F, sep = '\t', row.names = T, col.names = T)
      }
    }
  }
  
  # plotting fit lines
  
  if(FITPLOT)
  {
    dt_ = list()                # a matrix with these columns: MFI, plate and channel
    k_ = 1
    for(plate_ in plates_)      # for each plate
    {
      for (chnl_ in chnls_)     # for each channel
      {
        tmp = if(drctn_ == 'column'){ matrix(data = MFI_mats[[plate_]][[chnl_]], ncol = 1, nrow = length(MFI_mats[[plate_]][[chnl_]]), byrow = F)[,,drop = T]/1e4 }else     # reading MFI matrix column-wise
                                    { matrix(data = MFI_mats[[plate_]][[chnl_]], ncol = 1, nrow = length(MFI_mats[[plate_]][[chnl_]]), byrow = T)[,,drop = T]/1e4 }         # reading MFI matrix row-wise
        dt_[[k_]] = data.frame(MFI = tmp, plate = as.character(plate_), channel = chnl_, stringsAsFactors = T)
        k_ = k_ + 1
      }
    }
    dt_ = do.call(what = rbind, args = dt_)
    
    # plotting all dot plots on one PDF page per channel
    
    if(CORRECT) { pdf(file = paste0(outURL,'/MFI_fit_inter-plate_corrected.pdf')) }else{ pdf(file = paste0(outURL,'/MFI_fit_inter-plate_no_correction.pdf')) }
    for(chnl_ in chnls_)      # for each channel 
    {
      dt_tmp = dt_[dt_$channel %in% chnl_,]
      p_ = ggplot(data = dt_tmp, aes(x = 1:nrow(dt_tmp), y = MFI, color = plate))+
              theme(axis.line = element_line(color = 'black'), panel.background = element_blank(), legend.key = element_rect(fill = NA),
                    text = element_text(face = 'bold',size = 20), plot.title = element_text(vjust = 0.5, hjust = 0.5), aspect.ratio = 1)+
              guides(color = guide_legend(override.aes = list(size = 5)))+
              labs(x = 'Pseudo-time', y = expression('MFI ('%*%10^-4*')'), title = chnl_, color = 'Plate')+
              geom_point(pch = 20, size = 2)+
              geom_hline(yintercept = median(dt_tmp$MFI), color = 'darkred')+
              scale_x_continuous(expand = c(0.01, 0.01)) + scale_y_continuous(expand = c(0.01, 0.01))
      plot(p_)
    }
    graphics.off()
  }
  
  # STEP 4: Plotting plate heatmaps ####
  
  if(HEATPLOT)
  {
    numOfPlates = length(plates_)                 # number of plates in this assay
    nrow_ = floor(sqrt(numOfPlates))              # number of rows per page
    ncol_ = ceiling(numOfPlates/nrow_)            # number of columns per page
    width_ = height_ = max(nrow_, ncol_)*8.45     # every heatmap takes up this width/height in inches
    
    if(!CORRECT)      # if no correction was requested
    {
      fl_ = paste0(outURL,'/plates_heatmap_no_correction.pdf')
      pdf(file = fl_, width = width_, height = height_)

      for(chnl_ in chnls_)      # plots plates for each channel on one PDF page
      {
        plot_list = list()
        
        # finding max across all plates
        
        max_ = 0
        for(plate_ in 1:numOfPlates)
        {
          max_tmp = max(MFI_mats[[plate_]][[chnl_]]/1e4)
          if(max_ < max_tmp) { max_ = max_tmp }
        }
        
        for(plate_ in 1:numOfPlates)      # plotting each plate for current channel
        {
          dt_mat = MFI_mats[[plate_]][[chnl_]]/1e4
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
            if(u1_ == 0) { u1_ = NULL}      # to avoid redundant legend key on the heatmap
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
                        main = paste0('Plate: ',plate_,', Channel: ', chnl_),
                        fontsize = 7,
                        border_color = 'grey90', color = cols_,
                        breaks = c(min(dt_)-col_step,dt_),     # in pheatmap color intervals (showing lower and uper bounds) are open on the left and closed on the right
                        legend = T,
                        legend_breaks = round(c(0,u1_,u2_, u3_, u4_, u5_),2),
                        annotation_col = NA, silent = T)
          plot_list[[plate_]] = p_[[4]]
        }
        grid.arrange(arrangeGrob(grobs= plot_list, nrow = nrow_, ncol = ncol_))
      }
      graphics.off()
      
    }else     # if correction was requested
    {
      fl_ = paste0(outURL,'/plates_heatmap_correction.pdf')
      pdf(file = fl_, width = width_, height = height_)
      
      for(chnl_ in chnls_)      # plots plates for each channel on one PDF page
      {
        plot_list = list()
    
        # finding max across all plates
        max_ = 0
        dt_ = NULL
        for(plate_ in 1:numOfPlates)
        {
          dt_ = c(dt_, as.numeric(MFI_mats[[plate_]][[chnl_]]/1e4))     # palte plate_ showing MFI of channel chnl_
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
          if(u1_ == 0) { u1_ = 0}     # to avoid redundant legend key on heatmap
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
    
        for(plate_ in 1:numOfPlates)
        {
          LEGEND_SHOW = F
          if(plate_ == 1) { LEGEND_SHOW = T }
          
          dt_mat = MFI_mats[[plate_]][[chnl_]]/1e4      # plate plate_ showing MFI of channel chnl_
          dt_mat = rbind(dt_mat, MAX = c(max_, rep(0,ncol(dt_mat)-1)))
          
          p_ = pheatmap(mat = dt_mat,
                        cluster_cols = F,
                        cluster_rows = F,
                        show_rownames = T, cellheight = 20,
                        show_colnames = T, cellwidth = 20,
                        main = paste0('Plate: ',plate_,', Channel: ', chnl_),
                        fontsize = 7,
                        border_color = 'grey90', color = cols_,
                        breaks = c(min(dt_)-col_step,dt_),     # in pheatmap color intervals (showing lower and uper bounds) are open on the left and closed on the right
                        legend = LEGEND_SHOW,
                        legend_breaks = round(c(0,u1_,u2_, u3_, u4_, u5_),2),
                        annotation_col = NA, silent = T)
          plot_list[[plate_]] = p_[[4]]
        }
        grid.arrange(arrangeGrob(grobs= plot_list, nrow = nrow_, ncol = ncol_))
      }
      graphics.off()
    }
  }

  return(NULL)
}
