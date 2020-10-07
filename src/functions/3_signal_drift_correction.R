
require(flowCore)
require(pheatmap)
require(gridExtra)
require(ggplot2)

# chnls_
# CORRECT
# FITPLOT
# HEATPLOT
# inURL
# outURL
#
step3_signal_drift_correction = function(chnls_, CORRECT, FITPLOT, HEATPLOT, inURL = '../../data/', outURL = '../../out/')
{
  # STEP 1: Computing MFIs ####
  
  # reading in annotation file
  
  annot_ = read.table(file = paste0(inURL,'/Annotations.txt'), header = T, sep = '\t',as.is = T, check.names = F, stringsAsFactors = F)
  plates_ = sort(unique(annot_$plate))
  
  # reading fcs files of each plate
  
  MFI_mats = list()
  for(plate_ in plates_)
  {
    # reading fcs files
    message('Reading fcs files of plate #',plate_,':')
  
    annot_tmp = annot_[annot_$plate %in% plate_,]
    rows_ = sort(unique(annot_tmp$row))
    m_ = length(rows_)
    cols_ = sort(unique(annot_tmp$column))
  
    fcs_dt = list()
    MFI_mat = list()
    fcs_flNms = character()
    for(chnl_ in chnls_)
    {
      mat_ = matrix(data = 0, nrow = length(rows_), ncol = length(cols_))
      rownames(mat_) = rows_
      colnames(mat_) = cols_
      MFI_mat[[chnl_]] = mat_
    }
    for(rw_ in 1:nrow(annot_tmp))
    {
      fcs_dt_tmp = read.FCS(filename = paste0(inURL,annot_tmp$file[rw_],'.fcs'), transformation = F, truncate_max_range = F)
  
      # computing offset from beginning of plate matrix
      i_ = which(rows_ == annot_tmp$row[rw_])
      j_ = which(cols_ == annot_tmp$column[rw_])
      offset_ = m_*(j_-1)+i_                        # column-wise offset
      fcs_dt[[offset_]] = fcs_dt_tmp
      fcs_flNms[offset_] = annot_tmp$file[rw_]      # fcs file name in this offset of plate
  
      # computing MFIs: because here we are not assigning negatives to 0, it is fine to use IQR for removing outliers
      for(chnl_ in chnls_)
      {
        # MFIs
        dt_tmp = fcs_dt[[offset_]]@exprs[,chnl_, drop = T]
        dt_tmp = dt_tmp[which(0 <= dt_tmp)]                           # N.B.: which igonores NA and NaNs. Non-positives are always non-positives even in presence of drift
        IQR_ = IQR(dt_tmp)                                            # inter-quantile region
        quartiles_ = quantile(dt_tmp, probs = c(.25, .75))            # 25th and 75th percentiles
        lowWhisker_ = max(min(dt_tmp), quartiles_[1] - IQR_*1.5)      # lower whisker
        upWhisker_ = min(max(dt_tmp), quartiles_[2] + IQR_*1.5)
        MFI_mat[[chnl_]][i_,j_] = median(dt_tmp[lowWhisker_ < dt_tmp & dt_tmp < upWhisker_])
      }
    }
  
    # STEP 2: Intra-plate correcting fcs file ####
    
    if(CORRECT)
    {
      message('intra-plate correction of plate #',plate_,':')
  
      # computing intra-plate correction coefficients for each channels of current plate!
  
      for(chnl_ in chnls_)
      {
        y_ = matrix(data = MFI_mat[[chnl_]], ncol = 1, nrow = length(MFI_mat[[chnl_]]), byrow = F)[,,drop = T]      # reading matrix columnwise
        x_ = 1:length(y_)                                                                                           # x_ then is offset of y_
  
        # removing outliers
  
        IQR_ = IQR(y_)
        quartiles_ = quantile(y_, probs = c(.25, .75))
        lowWhisker_ = max(min(y_), quartiles_[1] - IQR_*1.5)
        upWhisker_ = min(max(y_), quartiles_[2] + IQR_*1.5)
        y_tmp = y_[lowWhisker_ < y_ & y_ < upWhisker_]
  
        # fitting regressed line
  
        fit_ = lm(data = data.frame(offset = 1:length(y_tmp), chnl = y_tmp), formula = chnl~offset)
        a_ = fit_$coefficients["offset"]      # slope
        if(a_ <= 0) { warning(paste0('Slope for plate #',plate_,' channel #',chnl_,' was not positive, no intra-plate correction is necessary!')); next() }
        b_ = fit_$coefficients["(Intercept)"]
        alpha_ = (a_*x_)/(a_*x_ + b_)         # correction factors
  
        # correction
  
        for(offset_ in 1:length(fcs_dt))                # fcs files in fcs_dt are already ordered according to their offset in MFI matrix (MFI_mat)
        {
          if(is.null(fcs_dt[[offset_]])) { next() }     # some wells could be empty on the plate!
  
          # correcting fcs file
          fcs_dt[[offset_]]@exprs[,chnl_] = fcs_dt[[offset_]]@exprs[,chnl_]*(1-alpha_[offset_])
  
          # correcting MFIs of current fcs (to avoid recomputing IQR)
          p_ = offset_/m_
          i_ = round((p_-ceiling(p_)+1)*m_)             # round is added because arithmatic on integer and double is not exact in R
          j_ = round(ceiling(p_))
          MFI_mat[[chnl_]][i_,j_] = MFI_mat[[chnl_]][i_,j_]*(1-alpha_[offset_])
        }
      }
  
      # writing intra-plate corrected fcs files temporarily in "data" directory
  
      for(offset_ in 1:length(fcs_dt))
      {
        if(is.null(fcs_dt[[offset_]])) { next() }
        write.FCS(x = fcs_dt[[offset_]], filename = paste0(inURL,fcs_flNms[offset_],'.fcs'))      # matrix of events in this well
      }
    }
    rm(fcs_dt)
    MFI_mats[[plate_]] = MFI_mat
  }
  
  if(FITPLOT)
  {
    if(CORRECT) { pdf(file = paste0(outURL,'/MFI_fit_intra-plate_corrected.pdf')) }else{ pdf(file = paste0(outURL,'/MFI_fit_intra-plate_no_correction.pdf')) }
    
    for(plate_ in plates_)
    {
      for (chnl_ in chnls_)
      {
        y_ = matrix(data = MFI_mats[[plate_]][[chnl_]], ncol = 1, nrow = length(MFI_mats[[plate_]][[chnl_]]), byrow = F)[,,drop = T]/1e4
        x_ = 1:length(y_)                                                                                           # x_ then is offset of y_
        
        # removing outliers
        IQR_ = IQR(y_)
        quartiles_ = quantile(y_, probs = c(.25, .75))
        lowWhisker_ = max(min(y_), quartiles_[1] - IQR_*1.5)
        upWhisker_ = min(max(y_), quartiles_[2] + IQR_*1.5)
        y_tmp = y_[lowWhisker_ < y_ & y_ < upWhisker_]
        
        # fitting regressed line
        
        fit_ = lm(data = data.frame(offset = 1:length(y_tmp), chnl = y_tmp), formula = chnl~offset)
        a_ = if(CORRECT) { 0 }else { fit_$coefficients["offset"] }     # slope
        b_ = fit_$coefficients["(Intercept)"]
  
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
  
  # STEP 3: Correcting for inter-plate batch effect ####
  
  if(CORRECT)
  {
    message('inter-plate correction')
    
    # computing inter-plate correction coefficient
    new_intcpts = list()
    old_intcpts = matrix(data = NA, nrow = length(plates_), ncol = length(chnls_))
    rownames(old_intcpts) = plates_
    colnames(old_intcpts) = chnls_
    for(chnl_ in chnls_)
    {
      dt_ = numeric()
      for(plate_ in plates_)
      {
        dt_[plate_] = median(MFI_mats[[plate_]][[chnl_]])
        old_intcpts[plate_,chnl_] = dt_[plate_]
      }
      new_intcpts[[chnl_]] = median(dt_)
    }
    
    for(plate_ in plates_)
    {
      # reading fcs files
      
      annot_tmp = annot_[annot_$plate %in% plate_,]
      rows_ = sort(unique(annot_tmp$row))
      cols_ = sort(unique(annot_tmp$column))
      for(rw_ in 1:nrow(annot_tmp))
      {
        fcs_dt = read.FCS(filename = paste0(inURL,annot_tmp$file[rw_],'.fcs'), transformation = F, truncate_max_range = F)
        
        # computing offset from beginning of plate matrix
        i_ = which(rows_ == annot_tmp$row[rw_])
        j_ = which(cols_ == annot_tmp$column[rw_])
        
        # inter-plate correction
        
        for(chnl_ in chnls_)
        {
          # computing correction factors
          alpha_ = (new_intcpts[[chnl_]] - old_intcpts[plate_,chnl_])/old_intcpts[plate_,chnl_]     # correction factors (fold changes)
  
          # correcting fcs file
          fcs_dt@exprs[,chnl_] = fcs_dt@exprs[,chnl_]*(1+alpha_)
  
          # correcting MFIs of currect fcs (to avoid recomputing IQR)
          MFI_mats[[plate_]][[chnl_]][i_,j_] = MFI_mats[[plate_]][[chnl_]][i_,j_]*(1+alpha_)
        }
        
        # writing corrected fcs files
        keyword(fcs_dt)[['$FIL']] = paste0(annot_tmp$file[rw_],'_compensated_corrected')
        write.FCS(x = fcs_dt, filename = paste0(inURL,annot_tmp$file[rw_],'.fcs'))      # matrix of events in this well
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
    dt_ = list()
    k_ = 1
    for(plate_ in plates_)
    {
      for (chnl_ in chnls_)
      {
        tmp = matrix(data = MFI_mats[[plate_]][[chnl_]], ncol = 1, nrow = length(MFI_mats[[plate_]][[chnl_]]), byrow = F)[,,drop = T]/1e4
        dt_[[k_]] = data.frame(MFI = tmp, plate = as.character(plate_), channel = chnl_, stringsAsFactors = T)
        k_ = k_ + 1
      }
    }
    dt_ = do.call(what = rbind, args = dt_)
    
    if(CORRECT) { pdf(file = paste0(outURL,'/MFI_fit_inter-plate_corrected.pdf')) }else{ pdf(file = paste0(outURL,'/MFI_fit_inter-plate_no_correction.pdf')) }
    for(chnl_ in chnls_)
    {
      dt_tmp = dt_[dt_$channel %in% chnl_,]
      p_ = ggplot(data = dt_tmp, aes(x = 1:nrow(dt_tmp), y = MFI, color = plate))+
              theme(axis.line = element_line(color = 'black'), panel.background = element_blank(), legend.key = element_rect(fill = NA),
                    text = element_text(face = 'bold',size = 20), plot.title = element_text(vjust = 0.5, hjust = 0.5), aspect.ratio = 1)+
              guides(color = guide_legend(override.aes = list(size = 5)))+
              labs(x = 'Pseudo-time', y= expression('MFI ('%*%10^-4*')'), title = chnl_, color = 'Plate')+
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
    numOfPlates = length(plates_)
    nrow_ = floor(sqrt(numOfPlates))                          # number of rows per page
    ncol_ = ceiling(numOfPlates/floor(sqrt(numOfPlates)))     # number of columns per page
    width_ = height_ = max(nrow_, ncol_)*8.45                 # every heatmap takes up this width/height in inches
    
    if(!CORRECT)
    {
      fl_ = paste0(outURL,'/plates_heatmap_no_correction.pdf')
      pdf(file = fl_, width = width_, height = height_)
      for(chnl_ in chnls_)
      {
        plot_list = list()
        
        # finding max in all plates
        max_ = 0
        for(plate_ in 1:numOfPlates)
        {
          max_tmp = max(MFI_mats[[plate_]][[chnl_]]/1e4)     # palte plate_ showing MFI of channel chnl_
          if(max_ < max_tmp) { max_ = max_tmp }
        }
        
        for(plate_ in 1:numOfPlates)
        {
          dt_mat = MFI_mats[[plate_]][[chnl_]]/1e4      # palte plate_ showing MFI of channel chnl_
          dt_mat = rbind(dt_mat, MAX = c(max_, rep(0,ncol(dt_mat)-1)))
          
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
      
    }else
    {
      fl_ = paste0(outURL,'/plates_heatmap_correction.pdf')
      pdf(file = fl_, width = width_, height = height_)
      for(chnl_ in chnls_)
      {
        plot_list = list()
    
        # finding max in all plates
        max_ = 0
        dt_ = NULL
        for(plate_ in 1:numOfPlates)
        {
          dt_ = c(dt_, as.numeric(MFI_mats[[plate_]][[chnl_]]/1e4))     # palte plate_ showing MFI of channel chnl_
        }
        max_ = max(dt_)
        
        ###
        
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
          
          dt_mat = MFI_mats[[plate_]][[chnl_]]/1e4      # palte plate_ showing MFI of channel chnl_
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
