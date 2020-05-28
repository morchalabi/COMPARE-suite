require(flowCore)
require(pheatmap)
require(gridExtra)
require(ggplot2)

args_ = commandArgs(trailingOnly = T)

###### Manual Debugging ########
# args_ = c('-correct', 'T',
#           '--fit-plot', 'T',
#           '--heat-plot', 'T')
################################

# STEP 0: Options control ####

INVALID = F

CORRECT = as.logical(args_[which(grepl(x = args_, pattern = '^-correct', fixed = F))+1])
CORRECT = CORRECT[!is.na(CORRECT)]
if(length(CORRECT) == 0) { INVALID = T }

FITPLOT = as.logical(args_[which(grepl(x = args_, pattern = '^--fit-plot', fixed = F))+1])
FITPLOT = FITPLOT[!is.na(FITPLOT)]
if(length(FITPLOT) == 0) { INVALID = T }

HEATPLOT = as.logical(args_[which(grepl(x = args_, pattern = '^--heat-plot', fixed = F))+1])
HEATPLOT = HEATPLOT[!is.na(HEATPLOT)]
if(length(HEATPLOT) == 0) { INVALID = T }

if(INVALID)
{
  message('\nInvalid call. Usage:\n',
          'Rscript 4_cytotoxicity_correction.R \\\n',
          '-correct TRUE \\\n',
          '--fit-plot TRUE \\\n',
          '--heat-plot TRUE\n')
  quit(save = 'no')
}

message('You set:',
        '\ncytotoxicity drift correction to: ', CORRECT,
        '\nplotting regressed lines to: ',FITPLOT,
        '\nplotting plate heatmap to:',   HEATPLOT,'\n')

options(nwarnings = 10000)      # shows all warnings (default is last 50)

# STEP 1: Computing toxicities ####

# reading in annotation file

annot_ = read.table(file = '../data/Annotations.txt',header = T, sep = '\t',as.is = T, check.names = F, stringsAsFactors = F)
plates_ = sort(unique(annot_$plate))

# reading fcs files of each plate

toxicity_mats = list()
for(plate_ in plates_)
{
  # Reading fcs files ####
  message('Reading fcs files of plate #',plate_,':')
  
  annot_tmp = annot_[annot_$plate %in% plate_,]
  rows_ = sort(unique(annot_tmp$row))
  m_ = length(rows_)
  cols_ = sort(unique(annot_tmp$column))
  
  fcs_flNms = list()
  
  toxicity_mat = matrix(data = 0, nrow = length(rows_), ncol = length(cols_))
  rownames(toxicity_mat) = rows_
  colnames(toxicity_mat) = cols_
  for(rw_ in 1:nrow(annot_tmp))
  {
    fcs_dt = read.FCS(filename = paste0('../data/',annot_tmp$file[rw_],'.fcs'), transformation = F)

    # computing offset from beginning of plate matrix
    
    i_ = which(rows_ == annot_tmp$row[rw_])
    j_ = which(cols_ == annot_tmp$column[rw_])
    offset_ = m_*(j_-1)+i_                          # column-wise offset
    toxicity_mat[i_,j_] = nrow(fcs_dt@exprs)        # computing toxicities
    fcs_flNms[[offset_]] = annot_tmp$file[rw_]      # fcs file name in this offset of plate
  }

  # STEP 2: Intra-plate correcting fcs file ####
  
  if(CORRECT)
  {
    message('intra-plate correction of plate #',plate_,':')

    y_ = matrix(data = toxicity_mat, ncol = 1, nrow = length(toxicity_mat), byrow = F)[,,drop = T]      # reading matrix columnwise
    x_ = 1:length(y_)                                                                                   # x_ then is offset of y_

    # removing outliers

    IQR_ = IQR(y_)
    quartiles_ = quantile(y_, probs = c(.25, .75))
    lowWhisker_ = max(min(y_), quartiles_[1] - IQR_*1.5)
    upWhisker_ = min(max(y_), quartiles_[2] + IQR_*1.5)
    y_tmp = y_[lowWhisker_ < y_ & y_ < upWhisker_]

    # fitting regressed line

    fit_ = lm(data = data.frame(offset = 1:length(y_tmp), live_cells = y_tmp), formula = live_cells~offset)
    a_ = fit_$coefficients["offset"]      # slope
    if(a_ <= 0)     # if slope is non-positive, no correction is carried out
    {
      warning(paste0('Slope for plate #',plate_,' was negative, no intra-plate correction is performed!'))
      toxicity_mats[[plate_]] = toxicity_mat
      next()
    }
    b_ = fit_$coefficients["(Intercept)"]
    alpha_ = (a_*x_)/(a_*x_ + b_)         # correction factors

    # step 2.2: correction

    for(offset_ in 1:length(fcs_flNms))      # fcs files in fcs_flNms are already ordered by their offset in toxicity matrix (toxicity_mat)
    {
      if(is.null(fcs_flNms[[offset_]])) { next() }

      # correcting toxicities of current fcs
      p_ = offset_/m_
      i_ = round((p_-ceiling(p_)+1)*m_)     # round is added because arithmatic on integer and double is not exact in R
      j_ = round(ceiling(p_))
      toxicity_mat[i_,j_] = toxicity_mat[i_,j_]*(1-alpha_[offset_])
    }

    # step 2.3: plotting fit lines if requested

    if(FITPLOT)
    {
      pdf(file = paste0('../out/toxicity_fit_plate',plate_,'.pdf'))
      
      # plotting uncorrected toxicities

      plot(x = x_, y = y_, col = 'red', main = 'Uncorrected', xlab = 'Pseudo-time', ylab = 'Live cells',
           ylim = range(y_), font = 2, font.lab = 2, pch = 20, cex = 1, cex.lab = 1.4, cex.axis = 1.5)
      abline(a = b_, b = a_, col = 'blue')

      # plotting corrected toxicities

      y_ = matrix(data = toxicity_mat, ncol = 1, nrow = length(toxicity_mat), byrow = F)[,,drop = T]      # corrected toxicities
      plot(x = x_, y = y_, col = 'red', main = 'Corrected', xlab = 'Pseudo-time', ylab = 'Live cells',
           ylim = range(y_), font = 2, font.lab = 2, pch = 20, cex = 1, cex.lab = 1.4, cex.axis = 1.5)
      abline(h = b_, col = 'blue')
      
      graphics.off()
    }
  }
  toxicity_mats[[plate_]] = toxicity_mat
}

# STEP 3: Correcting for inter-plate batch effect ####

if(CORRECT)
{
  message('inter-plate correction')
  
  # computing inter-plate correction coefficient
  old_intcpts = double()
  for(plate_ in plates_) { old_intcpts[plate_] = median(toxicity_mats[[plate_]]) }
  new_intcpts = median(old_intcpts)
  
  for(plate_ in plates_)
  {
    # Reading fcs files ####
    
    annot_tmp = annot_[annot_$plate %in% plate_,]
    rows_ = sort(unique(annot_tmp$row))
    cols_ = sort(unique(annot_tmp$column))
    for(rw_ in 1:nrow(annot_tmp))
    {
      fcs_dt = read.FCS(filename = paste0('../data/',annot_tmp$file[rw_],'.fcs'), transformation = F)     # matrix of events in this well
      
      # computing offset from beginning of plate matrix
      
      i_ = which(rows_ == annot_tmp$row[rw_])
      j_ = which(cols_ == annot_tmp$column[rw_])
      
      # inter-plate correction
      
      # computing correction factors
      alpha_ = (new_intcpts - old_intcpts[plate_])/old_intcpts[plate_]     # correction factors (fold changes)

      # correcting toxicities of currect fcs
      toxicity_mats[[plate_]][i_,j_] = toxicity_mats[[plate_]][i_,j_]*(1+alpha_)

      # correcting fcs file
      keyword(fcs_dt) = list(Toxicity = toxicity_mats[[plate_]][i_,j_])
      
      # writing corrected fcs files
      write.FCS(x = fcs_dt, filename = paste0('../data/',annot_tmp$file[rw_],'.fcs'))      # matrix of events in this well
    }
    
    # Writing toxicity matrices ####
    write.table(x = toxicity_mats[[plate_]], file = paste0('../out/toxicity_mat_P',plate_,'.txt'), quote = F, sep = '\t', row.names = T, col.names = T)
  }
  
  # plotting fit lines ####
  
  if(FITPLOT)
  {
    dt_ = list()
    for(plate_ in plates_)
    {
      tmp = matrix(data = toxicity_mats[[plate_]], ncol = 1, nrow = length(toxicity_mats[[plate_]]), byrow = F)
      dt_[[plate_]] = data.frame(toxicity = tmp, plate = as.character(plate_), stringsAsFactors = T)
    }
    dt_ = do.call(what = rbind, args = dt_)
    
    pdf(file = '../out/toxicity_fit_all_plates.pdf')
    p_ = ggplot(data = dt_, aes(x = 1:nrow(dt_), y = toxicity/1e4, color = plate))+
           theme(axis.line = element_line(color = 'black'), panel.background = element_blank(), legend.key = element_rect(fill = NA),
                  text = element_text(face = 'bold',size = 20), plot.title = element_text(vjust = 0.5, hjust = 0.5), aspect.ratio = 1)+
           guides(color = guide_legend(override.aes = list(size = 5)))+
           labs(x = 'Pseudo-time', y = expression(bold('Live cells ('%*%10^-4*')')), title = 'Inter-plate correction', color = 'Plate')+
           geom_point(pch = 20, size = 2)+
           geom_hline(yintercept = new_intcpts/1e4, color = 'darkred')+
           scale_x_continuous(expand = c(0.01, 0.01)) + scale_y_continuous(expand = c(0.01, 0.01))
    plot(p_)
    graphics.off()
  }
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
    fl_ = paste0('../out/plates_heatmap_no_correction_cyto.pdf')
    pdf(file = fl_, width = width_, height = height_)
  
    plot_list = list()
    
    # finding max in all plates
    max_ = 0
    for(plate_ in 1:numOfPlates)
    {
      max_tmp = max(toxicity_mats[[plate_]]/1e4)     # palte plate_ showing toxicity
      if(max_ < max_tmp) { max_ = max_tmp }
    }
    
    for(plate_ in 1:numOfPlates)
    {
      dt_mat = toxicity_mats[[plate_]]/1e4      # palte plate_ showing toxicity
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
    
  }else
  {
    fl_ = paste0('../out/plates_heatmap_correction_cyto.pdf')
    pdf(file = fl_, width = width_, height = height_)
  
    plot_list = list()

    # finding max in all plates
    max_ = 0
    dt_ = NULL
    for(plate_ in 1:numOfPlates)
    {
      dt_ = c(dt_, as.numeric(toxicity_mats[[plate_]]/1e4))     # palte plate_ showing toxicity
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

    for(plate_ in 1:numOfPlates)
    {
      LEGEND_SHOW = F
      if(plate_ == 1) { LEGEND_SHOW = T }
      
      dt_mat = toxicity_mats[[plate_]]/1e4      # palte plate_ showing toxicity
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

message('\nDone!\n')
summary(warnings())
