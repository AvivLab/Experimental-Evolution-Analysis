# Biolog
TalPlot <- function(eigenvalues, cut){
  sigma= eigenvalues
  n = length(eigenvalues)
  sigmak = log(c(1,sigma))
  sigma = log(c(sigma,1))
  st = sigmak-sigma
  std = st[1:n+1]
  pcs =c(1:(n))
  df = data.frame(PC = pcs, Logk = std)
  df = df[-n,]
  ggline(df, x = "PC",y="Logk",color = "royalblue",xlab = "Principal Component", ylab = TeX("log $\\frac{\\sigma_k}{\\sigma_{k+1}}$")
  ) + geom_vline(xintercept = cut,linetype="dotted", size =1)
}

# From somewhere else
save_pheatmap_pdf <- function(x, filename, width=1200, height=1000, res = 150) {
  pdf(filename)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


# Theme
theme_xp <- function () { 
  theme_minimal(base_size=14, base_family="") %+replace%
    theme(
      plot.title = element_text(face = "bold",
                                size = rel(1.2), hjust = 0.5),
      text = element_text(),
      panel.background = element_rect(colour = NA),
      #plot.background = element_rect(colour = NA),
      # panel.border = element_rect(colour = NA),
      axis.title = element_text(face = "bold",size = rel(1)),
      axis.title.y = element_text(angle=90,vjust =2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(), 
      axis.line = element_line(colour="black"),
      axis.ticks = element_line(),
      panel.grid.major = element_line(colour="#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.size= unit(0.2, "cm"),
      legend.margin = unit(0, "cm"),
      legend.title = element_text(face="italic"),
      plot.margin=unit(c(10,5,5,5),"mm"),
      strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
      strip.text = element_text(face="bold")
    )
}

# Color mutations
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  Non_Synonymous = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#A6CEE3", col = NA))
  },
  Small_indel = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#FB9A99", col = NA))
  },
  Intergenic_snp = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = '#FF7F00', col = NA))
  },
  Stop = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = '#6A3D9A', col = NA))
  },
  Synonymous = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#FDBF6F", col = NA))
  },
  Large_deletion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = '#E31A1C', col = NA))
  },
  Large_amplification = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#1F78B4", col = NA))
  }
)
col = c("Non_Synonymous" = "#A6CEE3", "Small_indel" = "#FB9A99", "Intergenic_snp" = '#FF7F00',
        'Stop'='#6A3D9A','Synonymous'='#FDBF6F','Large_deletion'='#E31A1C',"Large_amplification"="#1F78B4")