# # base R plot
# plot(Kcluster.obj$KL,
#      xlab='K: number of clusters',
#      ylab='KL',
#      main=sprintf('Case %d: fourier basis n from %d to %d\n select K = %d at threshold of %.1f', t, min(idxA),max(idxA), nCluster, thKL))
# abline(h=threshold, col='red')

# grid plot
kmeans_fitted.plot <- function(ls_sB,plotRow, subTitles,xText,yText){
  tb <- ls_sB$sBtb
  cutoff <- ls_sB$cutoff
  nClusters <- ls_sB$nClusters
  nPlot <- length(nClusters)
  plotCol <- nPlot/plotRow
  nCl <- nrow(tb)+1

  pushViewport(plotViewport(c(1,2,.5,1)))
  pushViewport(viewport(layout=grid.layout(nrow=plotRow,ncol=plotCol)))

  for(t in seq(nPlot)){
    pRow <- ceiling(t/plotCol)
    pCol <-  map_dbl(t%%plotCol, function(x) ifelse(x==0, plotCol, x))

    threshold <- cutoff[t]
    x_axis <- pretty(1:nCl)
    x_axis[1] <- 1
    x_rg <- range(x_axis)
    y_axis <- pretty(tb[,t],n=4)
    y_rg <- range(y_axis)

    pushViewport(viewport(layout.pos.row=pRow,layout.pos.col=pCol))
    pushViewport(plotViewport(c(3,3,3,1),
                              xscale=range(x_rg), yscale=y_rg))

    grid.text(sprintf('(%s) %s',letters[t], subTitles[t]),
              x= unit(-1,'lines'), y=unit(1, 'npc')+unit(2,'lines'), just='left')

    grid.points(2:nCl, tb[,t], default.units = 'native', size = unit(.5, "char"))
    grid.points(nClusters[t], tb[nClusters[t]-1,t], default.units = 'native', size = unit(.5, "char"), pch=19, gp=gpar(col='red'))
    grid.lines(x = x_rg,y = rep(threshold,2), default.units = 'native', gp= gpar(col='red',lty=2))
    grid.xaxis(x_axis,x_axis)
    grid.yaxis(y_axis,y_axis)
    if(pCol==1) grid.text(yText, rot=90, x= unit(-4, 'lines'))
    if(pRow==plotRow) grid.text(xText, y= unit(-3, 'lines'))

    popViewport()
    popViewport()
  }

  popViewport()
  popViewport()
}
