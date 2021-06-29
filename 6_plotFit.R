# dat %>%
#   bind_cols(mat_fitted) %>%
#   ggplot(aes(x=Time, y=Record, group=Sub, color='a')) +
#   geom_line() +
#   facet_wrap(~Sub)+
#   geom_line(aes(y=y_hatA1, color='c')) + # fits with all basis
#   geom_line(aes(y=y_hatA2, color='d1')) + # first half of basis
#   geom_line(aes(y=y_hatA3, color='d2')) + # second half of basis
#   geom_line(aes(y=y_hatA4, color='d3')) +
#   theme_bw() +
#   scale_colour_manual(
#     name = NULL,
#     values = c(
#       'a'='black',
#       'c'='red',
#       'd1'='green',
#       'd2'= 'blue',
#       'd3'='purple',
#       'd4'='pink'),
#     labels = c(
#       'a'='observation',
#       'c'=sprintf('fitted all %d basis',nFourier),
#       'd1'='fitted low frequency basis',
#       'd2'='fitted medium low frequency basis',
#       'd3'='fitted high frequency basis')

# grid plot
fitted.plot <- function(dat,mat_fitted,colNames, plotColors, legendLabels, plotRow,xText,yText, nPerGroup){
# plot the fitted curve of first nPerGroup subjects in each group

  dat1 <- data.frame(dat,mat_fitted)

  sub <-  as.numeric(dat$ID)
  dat2 <- dat %>%
    select(ID,Group) %>%
    unique %>%
    group_by(Group) %>%
    slice(1:4) %>%
    ungroup
  subUniq <- dat2 %>%
    select(ID) %>%
    unlist
  nSub <- length(subUniq)
  plotCol <- ceiling(nSub/plotRow)
  subTitles <- dat2$Group

  x_rg <- range(dat1[,'t'])
  x_axis <- pretty(x_rg)
  x_axis <- x_axis[between(x_axis, x_rg[1],x_rg[2])]
  x_axis <- x_axis[-length(x_axis)]
  y_rg <- range(dat1[,colNames])
  y_axis <- pretty(y_rg)
  y_axis <- y_axis[between(y_axis, y_rg[1],y_rg[2])]


  pushViewport(plotViewport(c(5,3,.5,.5)))

  grid.text(yText, rot=90, x= unit(-2.5, 'lines'))
  grid.text(xText, y= unit(-3, 'lines'))
  grid.legend(
    vp=viewport( y= unit(-4.5, 'lines')),
    labels= legendLabels, nrow=1, ncol= length(legendLabels),
    do.lines= T, lines.first= F,
    gp= gpar(col= plotColors, lwd=2),

  )
  pushViewport(viewport(layout=grid.layout(nrow=plotRow,ncol=plotCol)))

  for(t in seq(nSub)){
    pRow <- ceiling(t/plotCol)
    pCol <- map_dbl(t%%plotCol, function(x) ifelse(x==0, plotCol, x))

    pushViewport(viewport(layout.pos.row=pRow,layout.pos.col=pCol))
    pushViewport(plotViewport(c(0,0, 0,0),
                              xscale=range(x_rg), yscale=y_rg))

    grid.rect()
    grid.text(
      subTitles[t],
      x= unit(.5,'lines'), y=unit(1, 'npc')-unit(.5,'lines'), just='left')

    idx <- which(sub==subUniq[t])

    for(i in seq(colNames)){
      grid.lines(
        dat1[idx,'t'], dat1[idx, colNames[i]],
        default.units='native',
        gp= gpar(col= plotColors[i]))
    }

    if(pRow==plotRow) grid.xaxis(x_axis,x_axis)
    if(pCol==1) grid.yaxis(y_axis,y_axis)

    popViewport()
    popViewport()
  }

  popViewport()
  popViewport()
}
