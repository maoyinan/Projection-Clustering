ellipse.fun <- function(xLeft,xRight, gridStep=.01){
  ellCenter <- (xLeft+xRight)/2
  ellR <- (xRight-xLeft)/2
  grid <- seq(xLeft,xRight,gridStep)
  y <- 1/ellR*sqrt(ellR^2-(grid-ellCenter)^2)
  list(grid=grid,y=y)
}
place.lower <- function(y,w=1){
  unit(-y, 'native') + unit(-w,'lines')
}
prop.alpha <- function(prop, cutPoints){
  group <- cut(prop, breaks=sort(c(0,1,cutPoints)), include.lowest = T)
  allAlpha <- numeric(3)
  allAlpha[1] <- 0
  allAlpha[2] <- .2
  allAlpha[3] <- 1
  al <- allAlpha[as.numeric(group)]
  al
}
cluster.plot <- function(ls_prob,df_groups,ls_groups,cutPoints, fileName, subTitles, allColors, seed,thin, ls_labels){

  nSub <- nrow(ls_prob[[1]])
  x_axis <- seq(nSub)
  x_rg <- range(x_axis)
  if(is.null(ls_labels)) ls_labels <- rep(list(unique(df_groups$Group)),length(ls_prob))
  df_groups0 <- df_groups

  pushViewport(plotViewport(c(.5,.5,0,.5)))
  pushViewport(viewport(layout=grid.layout(nrow=length(ls_prob),ncol=1)))

  for(t in seq_along(ls_prob)){
    set.seed(seed)
    tb <- ls_prob[[t]]
    tb1 <- data.frame(as.matrix(tb,dimnames=list(1:nSub,1:nSub))) %>%
      bind_cols(row=1:nSub) %>%
      pivot_longer(-row,names_to='col',values_to='prob') %>%
      mutate(col=as.numeric(gsub('X','',col))) %>%
      filter(!is.na(prob)) %>%
      mutate(index=1:n()) %>%
      filter(index%in%sample(n(),thin*n(),replace=F)) %>%
      select(-index)

    pushViewport(viewport(layout.pos.row=t,layout.pos.col=1))
    pushViewport(plotViewport(c(4,0,2,0),
                              xscale=range(x_rg), yscale=c(0,1)))

    grid.text(sprintf('(%s) %s',letters[t], subTitles[t]),x= 0, y=unit(1, 'npc')+unit(1,'lines'), just='left')

    # x axis of subjects and label the groups in rectangle boxes
    groups <- ls_labels[[t]]
    df_groups <- data.frame(Group=groups) %>%
      left_join(df_groups0) %>%
      mutate(Sub1=1:n())
    for(i in seq_along(groups)){
      cat(groups[i],'\n')
      x <- df_groups %>% filter(Group==groups[i]) %>% .$Sub1
      grid.rect(
        x=mean(x), y=0,
        width = diff(range(x)), height= unit(1,'lines'),
        just = 'top',default.units = 'native')
      grid.text(label = groups[i], x=unit(mean(x),'native'), y=unit(-.5,'lines'))
    }
    {
      for(i in seq(nrow(tb1))){
        sub1 <- df_groups %>%
          filter(Sub==tb1$row[i]) %>%
          .$Sub1
        sub2 <- df_groups %>%
          filter(Sub==tb1$col[i]) %>%
          .$Sub1
        gi <- df_groups %>% filter(Sub1==sub1) %>% .$Group
        gj <- df_groups %>% filter(Sub1==sub2) %>% .$Group
        pij <- tb1$prob[i]
        probij <- prop.alpha(pij,cutPoints)
        if(probij>0){
          xLeft <- min(c(sub1,sub2))
          xRight <- max(c(sub1,sub2))
          ell <- ellipse.fun(xLeft,xRight,gridStep=.01)
          if(gi==gj)
            grid.lines(x = ell$grid, y= ell$y, default.units = 'native', gp= gpar(col=allColors[t], alpha=probij))
          else {
            if(!is.null(ls_groups)){
              if (gi%in%ls_groups[[gj]])
                grid.lines(x = ell$grid, y= place.lower(ell$y,w=1), default.units = 'native', gp= gpar(col=allColors[t], alpha=probij))
              else
                grid.lines(x = ell$grid, y= place.lower(ell$y,w=1), default.units = 'native', gp= gpar(col='black', alpha=probij))
            }
            else
              grid.lines(x = ell$grid, y= place.lower(ell$y,w=1), default.units = 'native', gp= gpar(col='black', alpha=probij))
          }
        }
      }
    }
    popViewport()
    popViewport()
  }
  popViewport()
  popViewport()
}
