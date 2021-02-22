## Marcin Imielinski
## The Broad Institute of MIT and Harvard / Cancer program.
## marcin@broadinstitute.org
##
## Weill-Cornell Medical College
## mai9037@med.cornell.edu
##
## New York Genome Center
## mimielinski@nygenome.org
##
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.

## You should have received a copy of the GNU Lesser General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

                                        #
                                        # General utility functions
                                        #

.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
        fn(base::get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
    {
        if (all(class(x) == 'Hits'))
            c(S4Vectors::queryLength(x), S4Vectors::subjectLength(x))
        else
            as.numeric(dim(x))[1:2]
    }))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, function(x) c(muffle(length(x)), NA)[1])[vec]    

    out <- data.table(name = names, obj.type, obj.size, prettyMem(obj.size), obj.dim)
    setnames(out, c("Name", "Type", "Size", "Pretty", "Rows", "Columns"))
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}


prettyMem = function(x, places = 3)
{
    power <- pmin(6, floor(log(abs(x), 1000)))
    units <- c("B", "kB", "MB", "GB", "TB", "PB", "EB")[power+1]
    x <- x/(1000^power)
    paste(prettyNum(signif(x,places)), units)
}


## shorthand listing largest objects in the workspace

#' @name lsos
#' @title lsos
#'
#' @description
#' returns largest object in workspace
#'
#' @param n num of objects to return
#' @author Stack Overflow Post 21442
#' @export
lsos <- function(..., n=10) {
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=FALSE)
}

#' @name du
#' @title du
#'
#' @description
#' du of directories one folder deep
#'
#' @param path folder
#' @param d max depth
#' @export
du = function(path, d = 1)
{
    if (!file.exists(path))
        stop('path does not exist')
    udf = fread(cmd = sprintf("cd %s; du -h --max-depth %s", path, d), sep = '\t', header = FALSE)
    setnames(udf, c("size", "path"))
    units = gsub('.*([A-Z])$', '\\1', udf$size)
    units = ifelse(grepl('[A-Z]', udf$size), units, 'B')
    units = c('K' = 1000, 'P' = 1e15, 'T' = 1e12, 'B' = 1, 'G' = 1e9, 'M' = 1e6)[units]
    sizenum = as.numeric(gsub('[A-Z]', '', udf$size))
    udf[, MB := sizenum*units/1e6]
    udf = udf[rev(order(MB)), ]
    udf[, path := gsub('^\\.\\/', '', path)]
    return(udf)
}


#' @name ldim
#' @title ldim
#'
#' @description
#' returns dimensions of all objects contained in list
#'
#' @param l list
#' @author Marcin Imielinski
#' @export
ldim = function(l)
{
    return(lapply(l, function(x) {if (!is.null(dim(x))) dim(x) else length(x)}))
}


#' @name file.name
#' @title file.name
#'
#' @description
#' parses filenames from character vector of paths
#'
#' @param paths character vector of full paths
#' @return character vector of just file names
#' @author Marcin Imielinski
#' @export
########
file.name = function(paths)
{
    return(gsub('(^|(.*\\/))?([^\\/]*)', '\\3', paths))
}

#' @name file.dir
#' @title file.dir
#'
#' @description
#' grabs file.dirs from list of paths
#'
#' @name file.dir
#' @param paths character vector of full paths
#' @return character vector of just file.names
#' @author Marcin Imielinski
#' @export
file.dir = function(paths)
{
    return(gsub('(^|(.*\\/))?([^\\/]*)$', '\\2', paths))
}

########
#' splits a single string according to fixed widths contained in fw (ie each components i of fw denotes the width of field i in string str
#'
########
strsplit.fwf = function(str, fw)
{
    if (length(str)>1)
    {
        warning('String should be of length 1, only taking first element')
        str = str[1];
    }

    cs = cumsum(fw);
    return(substr(rep(str, length(fw)), cs-fw+1,c(cs[1:(length(cs)-1)], nchar(str))))
}

########
                                        # Returns vector of line counts for each file in path
########
line.counts = function(paths)
{
    out = rep(NA, length(paths))
    ix = which(file.exists(paths))
    out[ix] = sapply(paths, function(x) { p = pipe(paste('cat ', x, ' | wc -l ')); as.numeric(readLines(p)); close(p)});
    return(out)
}



##############
#' @name border
#' @title border
#' @description
#' orders rows of a logical / binary matrix treating each row as binary number with digits encoded as TRUE / FALSE values of entries
#'
#' @name border
#' @param B input matrix logical format, or castable to logical
#' @param na.rm removes NA
#' @return B resorted using "binary" order
#' @export
##############
border = function(B, na.rm = TRUE)
{
    B = array(as.logical(B), dim = dim(B))
    tmp = vector(mode = "numeric", length = nrow(B));
    if (na.rm)
        B[is.na(B)] = FALSE;
    for (i in 1:ncol(B))
        tmp = tmp + 2^(ncol(B)-i)*as.numeric(B[,i]==1);
    return(order(tmp))
}


##############
#' @name fready
#' @title fread with name cleaning
#' @description
#' calls fread while cleaning names using provided or default pattern and sub
#'
#' @name fready
#' @param pattern character (default \\W)
#' @param sub character to sub in names (default _)
#' @return data.table
#' @export
##############
fready = function(..., pattern = "\\W+", sub = "_")
{
    tab = fread(...)
    nms = dedup(gsub(pattern, sub, names(tab), perl = TRUE), suffix = '.') %>% gsub('^[^A-Za-z]', '', ., perl = TRUE)
    setnames(tab, nms)
    return(tab)
}



################################
#' @name dedup
#' @title dedup
#'
#' @description
#' relabels duplicates in a character vector with .1, .2, .3
#' (where "." can be replaced by any user specified suffix)
#'
#' @param x input vector to dedup
#' @param suffix suffix separator to use before adding integer for dups in x
#' @return length(x) vector of input + suffix separator + integer for dups and no suffix for "originals"
#' @author Marcin Imielinski
#' @export
################################
dedup = function(x, suffix = '.')
{
    dup = duplicated(x);
    udup = setdiff(unique(x[dup]), NA)
    udup.ix = lapply(udup, function(y) which(x==y))
    udup.suffices = lapply(udup.ix, function(y) c('', paste(suffix, 2:length(y), sep = '')))
    out = x;
    out[unlist(udup.ix)] = paste(out[unlist(udup.ix)], unlist(udup.suffices), sep = '');
    return(out)
}


#' @name nona
#' @title nona
#'
#' @description
#'
#' Given data.frame / data.table dt outputs only the columns that are non.na 
#' in at least thresh fraction of the rows
#' 
#' @param dt
#' @param thresh
#' @return data.table with a subset of columns that have at least a given fraction of non NA entries
#' @author Marcin Imielinski
#' @export
################################
nona = function(dt, thresh = 0.1)
{
    na = colSums(is.na(dt))/nrow(dt)
    if (inherits(dt, 'data.table'))
        dt[, which(!na), with = FALSE]
    else
        dt[, which(!na), with = FALSE]
}



#' @name qq_pval
#' @title qq plot given input p values
#' @param obs vector of pvalues to plot, names of obs can be intepreted as labels
#' @param highlight optional arg specifying indices of data points to highlight (ie color red)
#' @param samp integer, optional specifying how many samples to draw from input data (default NULL)
#' @param lwd integer, optional, specifying thickness of line fit to data
#' @param pch integer dot type for scatter plot
#' @param cex integer dot size for scatter plot
#' @param conf.lines logical, optional, whether to draw 95 percent confidence interval lines around x-y line
#' @param max numeric, optional, threshold to max the input p values
#' @param repel logical flag whether to use ggrepel
#' @param label length(obs) charater vector of labels (NULL)
#' @param plotly toggles between creating a pdf (FALSE) or an interactive html widget (TRUE)
#' @param annotations named list of vectors containing information to present as hover text (html widget), must be in same order as obs input
#' @param gradient named list that contains one vector that color codes points based on value, must bein same order as obs input
#' @param titleText title for plotly (html) graph only
#' @author Marcin Imielinski, Eran Hodis, Zoran Gajic
#' @export
qq_pval = function(obs, highlight = c(), exp = NULL, lwd = 1, bestfit=T, col = NULL, col.bg='black', pch=18, cex=1, conf.lines=FALSE, max=NULL, max.x = NULL, max.y = NULL, qvalues=NULL, label = NULL, repel = FALSE, plotly = FALSE, annotations = list(), gradient = list(), titleText = "", subsample = NA, ...)
{
  if(!(plotly))
  {
    is.exp.null = is.null(exp)

    if (is.null(col))
      col = rep('black', length(obs))

    ix1 = !is.na(obs)
    if (!is.null(exp))
      if (length(exp) != length(obs))
        stop('length of exp must be = length(obs)')
      else
        ix1 = ix1 & !is.na(exp)

    if (is.null(highlight))
      highlight = rep(FALSE, length(obs))

    if (is.null(label))
      label = rep('', length(label))

    else if (is.logical(highlight))
    {
      if (length(highlight) != length(obs))
        stop('highlight must be either logical vector of same length as obs or a vector of indices')
    }
    else
      highlight = 1:length(obs) %in% highlight

    obs = -log10(obs[ix1])
    col = col[ix1]
    highlight = highlight[ix1]
    label = label[ix1]
    
    if (!is.null(exp))
      exp = -log10(exp[ix1])

    ix2 = !is.infinite(obs)
    if (!is.null(exp))
      ix2 = ix2 &  !is.infinite(exp)

    obs = obs[ix2]
    col = col[ix2]
    highlight = highlight[ix2]
    label = label[ix2]

    if (!is.null(exp))
      exp = exp[ix2]

    N <- length(obs)
    ## create the null distribution
    ## (-log10 of the uniform)

    if (is.null(exp))
      exp <- -log(1:N/N,10)
    else
      exp = sort(exp)

    if (is.null(max))
      max = max(obs,exp) + 0.5

    if (!is.null(max) & is.null(max.x))
      max.x = max

    if (!is.null(max) & is.null(max.y))
      max.y  = max

    if (is.null(max.x))
      max.x <- max(obs,exp) + 0.5

    if (is.null(max.y))
      max.y <- max(obs,exp) + 0.5

    if (is.exp.null)
    {
      tmp.exp = rev(seq(0, 7, 0.01))
      ix = 10^(-tmp.exp)*N
      c95 <-  qbeta(0.975,ix,N-ix+1)
      c05 <-  qbeta(0.025,ix,N-ix+1)

      if (conf.lines){
        ## plot the two confidence lines
        plot(tmp.exp, -log(c95,10), ylim=c(0,max.y), xlim=c(0,max.x), type="l", axes=FALSE, xlab="", ylab="")
        par(new=T)
        plot(tmp.exp, -log(c05,10), ylim=c(0,max.y), xlim=c(0,max.x), type="l", axes=FALSE, xlab="", ylab="")
        par(new=T)

        p1 <- rep(tmp.exp[1], 2)
        p2 <- c(-log(c95,10)[1], -log(c05,10)[1])

        lines(x=p1, y=p2)
        x.coords <- c(tmp.exp,rev(tmp.exp))
        y.coords <- c(-log(c95,10),rev(-log(c05,10)))
        polygon(x.coords, y.coords, col='light gray', border=NA)
        par(new=T)
      }
    }

    ord = order(obs)

                                        #colors = vector(mode = "character", length = length(obs)); colors[] = "black";

    colors = col
    colors[highlight] = "red";

    dat = data.table(x = sort(exp), y = obs[ord], colors = colors[ord], label = label[ord], pch = pch, cex = cex)

    if (!is.null(names(obs)))
    {
      dat$names = names(obs[ord])
      setkey(dat, names)
    }

    if (nrow(dat)>1e5) ## rough guide to subsampling the lower p value part of the plot
      subsample = 5e4/nrow(dat)

    lambda = lm(y ~ x-1, dat)$coefficients;

    if (is.na(subsample[1]))
      dat[, plot(x, y, xlab = expression(Expected -log[10](italic(P))), ylab = expression(Observed -log[10](italic(P))), xlim = c(0, max.x), col = colors, ylim = c(0, max.y), pch=pch, cex=cex, bg=col.bg, ...)]
    else
    {
      subsample = pmin(pmax(0, subsample[1]), 1)
      dat[ifelse(x<=2, ifelse(runif(length(x))<subsample, TRUE, FALSE), TRUE), plot(x, y, xlab = expression(Expected -log[10](italic(P))), ylab = expression(Observed -log[10](italic(P))), xlim = c(0, max.), col = colors, ylim = c(0, max.y), pch=pch, cex=cex, bg=col.bg, ...)]
    }
    
    if (!is.null(dat$label) && any(nchar(dat$label)>0, na.rm = TRUE))
    {
      dat[nchar(label)>0, text(x, y, labels=label, pos=3)];
    }
    
    lines(x=c(0, max(max.y, max.x)), y = c(0, max(max.x, max.y)), col = "black", lwd = lwd)
    
    if (!is.na(subsample))
      dat = dat[sample(nrow(dat), subsample*nrow(dat)), ]


    lines(x=c(0, max.x), y = c(0, lambda*max.x), col = "red", lty = 2, lwd = lwd);
    legend('bottomright',sprintf('lambda=\n %.2f', lambda), text.col='red', bty='n')
  }
  else{

                                        #browser()

        if(length(annotations) < 1){
            hover <- do.call(cbind.data.frame, list(p = obs))
        }
        else{
            hover <- do.call(cbind.data.frame, list(annotations, p = obs))
        }
        hover <- as.data.table(hover)


        is.exp.null = is.null(exp)
        if (is.null(col))
            col = "black"
        ix1 = !is.na(hover$p)
        if (!is.null(exp))
            if (length(exp) != length(hover$p))
                stop("length of exp must be = length(hover$obs)")
            else ix1 = ix1 & !is.na(exp)
        if (is.null(highlight))
            highlight = rep(FALSE, length(hover$p))
        else if (is.logical(highlight)) {
            if (length(highlight) != length(hover$p))
                stop("highlight must be either logical vector of same length as obs or a vector of indices")
        }
        else highlight = 1:length(hover$p) %in% highlight
        hover$obs = -log10(hover$p[ix1])
        hover = hover[ix1]
        highlight = highlight[ix1]
        if (!is.null(exp))
            exp = -log10(exp[ix1])
        ix2 = !is.infinite(hover$obs)
        if (!is.null(exp))
            ix2 = ix2 & !is.infinite(exp)
        hover = hover[ix2]
        highlight = highlight[ix2]
        if (!is.null(exp))
            exp = exp[ix2]
        N <- length(hover$obs)
        if (is.null(exp))
            exp <- -log(1:N/N, 10)
        else exp = sort(exp)
        if (is.null(max))
            max <- max(hover$obs, exp) + 0.5
        else max <- max
        if (is.exp.null) {
            tmp.exp = rev(seq(0, 7, 0.01))
            ix = 10^(-tmp.exp) * N
            c95 <- qbeta(0.975, ix, N - ix + 1)
            c05 <- qbeta(0.025, ix, N - ix + 1)
            if (FALSE) {   ##Don't need if not using conf.line (might put this in the future)
                plot(tmp.exp, -log(c95, 10), ylim = c(0, max), xlim = c(0, max),
                     type = "l", axes = FALSE, xlab = "", ylab = "")

                par(new = T)
                plot(tmp.exp, -log(c05, 10), ylim = c(0, max), xlim = c(0, max),
                     type = "l", axes = FALSE, xlab = "", ylab = "")

                par(new = T)
                p1 <- rep(tmp.exp[1], 2)
                p2 <- c(-log(c95, 10)[1], -log(c05, 10)[1])
                lines(x = p1, y = p2)
                x.coords <- c(tmp.exp, rev(tmp.exp))
                y.coords <- c(-log(c95, 10), rev(-log(c05, 10)))
                polygon(x.coords, y.coords, col = "light gray", border = NA)
                par(new = T)
            }
        }

                                        #creating the ploting data.table (dat) and organizing the annotations to create hover text
        ord = order(hover$obs)
        hover = hover[ord]
        dat = hover
        hover$obs = NULL

                                        #Creating the hover text
        if(length(colnames(hover)) > 1){
            annotation_names  = sapply(colnames(hover), paste0, " : ")
            annotation_names_wLineBreak  = paste("<br>", annotation_names[2:length(annotation_names)],
                                                 sep = "")
            annotation_names = c(annotation_names[1], annotation_names_wLineBreak)
        }
        else{
            annotation_names  = sapply(colnames(hover), paste0, " : ")
        }

                                        #Checking if there is a gradient and if so adding it to the plotting data.table (dat)
        gradient_control = FALSE
        if(length(gradient )!= 0){
            dat$grad = gradient[[1]][ord]
            gradient_control = TRUE
        }
        else {
            dat$grad = c()
        }


        dat$x = sort(exp)
        dat$y = dat$obs

                                        #declare so we can use in If statement
        p <- NULL

                                        #hacky subsampling but works really well, just maxing out the number of points at 8k
                                        #and removing the extra from the non-sig
                                        #(looks to be -logp of 2.6 here can make this more dynamic later )

        if (nrow(dat) <=  8000){

            dat4 = dat
            dat4$obs = NULL
            dat4$x = NULL
            dat4$y = NULL
            dat4$grad = NULL

            trans = t(dat4)
            hover_text = c()
            for (i in 1:dim(trans)[2]){
                outstr = paste(c(rbind(annotation_names, trans[,i])), sep = "", collapse = "")
                hover_text = c(hover_text,outstr)
            }
                                        #            browser()
            if(gradient_control){
                p <- dat[, plot_ly(data = dat, x=x, y=y, hoverinfo = "text",text = hover_text, color = grad,
                                   colors = c("blue2","gold"),marker = list(colorbar = list(title = names(gradient[1]))),
                                   mode = "markers",type = 'scatter')
                         %>% layout(xaxis = list(title = "<i>Expected -log<sub>10</sub>(P)</i>"),
                                    yaxis = list(title = "<i>Observed -log<sub>10</sub>(P)</i>")) ]
            }
            else{
                p <- dat[, plot_ly(data = dat, x=x, y=y, hoverinfo = "text",text = hover_text,
                                   mode = "markers",type = 'scatter')
                         %>% layout(xaxis = list(title = "<i>Expected -log<sub>10</sub>(P)</i>"),
                                    yaxis = list(title = "<i>Observed -log<sub>10</sub>(P)</i>")) ]
            }
        }


        else {

            dat$ID = c(1:nrow(dat))
            dat2 = dat[ y < 2.6,]
            dat3 = as.data.frame(dat2)
            dat3 = as.data.table(dat3[ sample(nrow(dat3), min(4000,nrow(dat3))), ])
            dat2 = rbind(dat3,dat[!(ID%in%dat2$ID),])
            dat2$ID = NULL

            dat4 = dat2
            dat4$obs = NULL
            dat4$x = NULL
            dat4$y = NULL
            dat4$grad = NULL

            trans = t(dat4)
            hover_text = c()
            for (i in 1:dim(trans)[2]){
                outstr = paste(c(rbind(annotation_names, trans[,i])), sep = "", collapse = "")
                hover_text = c(hover_text,outstr)
            }

            if(gradient_control){
                p <- dat2[, plot_ly(data = dat2, x=x, y=y,hoverinfo = "text", text = hover_text, color = grad,
                                    colors = c("blue2","gold"),marker = list(colorbar = list(title = names(gradient[1]))),
                                    mode = "markers",type = 'scatter')
                          %>% layout(xaxis = list(title = "<i>Expected -log<sub>10</sub>(P)</i>"),
                                     yaxis = list(title = "<i>Observed -log<sub>10</sub>(P)</i>")) ]
            }
            else{
                p <- dat2[,  plot_ly(data = dat2, x=x, y=y,hoverinfo = "text", text = hover_text,
                                     mode = "markers",type = 'scatter')
                          %>% layout(xaxis = list(title = "<i>Expected -log<sub>10</sub>(P)</i>"),
                                     yaxis = list(title = "<i>Observed -log<sub>10</sub>(P)</i>")) ]
            }

        }

                                        #       browser()

                                        #Calculating lambda, Note that this is using the whole data set not the subsampled one
        lambda = lm(y ~ x - 1, dat)$coefficients
        lambda_max = max*as.numeric(lambda)


        ##adding shapes (lines) + title  note that html <b></b> style is used for mods and plotting lines
        ##is done by specifying two points on the line (x0/y0 and x1/y1)
        p <- layout(p,title = sprintf("<b>%s</b>" ,titleText),titlefont = list(size = 24),
                    shapes = list(list(type = "line",line = list(color = 'black'),
                                       x0 = 0, x1  = max, xref = "x", y0 = 0, y1 = max,yref ="y"),
                                  list( type = "line", line = list(color = "red"),
                                       x0 = 0, x1 = max, xref = "x", y0 = 0, y1 = lambda_max, yref = "y")),
                    annotations = list(
                        x = (0.9 * max),
                        y = (0.03 * max),
                        text = paste("lambda =",sprintf("%.2f", signif(lambda,3)), collapse = " "),
                        font = list(
                            color = "red",
                            size = 20
                        ),
                        showarrow = FALSE,
                        xref = "x",
                        yref = "y"
                    ),
                    margin = list(
                        t = 100

                    ),
                    hovermode = "compare")
    }
}

#' @name qq_repel
#' @title QQ plot with repel
#'
#' @param ps vector of p values
#' @param label length(ps) vector of labels
#' @export
#' @author Ashley Doane
qq_repel <- function(ps, label = rep('', length(ps)), conf.lines = FALSE, ci = 0.95, print = TRUE) {
  N  <- length(ps)
  df <- data.table(
    p = ps,
    observed = -log10(ps),
    label = label)
  df = df[order(p), ][, ":="(
           expected = -log10(1:N / N),
           clower   = -log10(qbeta(ci,     1:N, N - 1:N + 1)),
           cupper   = -log10(qbeta(1 - ci, 1:N, N - 1:N + 1))
         )]
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))

  lambda = lm(observed ~ expected-1, df)$coefficients

  p= ggplot(df, aes(expected,observed,label=label)) +
    geom_point(aes(expected, observed), color = ifelse(df$label=="", "black", "red")) +
    geom_abline(intercept = 0,
                slope = 1,
                alpha = 0.5)

  if (conf.lines)
    {
     p = p + geom_line(aes(expected, cupper), linetype = 2) +
       geom_line(aes(expected, clower),
              linetype = 2,
              color = 'red') 
       }

  p = p + 
    xlab(log10Pe) +
    ylab(log10Po) +
                                        #coord_cartesian(xlim = c(0,5), ylim=c(0,25)) +
                                        #xlim(0, 10) +
    ggrepel::geom_text_repel(
      data = subset(df, (expected > 3 )),
      nudge_x      =  4.5 - subset(df, (expected > 3 ))$expected,
      direction    = "y",
      hjust        = 0,
      segment.size = 0.2) +
    ggrepel::geom_text_repel(
      data = subset(df, (expected <= 3 )),
      nudge_x      =  1.5 - subset(df, (expected <= 3 ))$expected,
      direction    = "y",
      hjust        = 1,
      segment.size = 0.2) +
    cowplot::theme_cowplot(font_size = 12) +
    annotate(geom = 'text', label = sprintf('lambda == %.2f', lambda), size = 6, x = max(df$expected), y = 0, hjust = 1, vjust = 0, parse = TRUE)
                                        #geom_text_repel(size=4,box.padding = 0.25,segment.size = .25,
                                        #   max.iter = 2000
                                        #)
  if (print)
    print(p)
  else
    p  
}


#' @name wfplot
#' @title Quick waterfall plot
#' @description Quick waterfall plot
#'
#' data is a numeric vector
#' labels are text labels of the same length as data
#' col is either (1) an unamed list of unique colors (2) a named list mapping unique labels to colors
#'
#' @param data length n numeric vector to be drawn and sorted on y axis
#' @param labels length n character vector categorical labels of data
#' @param names.arg length n character vector, optional, of individual labels to be drawn verticallyon x axis
#' @param col optional named character vector mapping unique category labels to colors
#' @param las optional integer vector specifying orientation of labels on barplot
#' @param cex numeric value specifying size of names.arg data labels
#' @param ... additional arguments to barplot
#' @param leg.pos NULL
#' @return plot
#' @author Marcin Imielinski
#' @export
wfplot = function(data, labels = NULL, names.arg = NULL, col = NULL, las = 2, cex = 1, leg.pos = NULL, ...)
{
    ix = order(data);
    labels = as.character(labels)
    ulab = unique(labels)

    par(mar = c(12.1, 8.1, 4.1, 2.1))


    if (is.null(col))
    {
        if (length(ulab)>2)
            col = RColorBrewer::brewer.pal(length(ulab), 'Set3')
        else
            col = c('gray', 'red')
    }

    if (is.null(names(col)))
    {
        col = col[match(labels, ulab)]
        names(col) = labels;
    }
    else
    {
        og.col = col;
        col = col[labels];
        ulab = intersect(names(og.col), names(col))
    }

    barplot(data[ix], col = col[ix], names.arg = names.arg[ix], las = las, cex.names = 0.5*cex, border = FALSE, ...)

    if (is.null(leg.pos))
        leg.pos = c(mean(par('usr')[1:2])/4, 3*data[which.max(abs(data))]/4)

    legend(leg.pos[1], leg.pos[2], legend = ulab, fill = col[ulab])
}





#' @name list.expr
#' @title list.expr
#' @description
#'
#' Takes a character or numeric vector and makes an expression for re-creating that character in source code
#'
#' @name list.expr
#' @author Marcin Imielinski
#' @param x input vector
#' @return character vector of command to create the input vector
#' @export
############
list.expr = function(x)
{
    y = ''
    if (!is.null(names(x)))
        y = paste0('"', names(x), '"=')

    if (is.character(x) | is.factor(x))
        out =paste("c('", paste(y, x, sep = "", collapse = "', '"), "')", sep = "")
    else
        out = paste("c(", paste(y, x, sep = "", collapse = ", "), ")", sep = "")
    writeLines(out)
}

########
#' @name fuckr
#' @title fuckr
#'
#' @description
#' ... what you feel when R is getting on your nerves. Toggles options(error = ) to enable / disable debugging mode.
#'
#' toggles options error recover / NULL
#'
#' @export
#' @author Marcin Imielinski
########
fuckr = function()
{
    if (!is.null(options()$error))
    {
        options(error = NULL);
        print('Options error set to NULL');
    }
    else
    {
        options(error = recover);
        print('Options error set to recover');
    }
}

#################
## flatten
##
## flattens 3rd dim of 3D array along cdim 1 (ie rows) or cdim 2 (ie cols) pasting together the appropriate combinations of dimnames with sep "sep"
## or if sep = NULL, then just dropping the 3rd dimension names
##
#################
flatten = function(A, cdim = 2, sep = "_")
{
    if (!(cdim==1 | cdim ==2))
        stop('cdim must be 1 or 2')

    ind = order(rep(c(1:dim(A)[cdim]), dim(A)[3]));

    out = A[,,1];

    if (cdim == 2)
    {
        if (dim(A)[3]>1)
            for (i in 2:dim(A)[3])
                out = cbind(out, A[,,i]);
        dimnames(A)[[1]] = dimnames(A)[[1]]
    }

    if (cdim == 1)
    {
        if (dim(A)[3]>1)
            for (i in 2:dim(A)[3])
                out = rbind(out, A[,,i]);
        dimnames(A)[[2]] = dimnames(A)[[2]]
    }

    out = out[,ind]; #reshuffle to get desired ordering
    newdimnames = rep(dimnames(A)[[cdim]], each = dim(A)[3]);
    if (!is.null(sep))
        newdimnames = paste(newdimnames, dimnames(A)[[3]], sep = sep);
    dimnames(out)[[cdim]] = newdimnames;

    return(out)
}

##################
#' @name bsub_cmd
#' @title bsub_cmd
#' @description
#'
#' Makes bsub command that wraps shell command "cmd" to send to queue "queue"
#' redirebmccting output / error etc streams to path prefixed by "jname",
#' optional_args: maximum memory requirements "mem", "jlabel" job label
#'b
#' @param cmd length n  vector of shell commands, optionally named, one per job
#' @param queue optional length n or length 1 character specifying queue to send jobs to (default hour)
#' @param jname optional length n character specifying names of jobs, this will be the root of the output files generated by the job
#' @param jlabel optional length n character specifying labels of jobs, this the string
#' @param jgroup optional length n character specifying job group name
#' @param mem length n or length 1 integer specifying GB of memory to be used by jobs
#' @param group character specifying job group (default cgafolk)
#' @param cwd character specifying which working directory to launch jobs from (default is current working directory of R session)
#' @param mc.cores length n or 1 integer specifying how many cores to assign to each job
#' @param deadline logical flag whether to send jobs to deadline queue
#' @return character vector of bsub commands, which can run using system or dumped to a shell script
#' @export
#' @author Marcin Imielinski
##################
bsub_cmd = function(cmd, queue, jname = NULL, jlabel=NULL, jgroup = NULL, mem=NULL, group = "cgafolk", cwd = NULL, mc.cores = NULL, deadline = F)
{
    if (is.null(jname) & is.null(names(cmd)))
        jname = 'job'

    if (length(jname) != length(cmd))
        jname = rep(jname, length(cmd))

    if (!is.null(jname))
        names(cmd) = dedup(jname)

    qjname = paste( "\"", names(cmd), "\"", sep="" )
    qjout = paste( "\"", names(cmd), ".bsub.out", "\" ", sep="" )
    qjerr = paste( "\"", names(cmd), ".bsub.err", "\" ", sep="" )
    qjrout = paste( "\"", names(cmd), ".R.out", "\" ", sep="" )
    out_cmd = paste( "bsub -q ", queue, " -o ", qjout, " -e ",  qjerr, " -P ", group);
    if (!is.null(mem)) out_cmd = paste(out_cmd, " -R \"rusage[mem=", mem, "]\" ", sep = "");
    if (!is.null(jlabel)) out_cmd = paste(out_cmd, " -J ", jlabel )
    if (!is.null(jgroup)) out_cmd = paste(out_cmd, " -g ", sub('^\\/*', '/', jgroup))
    if (!is.null(cwd)) out_cmd = paste(out_cmd, " -cwd ", cwd )
    if (!is.null(mc.cores)) out_cmd = paste(out_cmd, sprintf(" -n %d,%d -R 'span[hosts=1]'", mc.cores, mc.cores))
    if (deadline) out_cmd = paste(out_cmd, '-sla DEADLINEsla')
    out_cmd = paste(out_cmd," \"",  cmd, "\"", sep = "")
    names(out_cmd)= names(cmd)
    return(out_cmd)
}


##############
#' @name query_lsf_out
#' @title query_lsf_out
#' @description
#'
#' parses "out" and "err" files of jobs with jname root to identify exit status and error codes of jobs
#' @param dir character specifyhing directory where .out and .err files are located
#' @param jname character vector names of jobs (as specified in bsub_cmd
#' @param detailed logical flag specifying whether to return "detailed" information
#' @param mc.cores integer specifying how many cores to use to parse the output data
#' @return data.frame of job info
#' @export
#' @author Marcin Imielinski
##############
lsf_query = lsf_out_query = query_lsf_out = function(dir = NULL, jname = NULL, detailed = F, mc.cores = 1)
                            {
                                if (!is.null(dir))
                                    dir = paste(dir, '/', sep = '')
                                else
                                {
                                    if (!is.null(jname))
                                    {
                                        dir = dirname(jname)
                                        jname = basename(jname)
                                    }
                                    else
                                        dir = ''
                                }

                                input.jname = jname
                                jname = gsub('\\.bsub\\.out$', '', gsub('\\.bsub\\.err$', '', jname))
                                names(input.jname) = jname

                                tmp.run = paste(normalizePath(dir), '/', jname, sep = '')

                                if (length(jname)==0)
                                    outs = data.frame(jname = NA,
                                                      out.file = NA,
                                                      err.file = NA,
                                                      run.file = NA,
                                                      exit_flag = NA, term_flag = NA, started = NA, reported = NA, hours_elapsed = NA, max_mem = NA, cpu_time = NA, stringsAsFactors = F)
                                else
                                {
                                    outs = data.frame(jname = gsub('\\.R$', '', jname),
                                                      out.file = paste(normalizePath(dir),'/', jname, '.bsub.out', sep = ''),
                                                      err.file = paste(normalizePath(dir), '/', jname, '.bsub.err', sep = ''),
                                                      run.file = ifelse(file.exists(tmp.run), tmp.run, NA),
                                                      exit_flag = NA, term_flag = NA, started = NA, reported = NA, hours_elapsed = NA, max_mem = NA, cpu_time = NA, stringsAsFactors = F);


                                    fn = paste(dir, jname, '.bsub.out', sep = '')
                                    fn.ex = file.exists(fn);
                                    if (!any(fn.ex))
                                        break

                                    tmp = matrix(unlist(parallel::mclapply(fn[fn.ex],
                                                                           function(x)
                                                                           {
                                                                               y = readLines(x);
                                                                               y = split(y, cumsum(grepl('^Sender', y)))
                                                                               y = y[[length(y)]]  ## picks "last" dump from lsf to this out file
                                                                               return(c(c(grep('^Exited with', y, value = T), grep('^Successfully completed', y, value = T), '')[1],
                                                                                        c(grep('^TERM', y, value = T), '')[1],
                                                                                        c(gsub('Started at ', '', grep('^Started at', y, value = T)), '')[1],
                                                                                        c(gsub('Results reported at ', '', grep('^Results reported at', y, value = T)), '')[1],
                                                                                        c(gsub('[ ]+CPU time[ ]+\\:[ ]+(.*)[ ]+\\S+', '\\1', grep('^[ ]+CPU time', y, value = T)), '')[1],
                                                                                        c(gsub('[ ]+Max Memory[ ]+\\:[ ]+(.*)[ ]+\\S+', '\\1', grep('^[ ]+Max Memory', y, value = T)), '')[1],
                                                                                        c(gsub('[ ]+Max Swap[ ]+\\:[ ]+(.*)[ ]+\\S+', '\\1', grep('^[ ]+Max Swap', y, value = T)), '')[1],
                                                                                        c(gsub('[ ]+Max Processes[ ]+\\:[ ]+(.*)\\S*', '\\1', grep('^[ ]+Max Processes', y, value = T)), '')[1],
                                                                                        c(gsub('[ ]+Max Threads[ ]+\\:[ ]+(.*)\\S*', '\\1', grep('^[ ]+Max Threads', y, value = T)), '')[1]
                                                                                        ))

                                                                           }, mc.cores = mc.cores)), ncol = 9, byrow = T)
                                    colnames(tmp) = c('exit.flag', 'term.flag', 'started', 'reported', 'cpu.time', 'max.memory', 'max.swap', 'max.cpu', 'max.thr')

                                    TIME.FORMAT = '%a %b %d %H:%M:%S %Y';
                                    outs$exit_flag[fn.ex] = tmp[, 'exit.flag']
                                    outs$term_flag[fn.ex] = tmp[, 'term.flag']
                                    outs$started[fn.ex] = as.character(as.POSIXct(strptime(tmp[, 'started'], TIME.FORMAT)))
                                    outs$reported[fn.ex] = as.character(as.POSIXct(strptime(tmp[, 'reported'], TIME.FORMAT)))
                                    outs$hours_elapsed = round(as.numeric((as.POSIXct(outs$reported)-as.POSIXct(outs$started))/60), 2)
                                    outs$cpu_time[fn.ex] = as.numeric(tmp[, 'cpu.time'])
                                    outs$max_mem[fn.ex] = as.numeric(tmp[, 'max.memory'])

                                    if (detailed)
                                    {
                                        outs$max_swap[fn.ex] = tmp[, 'max.swap']
                                        outs$max_processes[fn.ex] = tmp[, 'max.processes']
                                        outs$max_threads[fn.ex] = tmp[, 'max.threads']
                                    }
                                    outs$success = grepl('Success', outs$exit_flag)
                                    rownames(outs) = dedup(outs$jname)
                                }

                                outs = as.data.table(outs)

                                if (!is.null(input.jname))
                                    outs = outs[, key := input.jname[jname]]
                                else
                                    outs = outs[, key := jname]

                                setkey(outs, 'key')
                                return(outs)
                            }


###################
#' @name chunk
#' @title  chunk
#'
#' @description
#' takes same input as seq (from, to, by, length.out) and outputs a 2 column matrix of indices
#' corresponding to "chunks"
#'
#' @param from integer where to begin sequence
#' @param to integer to end sequence
#' @param by interval to space sequence
#' @param length.out number of desired chunks, i.e. nrows of output matrix
#' @return 2 column matrix of indices, each row representing chunk
#' @export
#' @author Marcin Imielinski
###################
chunk = function(from, to = NULL, by = 1, length.out = NULL)
{
    if (is.null(to))
    {
        to = from;
        from = 1;
    }

    if (is.null(length.out))
        tmp = c(seq(from = from, to = to, by = by), to + 1)
    else
        tmp = c(seq(from = from, to = to, length.out = length.out), to + 1)

    out = floor(cbind(tmp[-length(tmp)], tmp[-1]-1))

    return(out)
}


##################
#' @name func_code
#' @title func_code
#'
#' @description
#' Produces (simple) R code calling function named "func" with args in list "argv", prepending with
#' source() call to directories in the vector "sources" if specified.
#'
#' NOTE: args in ... can be lists or vectors consisting of numerical values or characters.  Lists can have named fields.
#' These will be assigned in a "hard coded" way in the Rcode, so these should be ideally scalars or
#' pretty short vectors / lists.
#'
#' For code to run properly, the names of "argv" must correspond to argument names of "func", or
#' if the list has unnamed fields then they must be ordered in the order of the function args.
#'
#' Useful for dumping tmp code files for farming when there are many arguments being passed
#'
#' @param func function to call
#' @param sources files to source
#' @param ... additional arguments to function which should contain numerical or character vectors or lists of such vectors
#' @return character data of source file containing call to function with arguments hard coded
#' @export
#' @author Marcin Imielinski
##################
func_code = function(func, sources = c(), ...)
{
    out = "";

    argv = list(...);

    if (length(sources)>0)
    {
        out = sprintf('%s%s\n', out, paste("source(\"", sources, "\")", sep = "", collapse = "\n"));
    }

    argv_strings = vector(mode="character", length = length(argv));

    for (i in 1:length(argv))
    {
        this_arg = eval(argv[[i]]); # need to eval if data frame slice passed down as vector (i.e. as "call")

        if (is.list(this_arg) & is.null(dim(this_arg))) # checks we have a bona fide list ie not a data frame
        {
            if (max(unlist(lapply(this_arg, length)))>1)
            {
                print("Error: nested list arguments not allowed in argv");
                return( NA );
            }

            list_strings = as.vector(this_arg)

            chars = unlist(lapply(this_arg, is.character));  # put quotes around char list items
            list_strings[chars] = paste('\"', list_strings[chars], '\"', sep = "");

            if (!is.null(names(this_arg))) # take care of named list items if exists
            {
                named = names(this_arg) != "";
                list_strings[named] = paste(names(this_arg)[named], " = ", list_strings[named],  sep = "");  # prepend "name=" to named items of list
            }

            argv_strings[[i]] = sprintf("list(%s)", paste(list_strings, collapse = ", "));  # pre-pend list constructor and comma concat list items
        }
        else if (is.vector(this_arg) & is.null(dim(this_arg))) # make sure we have vector and not an array
        {
            vec_strings = this_arg;
            if (is.character(this_arg))
                vec_strings = paste('\"', vec_strings, '\"', sep = "");

            if (length(vec_strings)>1) # use c() if we have a vector greater than length 1
                argv_strings[i] = sprintf("c(%s)", paste(vec_strings, collapse = ","))
            else
                argv_strings[i] = vec_strings;
        }
        else if (is.null(this_arg))
            argv_strings[i] = 'NULL'
        else
        {
            print("Error: unsupported data type in argv");
            return( NA );
        }
    }

    if (!is.null(names(argv))) # take care of named args if exist
    {
        named = names(argv) != "";
        argv_strings[named] = paste(names(argv)[named], " = ", argv_strings[named], sep = "");
    }

    out = sprintf('%s\n%s(%s)\n', out, func, paste(argv_strings, collapse = ",\n "));
    out
}


#############################
#' @name read.delim.cat
#' @title read.delim.cat
#'
#' @description
#'
#' takes a vector of tab delimited file paths and concatenates them into a
#' single data frame (takin union of identically named / numbered columns as a default)
#'
#' @param paths length n character vector of paths to tsv files
#' @param skip optional length n or length 1 integer specifying  how many lines to skip
#' @param cols optional character vector of which cols to keep (by default union of all columns)
#' @param include.paths optional logical flag whether to include paths to files as column $source.path column
#' @param include.index optional logical flag whether to include source rownames if exist as $source.id column
#' @param cores optional integer specifying number of cores to use (def 1)
#' @param ... additional args to read.delim
#' @export
#' @author Marcin Imielinski
############################
read.delim.cat = function(paths, skip = NULL, cols = NULL, include.paths = T, include.index = TRUE, cores = NULL, ...)
{
    if (is.null(skip))
        skip = rep(0, length(paths))

    paths[is.na(paths)] = "";
    does.not.exist = !file.exists(paths);

    if (any(does.not.exist))
        warning(sprintf('Ignoring %s paths that do not exist on the file system.', length(which(does.not.exist))))
    paths = paths[!does.not.exist];

    if (is.null(names(paths)))
    {
        names(paths) = 1:length(paths)
    }

    if (length(skip) ==1)
        skip = rep(skip, length(paths))
    else
        skip = skip[!does.not.exist]

                                        # scope out files to filter out those with 0 rows and find common columns
    if (!is.null(cores))
        dfs = parallel::mclapply(1:length(paths),
                                 function(x) {tmp.df = read.delim(paths[x], skip = skip[x], ...);
                                     if (nrow(tmp.df) != 0) cbind(data.frame(source.path = paths[x]), source.id = names(paths)[x], tmp.df) else data.frame() }, mc.cores = cores)
    else
        dfs = lapply(1:length(paths),
                     function(x) {tmp.df = read.delim(paths[x], skip = skip[x], ...);
                         if (nrow(tmp.df) != 0) cbind(data.frame(source.path = paths[x]), source.id = names(paths)[x], tmp.df) else data.frame() })


    dfs = dfs[sapply(dfs, nrow)!=0];

    if (length(dfs)==0)
        return(NULL);

    out = do.call('rrbind', dfs)

    if (!is.null(cols))
        out = cbind(out[,1], out[, cols]);

    if (include.paths)
        names(out)[1] = 'source.path'
    else
        out$source.path = NULL

    if (include.index)
        names(out)[1] = 'source.path'
    else
        out$source.id = NULL

    return(out)
}

#' @name fisher.plot
#' @title Plots fisher contingency table with p value
#'
#' @description
#' Plots fisher contingency table with p value
#'
#' @param O observed matrix of counts
#' @export
fisher.plot = function(O)
{
    fish = stats::fisher.test(O)
    plot.new();
    par(usr = c(0, 1, 0, 1));
    plotrix::addtable2plot(0,0.5, O, display.colnames = TRUE, display.rownames = TRUE);
    text(0.5, 0.44, sprintf(paste('P = %0.', floor(-log10(fish$p.value))+2, 'f\nOR = %0.2f [%0.2f-%0.2f]', sep  = ""), fish$p.value, fish$estimate, fish$conf.int[[1]], fish$conf.int[[2]]))
}

###############################
#' @name fisher.combined
#' @title fisher.combined
#' @description
#'
#' Computes fisher combined p value for a matrix of p values where the columns correspond to individual (independent) tests
#' rows correspond to hypotheses.
#'
#' @author Marcin Imielinski
#' @export
#' @param Ps n x k matrix of p values from k different independent tests
#' @return length n numeric vector of p values
###############################
fisher.combined = function(Ps)
{
    if (is.vector(Ps))
        return(Ps)

    return(pchisq(rowSums(-2*log(Ps)), 2*ncol(Ps), lower.tail = F))
}

get.fwf.widths = function(file, skip=0)
{
    l = readLines(file,skip+1);

    w = get.field.widths(l[[length(l)]]);
}


############################
#' @name dev.all.off
#' @title dev.all.off
#'
#' @description
#' kills all plot devices
#'
#' @author Marcin Imielinski
#' @export
#############################
dev.all.off = function()
{
    sapply(dev.list(), dev.off)
}


get.field.widths = function(str)
{
    spl = strsplit(str, "")
    non.space = is.na(match(spl[[1]], " "));

    runs = run.lengths(non.space);

    out = c();

    if (dim(runs)[1]==1)
    {
        out = length(spl[[1]]);
    }
    else if (dim(runs)[1]>1)
    {
        out = c(runs$start[2:dim(runs)[1]], length(spl[[1]])+1) -
            c(1, runs$start[2:dim(runs)[1]]);
    }

    out
}


#################
#' @name write.tab
#' @title write.tab
#' writes tab delimited no quotes without row names table (passes remaining arguments to write.table)
#'
#' equivalent to write.table(sep = TAB.DELIM, quote = F, row.names = F)
#' @param x data.frame to dump
#' @param ... additional arguments to write.table
#' @export
#################
write.tab = function(x, ..., sep = "\t", quote = F, row.names = F)
{
    if (!is.data.frame(x))
        x = as.data.frame(x)

    write.table(x, ..., sep = sep, quote = quote, row.names = row.names)
}

footprint = function(gr)
    cat(prettyNum(sum(as.numeric(width(reduce(gr)))), big.mark = ','), '\n')

#' Dump GRanges to GATK file
#'
#' Dumps gr object into gatk intervals in file path "file"
#' @param gr GRanges
#' @param file file
#' @param add.chr Flag to add "chr" to seqnames. Default FALSE
#' @return returns 0 if completed
#' @export
gr2gatk = function(gr, file, add.chr = F)
{
    sn = as.character(seqnames(gr));
    if (add.chr)
        sn = paste('chr', sn, sep = '');

    writeLines(paste(sn, ':', start(gr), '-', end(gr), sep = ''), con = file)
    return(0)
}

#' gstring
#'
#' quick function to parse gr from character vector IGV / UCSC style strings of format gr1;gr2;gr3 where each gr is of format chr:start-end[+/-]
#'
#' @name gstring
#' @export
gstring = function(...)
{
    return(unlist(parse.grl(...)))
}



#' Find peaks in a \code{GRanges} over a given meta-data field
#'
#' Finds "peaks" in an input GRanges with value field y.
#' first piles up ranges according to field score (default = 1 for each range)
#' then finds peaks.  If peel > 0, then recursively peels segments
#' contributing to top peak, and recomputes nextpeak "peel" times
#' if peel>0, bootstrap controls whether to bootstrap peak interval nbootstrap times
#' if id.field is not NULL will peel off with respect to unique (sample) id of segment and not purely according to width
#' if FUN preovided then will complex aggregating function of piled up values in dijoint intervals prior to computing "coverage"
#' (FUN must take in a single argument and return a scalar)
#' if id.field is not NULL, AGG.FUN is a second fun to aggregate values from id.field to output interval
#'
#' @param gr \code{GRanges} with some meta-data field to find peaks on
#' @param field character field specifying metadata to find peaks on, default "score, can be NULL in which case the count is computed
#' @param minima logical flag whether to find minima or maxima
#' @param id.field character denoting field whose values specifyx individual tracks (e.g. samples)
#' @param bootstrap logical flag specifying whether to bootstrap "peel off" to statistically determine peak boundaries
#' @param na.rm remove NA from data
#' @param pbootstrap  quantile of bootstrap boundaries to include in the robust peak boundary estimate (i.e. essentially specifies confidence interval)
#' @param nboostrap   number of bootstraps to run
#' @param FUN  function to apply to compute score for a single individual
#' @param AGG.FUN function to aggregate scores across individuals
#' @export
#' @examples
#'
#' ## outputs example gene rich hotspots from example_genes GRanges
#' pk = gr.peaks(example_genes)
#'
#' ## now add a numeric quantity to example_genes and compute
#' ## peaks with respect to a numeric scores, e.g. "exon_density"
#' example_genes$exon_density = example_genes$exonCount / width(example_genes)
#' pk = gr.peaks(example_genes, field = 'exon_density')
#'
#' ## can quickly find out what genes lie in the top peaks by agggregating back with
#' ## original example_genes
#' pk[1:10] %$% example_genes[, 'name']
#'
#'
gr.peaks = function(gr, field = 'score',
                    minima = FALSE,
                    peel = 0,
                    id.field = NULL,
                    bootstrap = TRUE,
                    na.rm = TRUE,
                    pbootstrap = 0.95,
                    nbootstrap = 1e4,
                    FUN = NULL,
                    AGG.FUN = sum,
                    peel.gr = NULL, ## when peeling will use these segs instead of gr (which can just be a standard granges of scores)
                    score.only = FALSE,
                    verbose = peel>0)
{

    if (!is(gr, 'GRanges'))
        gr = seg2gr(gr)

    if (is.null(field))
        field = 'score'

    if (!(field %in% names(values(gr))))
        values(gr)[, field] = 1

    if (is.logical(values(gr)[, field]))
        values(gr)[, field] = as.numeric(values(gr)[, field])

    if (peel>0 & !score.only)
    {
        if (verbose)
            cat('Peeling\n')
        out = GRanges()

        if (bootstrap)
            pbootstrap = pmax(0, pmin(1, pmax(pbootstrap, 1-pbootstrap)))

        ## peel.gr are an over-ride if we have pre-computed the score and only want to match peaks to their supporting segments
        if (is.null(peel.gr))
            peel.gr = gr

        for (p in 1:peel)
        {
            if (verbose)
                cat('Peel', p, '\n')
            if (p == 1)
                last = gr.peaks(gr, field, minima, peel = 0, FUN = FUN, AGG.FUN = AGG.FUN, id.field = id.field)
            else
            {
                ## only need to recompute peak in region containing any in.peak intervals
                in.peak = gr.in(gr, peak.hood)

                tmp = NULL
                if (any(in.peak))
                    tmp = gr.peaks(gr[in.peak, ], field, minima, peel = 0, FUN = FUN, AGG.FUN = AGG.FUN, id.field = id.field)
                last = grbind(last[!gr.in(last, peak.hood)], tmp)
                names(values(last)) = field
            }

            ## these are the regions with the maximum peak value
            mix = which(values(last)[, field] == max(values(last)[, field]))

            ## there can be more than one peaks with the same value
            ## and some are related since they are supported by the same gr
            ## we group these peaks and define a tmp.peak to span all the peaks that are related
            ## to the top peak
            ## the peak is the span beteween the first and last interval with the maximum
            ## peak value that are connected through at least one segment to the peak value

            ##
            tmp.peak = last[mix]

            if (length(tmp.peak)>1)
            {
                tmp.peak.gr = gr[gr.in(gr, tmp.peak)]
                ov = gr.findoverlaps(tmp.peak, tmp.peak.gr)
                ed = rbind(ov$query.id, ov$subject.id+length(tmp.peak))[1:(length(ov)*2)]
                cl = igraph::clusters(igraph::graph(ed), 'weak')$membership
                tmp = tmp.peak[cl[1:length(tmp.peak)] %in% cl[1]]
                peak = GRanges(seqnames(tmp)[1], IRanges(min(start(tmp)), max(end(tmp))))
                values(peak)[, field] = values(tmp.peak)[, field][1]
            }
            else
                peak = tmp.peak
            ## tmp.peak is the interval spanning all the top values in this region

            in.peak1 =  gr.in(peel.gr, gr.start(peak))
            in.peak2 = gr.in(peel.gr, gr.end(peak))
            in.peak = in.peak1 | in.peak2

            ## peak.gr are the gr supporting the peak
            peak.gr = peel.gr[in.peak1 & in.peak2] ## want to be more strict with segments used for peeling
            peak.hood = reduce(peak.gr) ## actual peak will be a subset of this, and we can this in further iterations to limit peak revision

            in.peak = rep(FALSE, length(gr))
            if (bootstrap && length(peak.gr))
            {
                ## asking across bootstrap smaples how does the intersection fluctuate
                ## among segments contributing to the peak

                if (!is.null(id.field))
                {
                    peak.gr = seg2gr(gr2dt(peak.gr)[, list(seqnames = seqnames[1], start = min(start),
                                                           eval(parse(text = paste(field, '= sum(', field, '*(end-start))/sum(end-start)'))),end = max(end)),
                                                    by = eval(id.field)])
                    names(values(peak.gr))[ncol(values(peak.gr))] = field ## not sure why I need to do this line, should be done above
                }

                B = matrix(sample(1:length(peak.gr), nbootstrap * length(peak.gr), prob = abs(values(peak.gr)[, field]), replace = TRUE), ncol = length(peak.gr))
                ## bootstrap segment samples
                ## the intersection is tha max start and min end among the segments in each
                st = apply(matrix(start(peak.gr)[B], ncol = length(peak.gr)), 1, max)
                en = apply(matrix(end(peak.gr)[B], ncol = length(peak.gr)), 1, min)

                ## take the left tail of the start position as the left peak boundary
                start(peak) = quantile(st, (1-pbootstrap)/2)

                ## and the right tail of the end position as the right peak boundary
                end(peak) = quantile(en, pbootstrap + (1-pbootstrap)/2)

                in.peak =  gr.in(gr, peak)
            }
            gr = gr[!in.peak]
            peak$peeled = TRUE
            out = c(out, peak)
            if (length(gr)==0)
                return(out)
        }
        last$peeled = FALSE
        return(c(out, last[-mix]))
    }

    if (na.rm)
        if (any(na <- is.na(values(gr)[, field])))
            gr = gr[!na]

    if (!is.null(FUN))
    {
        agr = GenomicRanges::disjoin(gr)
        values(agr)[, field] = NA
        tmp.mat = cbind(as.matrix(values(gr.val(agr[, c()], gr, field, weighted = FALSE, verbose = verbose, by = id.field, FUN = FUN, default.val = 0))))
        values(agr)[, field] = apply(tmp.mat, 1, AGG.FUN)
        gr = agr
    }

    cov = as(GenomicRanges::coverage(gr, weight = values(gr)[, field]), 'GRanges')

    if (score.only)
        return(cov)

    dcov = diff(cov$score)
    dchrom = diff(as.integer(seqnames(cov)))

    if (minima)
        peak.ix = (c(0, dcov) < 0 & c(0, dchrom)==0) & (c(dcov, 0) > 0 & c(dchrom, 0)==0)
    else
        peak.ix = (c(0, dcov) > 0 & c(0, dchrom)==0) & (c(dcov, 0) < 0 & c(dchrom, 0)==0)

    out = cov[which(peak.ix)]

    if (minima)
        out = out[order(out$score)]
    else
        out = out[order(-out$score)]

    names(values(out))[1] = field

    return(out)
}

############################################
#' ra_breaks
#'
#' takes in either file or data frame from dranger or snowman or path to BND / SV type vcf file
#' and returns junctions in VCF format.
#'
#' The default output is GRangesList each with a length two GRanges whose strands point AWAY from the break.  If get.loose = TRUE (only relevant for VCF)
#'
#' @name ra_breaks
#' @import VariantAnnotation
#' @export
############################################
ra_breaks = function(rafile,
                     keep.features = T,
                     seqlengths = hg_seqlengths(),
                     chr.convert = T,
                     geno=NULL,
                     flipstrand = FALSE,
                     swap.header = NULL,
                     breakpointer = FALSE,
                     seqlevels = NULL,
                     force.bnd = FALSE,
                     skip = NA,
                     get.loose = FALSE){
    ## if TRUE will return a list with fields $junctions and $loose.ends
    if (is.character(rafile))
    {
        if (grepl('.rds$', rafile)){
            ra = readRDS(rafile)
            ## validity check written for "junctions" class
            return(junctions(ra))
        }
        else if (grepl('(.bedpe$)', rafile)){
            ra.path = rafile
            cols = c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'name', 'score', 'str1', 'str2')

            ln = readLines(ra.path)
            if (is.na(skip))
            {
                nh = min(c(Inf, which(!grepl('^((#)|(chrom))', ln))))-1
                if (is.infinite(nh)){
                    nh = 1
                }
            }
            else{
                nh = skip
            }

            if ((length(ln)-nh)==0){
                ## if (get.loose){
                ##     return(list(junctions = GRangesList(GRanges(seqlengths = seqlengths))[c()], loose.ends = GRanges(seqlengths = seqlengths)))
                ## }
                ## else{
                return(GRangesList(GRanges(seqlengths = seqlengths))[c()])
                ## }
            }

            if (nh ==0){
                rafile = fread(rafile, header = FALSE)
            }
            else
            {

                rafile = tryCatch(fread(ra.path, header = FALSE, skip = nh), error = function(e) NULL)
                if (is.null(rafile)){
                    rafile = tryCatch(fread(ra.path, header = FALSE, skip = nh, sep = '\t'), error = function(e) NULL)
                }

                if (is.null(rafile)){
                    rafile = tryCatch(fread(ra.path, header = FALSE, skip = nh, sep = ','), error = function(e) NULL)
                }

                if (is.null(rafile)){
                    stop('Error reading bedpe')
                }
            }
            setnames(rafile, 1:length(cols), cols)
            rafile[, str1 := ifelse(str1 %in% c('+', '-'), str1, '*')]
            rafile[, str2 := ifelse(str2 %in% c('+', '-'), str2, '*')]
        }
        else if (grepl('(vcf$)|(vcf.gz$)', rafile)){
            
            require(VariantAnnotation)
            vcf = readVcf(rafile, Seqinfo(seqnames = names(seqlengths), seqlengths = seqlengths))
            ## vgr = rowData(vcf) ## parse BND format
            vgr = skidb::read_vcf(rafile, swap.header = swap.header, geno=geno)
            if (!is.null(info(vcf)$SCTG))
                vgr$SCTG = info(vcf)$SCTG
            
            return(vgr2ra(vgr, force.bnd = force.bnd, get.loose = get.loose))
        }
        else{
            rafile = read.delim(rafile)
        }
    }

    if (is.data.table(rafile))
    {
        require(data.table)
        rafile = as.data.frame(rafile)
    }

    if (nrow(rafile)==0)
    {
        out = GRangesList()
        values(out) = rafile
        return(out)
    }

    if (flipstrand) ## flip breaks so that they are pointing away from junction
    {
        rafile$str1 = ifelse(rafile$strand1 == '+', '-', '+')
        rafile$str2 = ifelse(rafile$strand2 == '+', '-', '+')
    }

    if (!is.null(seqlevels)) ## convert seqlevels from notation in tab delim file to actual
    {
        rafile$chr1 = seqlevels[rafile$chr1]
        rafile$chr2 = seqlevels[rafile$chr2]
    }


    if (is.null(rafile$str1)){
        rafile$str1 = rafile$strand1
    }

    if (is.null(rafile$str2)){
        rafile$str2 = rafile$strand2
    }

    if (!is.null(rafile$pos1) & !is.null(rafile$pos2))
    {
        if (breakpointer)
        {
            rafile$pos1 = rafile$T_BPpos1
            rafile$pos2 = rafile$T_BPpos2
        }

        if (!is.numeric(rafile$pos1)){
            rafile$pos1 = as.numeric(rafile$pos1)
        }

        if (!is.numeric(rafile$pos2)){
            rafile$pos2 = as.numeric(rafile$pos2)
        }

        ## clean the parenthesis from the string

        rafile$str1 <- gsub('[()]', '', rafile$str1)
        rafile$str2 <- gsub('[()]', '', rafile$str2)

        ## goal is to make the ends point <away> from the junction where - is left and + is right
        if (is.character(rafile$str1) | is.factor(rafile$str1)){
            rafile$str1 = gsub('0', '-', gsub('1', '+', gsub('\\-', '1', gsub('\\+', '0', rafile$str1))))
        }

        if (is.character(rafile$str2) | is.factor(rafile$str2)){
            rafile$str2 = gsub('0', '-', gsub('1', '+', gsub('\\-', '1', gsub('\\+', '0', rafile$str2))))
        }


        if (is.numeric(rafile$str1)){
            rafile$str1 = ifelse(rafile$str1>0, '+', '-')
        }

        if (is.numeric(rafile$str2)){
            rafile$str2 = ifelse(rafile$str2>0, '+', '-')
        }

        rafile$rowid = 1:nrow(rafile)

        bad.ix = is.na(rafile$chr1) | is.na(rafile$chr2) | is.na(rafile$pos1) | is.na(rafile$pos2) | is.na(rafile$str1) | is.na(rafile$str2) | rafile$str1 == '*'| rafile$str2 == '*' | rafile$pos1<0 | rafile$pos2<0

        rafile = rafile[which(!bad.ix), ]

        if (nrow(rafile)==0){
            return(GRanges())
        }

        seg = rbind(data.frame(chr = rafile$chr1, pos1 = rafile$pos1, pos2 = rafile$pos1, strand = rafile$str1, ra.index = rafile$rowid, ra.which = 1, stringsAsFactors = F),
                    data.frame(chr = rafile$chr2, pos1 = rafile$pos2, pos2 = rafile$pos2, strand = rafile$str2, ra.index = rafile$rowid, ra.which = 2, stringsAsFactors = F))

        if (chr.convert){
            seg$chr = gsub('chr', '', gsub('25', 'M', gsub('24', 'Y', gsub('23', 'X', seg$chr))))
        }

        out = seg2gr(seg, seqlengths = seqlengths)[, c('ra.index', 'ra.which')];
        out = split(out, out$ra.index)
    }
    else if (!is.null(rafile$start1) & !is.null(rafile$start2) & !is.null(rafile$end1) & !is.null(rafile$end2))
    {
        ra1 = gr.flipstrand(GRanges(rafile$chr1, IRanges(rafile$start1, rafile$end1), strand = rafile$str1))
        ra2 = gr.flipstrand(GRanges(rafile$chr2, IRanges(rafile$start2, rafile$end2), strand = rafile$str2))
        out = grl.pivot(GRangesList(ra1, ra2))
    }



    if (keep.features){
        values(out) = rafile[, ]
    }

    if (!is.null(pad)){
        out = ra.dedup(out, pad = pad)
    }

    if (!get.loose){
        return(out)
    }
    else{
        return(list(junctions = out, loose.ends = GRanges()))
    }

    return(new("junctions", out))
}


vgr2ra = function(vgr, force.bnd = FALSE, get.loose = FALSE)
{    
    mc = data.table(as.data.frame(mcols(vgr)))
    
    if (!('SVTYPE' %in% colnames(mc))) {
        warning('Vcf not in proper format.  Is this a rearrangement vcf?')
        return(GRangesList());
    }
    
    if (any(w.0 <- (width(vgr)<1))){
        warning("Some breakpoint width==0.")
        ## right bound smaller coor
        ## and there's no negative width GR allowed
        vgr[which(w.0)] = gr.start(vgr[which(w.0)]) %-% 1
    }
    
    ## BND format doesn't have duplicated rownames
    if (any(duplicated(names(vgr)))) names(vgr) = NULL
    
    ## no events
    if (length(vgr) == 0)
        return (GRangesList())
    
    ## local function that turns old VCF to BND
    .vcf2bnd = function(vgr){
        if (!"END" %in% colnames(values(vgr)))
            stop("Non BND SV should have the second breakpoint coor in END columns!")
        
        if (!"CHR2" %in% colnames(values(vgr)) | any(is.na(vgr$CHR2)))
            vgr$CHR2 = as.character(seqnames(vgr))
        
        bp2 = data.table(as.data.frame(mcols(vgr)))
        bp2[, ":="(seqnames=CHR2, start=as.numeric(END), end=as.numeric(END))]
        bp2.gr = dt2gr(bp2)
        mcols(bp2.gr) = mcols(vgr)

        if (!is.null(names(vgr)) & !anyDuplicated(names(vgr))){
            jid = names(vgr)
        } else {
            jid = seq_along(vgr)
        }

        if (length(vgr)==0)
            return(vgr)

        names(vgr) = paste(paste0("exp", jid), "1", sep=":")
        names(bp2.gr) = paste(paste0("exp", jid), "2", sep=":")

        vgr=resize(c(vgr, bp2.gr), 1)

        if (all(grepl("[_:][12]$",names(vgr)))){
            ## row naming same with Snowman
            nm <- vgr$MATEID <- names(vgr)
            ix <- grepl("1$",nm)
            vgr$MATEID[ix] = gsub("(.*?)(1)$", "\\12", nm[ix])
            vgr$MATEID[!ix] = gsub("(.*?)(2)$", "\\11", nm[!ix])
            vgr$SVTYPE="BND"
        }
        return(vgr)
    }

    ## TODO: Delly and Novobreak
    ## fix mateids if not included
    if (!"MATEID" %in% colnames(mcols(vgr))) {
        ## TODO: don't assume every row is a different junction
        ## Novobreak, I'm looking at you.
        ## now delly...
        ## if SVTYPE is BND but no MATEID, don't pretend to be
        if (length(fake.bix <- which(values(vgr)$SVTYPE=="BND"))!=0){
            values(vgr[fake.bix])$SVTYPE = "TRA"
        }

        ## add row names just like Snowman
        if (all(names(vgr)=="N" | ## Novobreak
                is.null(names(vgr)) |
                all(grepl("^DEL|DUP|INV|BND", names(vgr)))) ## Delly
            ){
            ## otherwise if all "N", as Novobreak
            ## or starts with DEL|DUP|INV|BND, as Delly
            ## expand and match MATEID
            vgr=.vcf2bnd(vgr)
        }
    } else if (any(is.na(mid <- as.character(vgr$MATEID)))){
        ## like Lumpy, the BND rows are real BND but blended with non-BND rows
        ## treat them separately
        if (is.null(vgr$CHR2)){
            vgr$CHR2 = as.character(NA)
        }

        names(vgr) = gsub("_", ":", names(vgr))
        vgr$MATEID = sapply(vgr$MATEID, function(x) gsub("_", ":", x))

        values(vgr) = data.table(as.data.frame(values(vgr)))

        ## break up the two junctions in one INV line!
        if ("STRANDS" %in% colnames(mc) & any(ns <- sapply(vgr$STRANDS, length)>1)){
            ## first fix format errors, two strand given, but not comma separeted
            ## so you'd have taken them as single
            if (any(fuix <- sapply(vgr[which(!ns)]$STRANDS, str_count, ":")>1)){
                which(!ns)[fuix] -> tofix
                vgr$STRANDS[tofix] = lapply(vgr$STRANDS[tofix],
                                            function(x){
                                                strsplit(gsub("(\\d)([\\+\\-])", "\\1,\\2", x), ",")[[1]]
                                            })
                ns[tofix] = TRUE
            }

            ## for the one line two junction cases
            ## split into two lines
            vgr.double = vgr[which(ns)]
            j1 = j2 = vgr.double
            st1 = lapply(vgr.double$STRANDS, function(x)x[1])
            st2 = lapply(vgr.double$STRANDS, function(x)x[2])
            j1$STRANDS = st1
            j2$STRANDS = st2
            vgr.double = c(j1, j2)
            names(vgr.double) = dedup(names(vgr.double))
            vgr = c(vgr[which(!ns)], vgr.double)
        }

        mid <- as.logical(sapply(vgr$MATEID, length))
        vgr.bnd = vgr[which(mid)]
        vgr.nonbnd = vgr[which(!mid)]

        vgr.nonbnd = .vcf2bnd(vgr.nonbnd)

        mc.bnd = data.table(as.data.frame(values(vgr.bnd)))
        mc.nonbnd = data.table(as.data.frame(values(vgr.nonbnd)))
        mc.bnd$MATEID = as.character(mc.bnd$MATEID)

        vgr = c(vgr.bnd[,c()], vgr.nonbnd[,c()])
        values(vgr) = rbind(mc.bnd, mc.nonbnd)
    }

    ## sanity check
    if (!any(c("MATEID", "SVTYPE") %in% colnames(mcols(vgr))))
        stop("MATEID or SVTYPE not included. Required")

    vgr$mateid = vgr$MATEID
    ## what's this???
    vgr$svtype = vgr$SVTYPE

    if (force.bnd)
        vgr$svtype = "BND"

    if (any(is.na(names(vgr))))
    {
        stop('vgr names not provided, input likely malformed')
    }

    if (any(is.na(vgr$svtype)))
    {
        warning('rearrangements found with NA SVTYPE will assume BND for these')
        vgr$svtype[is.na(vgr$svtype)] = 'BND'
    }

    if (sum(vgr$svtype == 'BND')==0)
        warning('Vcf not in proper format.  Will treat rearrangements as if in BND format')

    if (!all(vgr$svtype == 'BND')){
        warning(sprintf('%s rows of vcf do not have svtype BND, treat them as non-BND!',
                        sum(vgr$svtype != 'BND')))

    }

    bix = which(vgr$svtype == "BND")
    vgr = vgr[bix]
    alt <- sapply(vgr$ALT, function(x) x[1])

    ## Determine each junction's orientation
    if ("CT" %in% colnames(mcols(vgr))){
        message("CT INFO field found.")
        if ("SVLEN" %in% colnames(values(vgr))){
            ## proceed as Novobreak
            ## ALERT: overwrite its orientation!!!!
            del.ix = which(vgr$SVTYPE=="DEL")
            dup.ix = which(vgr$SVTYPE=="DUP")
            vgr$CT[del.ix] = "3to5"
            vgr$CT[dup.ix] = "5to3"
        }

        ## also, Delly is like this
        ori = strsplit(vgr$CT, "to")
        iid = sapply(strsplit(names(vgr), ":"), function(x)as.numeric(x[2]))
        orimap = setNames(c("+", "-"), c("5", "3"))
        strd = orimap[sapply(seq_along(ori), function(i) ori[[i]][iid[i]])]
        strand(vgr) = strd
        vgr.pair1 = vgr[which(iid==1)]
        vgr.pair2 = vgr[which(iid==2)]
    }
    else if ("STRANDS" %in% colnames(mcols(vgr))){
        ## TODO!!!!!!!!!!!!!!!
        ## sort by name, record bp1 or bp2
        message("STRANDS INFO field found.")
        iid = sapply(strsplit(names(vgr), ":"), function(x)as.numeric(x[2]))
        vgr$iid = iid
        vgr = vgr[order(names(vgr))]
        iid = vgr$iid

        ## get orientations
        ori = strsplit(substr(unlist(vgr$STRANDS), 1, 2), character(0))
        orimap = setNames(c("+", "-"), c("-", "+"))

        ## map strands
        strd = orimap[sapply(seq_along(ori), function(i) ori[[i]][iid[i]])]
        strand(vgr) = strd

        vgr.pair1 = vgr[which(iid==1)]
        vgr.pair2 = vgr[which(iid==2)]
    }
    else if (any(grepl("[\\|[\\]]", alt, perl = TRUE))){
        message("ALT field format like BND")
        ## proceed as Snowman
        vgr$first = !grepl('^(\\]|\\[)', alt) ## ? is this row the "first breakend" in the ALT string (i.e. does the ALT string not begin with a bracket)
        vgr$right = grepl('\\[', alt) ## ? are the (sharp ends) of the brackets facing right or left
        vgr$coord = as.character(paste(seqnames(vgr), ':', start(vgr), sep = ''))
        vgr$mcoord = as.character(gsub('.*(\\[|\\])(.*\\:.*)(\\[|\\]).*', '\\2', alt))
        vgr$mcoord = gsub('chr', '', vgr$mcoord)

        ## add extra genotype fields to vgr
        if (all(is.na(vgr$mateid)))
            if (!is.null(names(vgr)) & !any(duplicated(names(vgr))))
            {
                warning('MATEID tag missing, guessing BND partner by parsing names of vgr')
                vgr$mateid = paste(gsub('::\\d$', '', names(vgr)),
                (sapply(strsplit(names(vgr), '\\:\\:'), function(x) as.numeric(x[length(x)])))%%2 + 1, sep = '::')
            }
            else if (!is.null(vgr$SCTG))
            {
                warning('MATEID tag missing, guessing BND partner from coordinates and SCTG')
                require(igraph)
                ucoord = unique(c(vgr$coord, vgr$mcoord))
                vgr$mateid = paste(vgr$SCTG, vgr$mcoord, sep = '_')

                if (any(duplicated(vgr$mateid)))
                {
                    warning('DOUBLE WARNING! inferred mateids not unique, check VCF')
                    bix = bix[!duplicated(vgr$mateid)]
                    vgr = vgr[!duplicated(vgr$mateid)]
                }
            }
            else{
                stop('Error: MATEID tag missing')
            }

        vgr$mix = as.numeric(match(vgr$mateid, names(vgr)))

        pix = which(!is.na(vgr$mix))

        vgr.pair = vgr[pix]

        if (length(vgr.pair)==0){
            stop('Error: No mates found despite nonzero number of BND rows in VCF')
        }

        vgr.pair$mix = match(vgr.pair$mix, pix)

        vix = which(1:length(vgr.pair)<vgr.pair$mix)
        vgr.pair1 = vgr.pair[vix]
        vgr.pair2 = vgr.pair[vgr.pair1$mix]

        ## now need to reorient pairs so that the breakend strands are pointing away from the breakpoint

        ## if "first" and "right" then we set this entry "-" and the second entry "+"
        tmpix = vgr.pair1$first & vgr.pair1$right
        if (any(tmpix))
        {
            strand(vgr.pair1)[tmpix] = '-'
            strand(vgr.pair2)[tmpix] = '+'
        }

        ## if "first" and "left" then "-", "-"
        tmpix = vgr.pair1$first & !vgr.pair1$right
        if (any(tmpix))
        {
            strand(vgr.pair1)[tmpix] = '-'
            strand(vgr.pair2)[tmpix] = '-'
        }

        ## if "second" and "left" then "+", "-"
        tmpix = !vgr.pair1$first & !vgr.pair1$right
        if (any(tmpix))
        {
            strand(vgr.pair1)[tmpix] = '+'
            strand(vgr.pair2)[tmpix] = '-'
        }

        ## if "second" and "right" then "+", "+"
        tmpix = !vgr.pair1$first & vgr.pair1$right
        if (any(tmpix))
        {
            strand(vgr.pair1)[tmpix] = '+'
            strand(vgr.pair2)[tmpix] = '+'
        }

        pos1 = as.logical(strand(vgr.pair1)=='+') ## positive strand junctions shift left by one (i.e. so that they refer to the base preceding the break for these junctions
        if (any(pos1))
        {
            start(vgr.pair1)[pos1] = start(vgr.pair1)[pos1]-1
            end(vgr.pair1)[pos1] = end(vgr.pair1)[pos1]-1
        }

        pos2 = as.logical(strand(vgr.pair2)=='+') ## positive strand junctions shift left by one (i.e. so that they refer to the base preceding the break for these junctions
        if (any(pos2))
        {
            start(vgr.pair2)[pos2] = start(vgr.pair2)[pos2]-1
            end(vgr.pair2)[pos2] = end(vgr.pair2)[pos2]-1
        }
    }

    ra = grl.pivot(GRangesList(vgr.pair1[, c()], vgr.pair2[, c()]))

    ## ALERT: vgr has already been subsetted to only include BND rows
    ## bix is the original indices, so NOT compatible!
    ## this.inf = values(vgr)[bix[pix[vix]], ]
    if (exists("pix") & exists("vix")) this.inf = values(vgr)[pix[vix], ]
    if (exists("iid")) this.inf = values(vgr[which(iid==1)])

    if (is.null(this.inf$POS)){
        this.inf = cbind(data.frame(POS = ''), this.inf)
    }
    if (is.null(this.inf$CHROM)){
        this.inf = cbind(data.frame(CHROM = ''), this.inf)
    }

    if (is.null(this.inf$MATL)){
        this.inf = cbind(data.frame(MALT = ''), this.inf)
    }

    this.inf$CHROM = seqnames(vgr.pair1)
    this.inf$POS = start(vgr.pair1)
    this.inf$MATECHROM = seqnames(vgr.pair2)
    this.inf$MATEPOS = start(vgr.pair2)
    this.inf$MALT = vgr.pair2$AL

    ## NOT SURE WHY BROKEN
    ## tmp = tryCatch(cbind(values(vgr)[bix[pix[vix]],], this.inf), error = function(e) NULL)
    ## if (!is.null(tmp))
    ##     values(ra) = tmp
    ## else
    ##     values(ra) = cbind(vcf@fixed[bix[pix[vix]],], this.inf)

    values(ra) = this.inf

    if (is.null(values(ra)$TIER)){
        ## baseline tiering of PASS vs non PASS variants
        ## ALERT: mind the naming convention by diff programs
        ## TODO: make sure it is compatible with Delly, Novobreak, Meerkat
        ## Snowman/SvABA uses "PASS"
        ## Lumpy/Speedseq uses "."
        values(ra)$tier = ifelse(values(ra)$FILTER %in% c(".", "PASS"), 2, 3)
    }
    else {
        values(ra)$tier = values(ra)$TIER
    }

    ra = ra.dedup(ra)

    if (!get.loose | is.null(vgr$mix)){
        return(ra)
    }
    else
    {
        npix = is.na(vgr$mix)
        vgr.loose = vgr[npix, c()] ## these are possible "loose ends" that we will add to the segmentation

        ## NOT SURE WHY BROKEN
        tmp =  tryCatch( values(vgr)[bix[npix], ],
                        error = function(e) NULL)
        if (!is.null(tmp)){
            values(vgr.loose) = tmp
        }
        else{
            values(vgr.loose) = cbind(vcf@fixed[bix[npix], ], info(vcf)[bix[npix], ])
        }

        return(list(junctions = ra, loose.ends = vgr.loose))
    }
}

#' @name score.walks
#' @title score.walks
#' @description
#'
#' Scores GRangesList of walks against GRanges of 10X reads with $BX tag.
#'
#' @param wks GRangesList of walks
#' @param bam bam file
#' @param reads GRanges of reads with BX tag
#' @param win  genomic window in which to score (default is just reduce(unlist(wks))))
#' @param wins tiles to chop up genome further (beyond walk segments)
#' @param raw  returns raw barcode by walk matrix of barcode scores
#'
#' @import gChain
#' @return scores of walks or (if raw == tRUE) raw barcode to walk maps
#' @export
#' @author Marcin Imielinski
score.walks = function(wks, bam = NULL, reads = NULL, win = NULL, wins = NULL, rthresh= 4, thresh = 1e5, pad = 1e4, raw = FALSE, allpaths = TRUE, verbose = TRUE)
{
    shift = data.table::shift
    rowSums = Matrix::rowSums
    colSums = Matrix::colSums
    if (is.null(wins))
    {

        tmp = unique(disjoin(gr.stripstrand(unlist(wks))))
        wins = sort(disjoin(c(gr.start(tmp, pad), gr.end(tmp, pad))))
        strand(wins) = '+'
    }

    ## add 1 unit of "padding" to any cyclic walks to adequately measure
    cyc.ix = values(wks)$is.cyc

    if (any(cyc.ix))
        wks[cyc.ix] = do.call(GRangesList, lapply(which(cyc.ix), function(x) c(wks[[x]], wks[[x]])))


    THRESH = thresh

    if (!is.null(win))
        wins = wins[gr.in(wins, win)]

    if (verbose)
        message('Total territory to analyze is ', round(sum(as.numeric(width(wins)))/1e6,2), 'MB')

    if (sum(as.numeric(width(wins)))/1e6==0)
        stop('No walk areas intersect with provided win')

    reads.dt = NULL;
    if (!is.null(bam))
    {
        if (verbose)
            message('Pulling out reads')

        reads = dt2gr(read.bam(bam, streduce(wins), tag = 'BX', as.data.table = TRUE))
    }

    if (!inherits(reads, 'GRanges'))
        reads = dt2gr(reads)


    if (verbose)
        message("Computing insert size distro for ", length(reads), " reads")

    reads = reads[!is.na(reads$BX), ]
    reads.dt = as.data.table(reads)

    ## rthresh is reads per barcode filter
    ## i.e. remove all barcodes with fewer than rthresh reads per barcode
    if (!is.na(rthresh))    {
        keep.bx = reads.dt[, length(start), keyby = "BX"][V1>=rthresh, BX]
        reads.dt = reads.dt[BX %in% keep.bx, ]
    }
    bxlev = unique(reads$BX)

    zthresh = 3
    
    reads.dt[, sn:= as.integer(seqnames)]
    reads.dt[, str := strand == '+']

    ## nullify isize for discordant pairs
    if (verbose)
        message("Identifying discordant pairs")

    reads.dt[, R1 := bamflag(flag)[, "isFirstMateRead"]==1]
    reads.dt[, both := any(R1) & any(!R1), by = qname]
    reads.dt[both == TRUE, ":="(sn.diff = any(diff(sn)!=0),
                                first.strand.pos = str[1],
                                other.strand.pos = any(str[R1!=R1[1]])
                                ), by = qname]

    reads.dt[first.strand.pos & !other.strand.pos & !sn.diff, insert.size := max(end[R1!=R1[1]])-start[1], by = qname]

                                        #reads.dt[, insert.sizez := scale(insert.size)]
    ithresh.high = quantile(reads.dt$insert.size, 0.99, na.rm = TRUE)
    ithresh.low = quantile(reads.dt$insert.size, 0.80, na.rm = TRUE)
    reads.dt[, count := length(start), by = qname]
    reads.dt[both == TRUE, discordant := insert.size > ithresh.high]
    init.disc = reads.dt[!duplicated(qname), sum(discordant, na.rm = TRUE)]
    
                                        #      ithresh.high = reads.dt[insert.sizez>zthresh, min(insert.size)]


    ## filter "short dup" read pairs from  discordants (- to the left of + read ...)
    ## which seems to be artifact mode in 10X data
    ## (i.e. no longer call them discordant)
    dthresh = 1e4
    reads.dt[discordant == TRUE, ddist := abs(end-mpos)]
    reads.dt[discordant == TRUE & ddist<dthresh & sn.diff == 0 & !first.strand.pos & other.strand.pos, discordant := NA]

    if (verbose)
  {
    final.disc = reads.dt[!duplicated(qname), sum(discordant, na.rm = TRUE)]
    message('Found ', final.disc, ' discordant pairs after removing ', init.disc - final.disc, ' small dup-like pairs')
  }

  if (verbose)
    message("Identifying barcode strobe width")

    setkeyv(reads.dt, c("seqnames", "start"))
    reads.dt[which(!discordant & R1 == TRUE), bx.diff := c((start-shift(end))[-1], NA), by = .(seqnames, BX)]
    reads.dt[which(!discordant & R1 == FALSE), bx.diff := c((start-shift(end))[-1], NA), by = .(seqnames, BX)]
    reads.dt[, bx.diffz := scale(log(pmax(0, bx.diff)+1))]
    bzthresh = 1.5
    bthresh = reads.dt[bx.diffz>bzthresh, min(bx.diff)]
    bmean = reads.dt[, mean(bx.diff, na.rm = TRUE)]

    ## concordant and discordant read pairs
    readsc = dt2gr(reads.dt)

    ## discordant read pairs --> strand flip secnod read in pair
    readsd = dt2gr(reads.dt[which(discordant), ][R1==FALSE, strand := c('+'='-', '-'='+')[strand]])

    if (verbose)
      message("Collapsing concordant linked reads by inferred strobe width ", bthresh)
    ## collapse / reduce concordant read pairs

#### ALT approach for read cloud generation given thresh
    .reads2clouds = function(reads, thresh = bthresh)
    {
      reads = gr2dt(reads)
      setkeyv(reads, c("seqnames", "start"))
#      if (is.null(reads$bx.diff))
      reads[, bx.diff := c((start-data.table::shift(end))[-1], NA), by = .(seqnames, BX)]
      reads[, rl := label.runs(bx.diff<thresh | is.na(bx.diff)), by = .(seqnames, BX)]
      reads[, rl.last := data.table::shift(rl), by = .(seqnames, BX)]
      reads[is.na(rl), rl := ifelse(is.na(rl.last), -(1:.N), rl.last)] ## label remaining loners
      reads[, rll := paste(seqnames, BX, rl, sep = '_')]
      reads = dt2gr(reads[, .(start = start[1], end = end[.N]), by = .(seqnames, BX, rll)])
      return(reads)
    }


    readsc = .reads2clouds(readsc, thresh = bthresh)

    readsc$BX = factor(readsc$BX, bxlev)
    readsd$BX = factor(readsd$BX, bxlev)

    wov = grl.unlist(wks)[, 'grl.ix'] %*% wins[, c()]

    ## matrix of base pair width overlap between walks and wins
    wovmat = sparseMatrix(as.integer(wov$grl.ix), wov$subject.id, x = as.numeric(width(wov)), dims = c(length(wks), length(wins)))

    if (length(reads)==0)
      stop("No reads with non NA BX provided, please check input")

    ## for discordant pairs ...
    ## we want to directly assess intersection with walks
    ## in a strand specific way
    wksu = grl.unlist(wks) ## these are unlisted
    wksur = gr.flipstrand(wksu) ## these are strand flipped unlisted

    qmap = as.data.table(readsd)[, .(qname, BX)][, BX[1], keyby = qname][, structure(V1, names = as.character(qname))]
    qlev = names(qmap)

    if (verbose)
      message("Lifting discordant reads onto walks")

    ## lift read onto walk coordinates using gChain
    wk.nm = names(wks)
    names(wks) = 1:length(wks)
    wks.chain = gChain::spChain(wks)

    ## now we want to ask what are the read pairs that now become concordant
    ## on each lifted walk??
    ## then score per BX and walk, how many discordant pairs are made concordant post-lift
    ## and finally note which barcodes have the maximum number of their discordant pairs
    ## lifted onto the walk
    readsdl = gChain::lift(wks.chain, readsd)
    readsdl.dt = as.data.table(readsdl)[order(seqnames, start), ]

    ## use similar criteria to above to identify discordant / concordant reads in "lifted coordinates"
    readsdl.dt = readsdl.dt[, both := any(R1) & any(!R1), by = .(seqnames, qname)][both == TRUE, ]
    readsdl.dt[both == TRUE, ":="(
                               first.strand.pos = str[1],
                               other.strand.pos = any(str[R1!=R1[1]])
                             ), by = qname]
    readsdl.dt[first.strand.pos & !other.strand.pos, insert.size := end[R1!=R1[1]][1]-start[1], by = qname]
    readsdl.dt[, insert.size := end[R1!=R1[1]][1]-start[1], by = qname]
    readsdl.dt[, concordant := insert.size<ithresh.low]
    readsdl.dt = readsdl.dt[concordant == TRUE, ][R1 == TRUE, ][, dup := duplicated(query.id), by = .(BX, seqnames)]
    readsdl.dt = readsdl.dt[dup==FALSE, ]
    bxstatsd = readsdl.dt[ , .(score = length(qname)), by = .(BX, seqnames)]
    bxstatsd[, max.score := max(score), by = .(BX)]
    bxstatsd = bxstatsd[score == max.score, ] ## only keep the max scoring BX, seqnames pairs

    ## we want to use these as votes for walk support, but then any other non matching discordant pairs as anti-matches
    ## first building a mat rix of qnames x walks
    rovd.mat = sparseMatrix(as.integer(factor(bxstatsd$BX, bxlev)), as.numeric(as.character(bxstatsd$seqnames)), x = bxstatsd$score,
                            dims = c(length(bxlev), length(wks)), dimnames = list(bxlev, 1:length(wks)))

    if (verbose)
      message("Lifting concordant linked read footprints onto walks")

    rovcl = gChain::lift(wks.chain, readsc)
    rovclb = gChain::lift(gChain::t(wks.chain),  rovcl)
    values(rovclb)$walk = as.integer(seqnames(rovcl)[rovclb$query.id]) ## without as.integer(), gr.findoverlaps fails below
    values(rovclb)$bx.walk = paste(values(rovclb)$BX, values(rovclb)$walk)

    ##   tmp = rep(readsc, length(wks))
    ##   tmp$walk = rep(1:length(wks), each = length(readsc))
    ##   leftover.ix = setdiff(1:length(tmp), gr.findoverlaps(tmp, rovclb, by = c("BX", "walk"))$query.id)
    ##   leftovers = tmp[leftover.ix]

    ##   ## count "leftover" clouds per BX per walk
    ## leftovers

    rovclbr = grl.reduce(split(rovclb, values(rovclb)$bx.walk))
    bxwid.lift = as.data.table(matrix(unlist(strsplit(names(rovclbr), ' ')), ncol= 2, byrow = TRUE))
    setnames(bxwid.lift, c('BX', 'seqnames'))
    bxwid.lift[, wid.lifted := grl.eval(rovclbr, sum(as.numeric(width)))]

    if (verbose)
      message("Analyzing concordant walk footprints on walks")

    ## for every barcode we want to ask how much of its width
    ## is "missing" post lift to that walk?  that's going to drive a negative score
    ## with regard to its match to a given walk
    bxwid = as.data.table(readsc)[, .(wid = sum(as.numeric(width))), keyby = BX]
    rsc = split(readsc, readsc$BX)
    bxwid.max = data.table(BX = names(rsc), width.max.og = grl.eval(rsc, max(width)))
    setkey(bxwid.max, BX)

    bxwid.lift = as.data.table(rovcl)[, .(wid.lifted = sum(width)), keyby = .(BX, seqnames)] 
### wid.left == wid not lifted 
    bxwid.lift[bxwid, wid.left := wid-wid.lifted, on = 'BX']

    ## neg.mat = negative overlap matrix from lift
    neg.mat = sparseMatrix(as.integer(factor(bxwid.lift$BX, bxlev)), as.numeric(as.character(bxwid.lift$seqnames)), x = bxwid.lift$wid.left, dims = c(length(bxlev), length(wks)), dimnames = list(bxlev, 1:length(wks)))

    ## reduce the footprint of each BX on each walk + bthresh pad
    ## rovcl.fp = as.data.table(grl.reduce(split(rovcl + bthresh, rovcl$BX)))[, BX := group_name]
    ##  rovcl.fp = as.data.table(.reads2clouds(rovcl, bthresh))[width>bthresh ,]
    rovcl.mfp = as.data.table(.reads2clouds(rovcl, bthresh))[, .(width.max.lifted = max(width)), by = .(seqnames, BX)]
    
#    setkeyv(rovcl.fp, c("seqnames", "BX"))
#    rovcl.fp[, gaps := start-shift(end),  by = .(seqnames, BX)]
#    rovcl.fp[is.na(gaps), gaps := 0]
#    rovcl.fp[, pgap := dexp(gaps, 1/bmean, log = TRUE)]

    ## 
    bxstats = (bxwid.lift %>% merge(rovcl.mfp, by = c('seqnames', 'BX')) %>% merge(bxwid.max, by = 'BX'))

    ## keep only those that have better width.max.lifted than width.max.og
    left.thresh = 0.01    
    bxstats = bxstats[wid.left<left.thresh*width.max.lifted & width.max.lifted>width.max.og, ]

    ## log.sum.exp =  function(x){
    ##   offset = max(x)
    ##   log(sum(exp(x - offset))) + offset
    ## }
   
    ## ## calculate widths of (largest vs next largest) footprints per walk
    ## bxstats = rovcl.fp[, .(jpgap = sum(pgap)), keyby = .(BX, seqnames)]
    ## bxstats[, lse := log.sum.exp(jpgap), by = BX]
    ## bxstats[, pw := exp(jpgap-lse), by = BX]
     
    
    if (verbose)
      message("Creating barcode x walk matrices")

    ## rovc.mat = convert bxstats to barcode x walk matrix
    ##  rovc.mat = sparseMatrix(as.integer(factor(bxstats$BX, bxlev)), as.numeric(as.character(bxstats$seqnames)), x = bxstats$wid.rel, dims = c(length(bxlev), length(wks)), dimnames = list(bxlev, 1:length(wks)))

    rovc.mat = sparseMatrix(as.integer(factor(bxstats$BX, bxlev)), as.numeric(as.character(bxstats$seqnames)), x = 1, dims = c(length(bxlev), length(wks)), dimnames = list(bxlev, 1:length(wks)))

    ## combine everything via logistic function into probability like score
    .logistic = function(x, x0) 1/(1+exp(-x))

    ## rescale all width based matrices by median bxwidth
                                        #  mbw = median(bxwid$wid)
                                        # rovc.mat = sweep(rovc.mat, 1, mbw, "/")
                                        #  neg.mat = sweep(neg.mat, 1, mbw, "/")

    ## instead just by strobe width for the neg.mat (overlap
    neg.mat = sweep(neg.mat, 1, bthresh/4, "/")

    if (verbose)
      message("Converting scores to quasi probabilities")

    ## transform everything by logistic and sweep rowSums
    ##neg.mat = sweep(neg.mat, 1, apply(neg.mat, 1, min), '-')

    provd = 2*(.logistic(rovd.mat)-0.5)
    ## provc = 2*(.logistic(rovc.mat)-0.5)
    provc = rovc.mat
                                        #  provc = t(apply(as.matrix(provc), 1, function(x) x == max(x))) + 0
    pneg = 2*(.logistic(neg.mat)-0.5)

    if (any( ix <- rowSums(provc)==0)) ## reset blank rows to flat uniform dist
      provc[ix, ] = 1/ncol(provc)

    if (any( ix <- rowSums(provd)==0)) ## reset blank rows to flat uniform dist
      provd[ix, ] = 1/ncol(provd)

    if (any(ix <- rowSums(pneg)==0)) ## reset blank rows to flat uniform dist
      pneg[ix, ] = 1/ncol(pneg)

##    sc = provd*provc*(1-pneg)
##    sc = sweep(sc, 1, rowSums(sc), '/')
    sc = rovc.mat

    ## ## NA all rows that are equivalently distributed across all walks
    ## sc[apply(sc, 1, function(x) all(diff(x)==0)), ] = NA

    ## if (!is.null(wk.nm))
    ##   colnames(sc) = wk.nm

    if (raw)
      return(list(sc = sc, rsc = rsc, provd = provd, provc = provc, pneg = pneg))

    scr = colSums(sc)

    return(scr)
}




#' @name write.htab
#' @title write.htab
#'
#' @description
#' writes data frame (or anything castable to data frame) to pretty HTML formatted table to a static location.
#' Useful for quick table inspection via Chrome or other web browser.
#'
#' Very similar syntax to write.tab
#'
#' @param tab data.frame, data.table, or GRanges
#' @param file optional output .html file (by default ~/public_html/htab.html but directory can be set using env variable HTAB.PATH)
#' @param title title to the page (=NULL)
#' @param footer footer to add to page (=NULL)
#' @param highlight optional integer vector specifying what rows to highlight
#' @param row.names logical flag whether to include row labels (=TRUE)
#' @param col.names logical flag whether to include col labels (=TRUE)
#' @param high.color chaoracter highlight color (= 'yellow')
#' @param row.colors length 2 character vector to shade data rows (= c('lightgray', 'white'))
#' @param header.colors length 2 character vector specifygin background and text for header row (= c('#4A4A4A', 'white'))
#' @param data.size integer font size in px for data, title, and footer (= 15)
#' @param title.size integer font size in px for title (= 15)
#' @param footer.size integer font size in px for footer (= 15)
#' @param header.size integer font size in px for header (= 15)
#' @author Marcin Imielinski
#' @export
write.htab = function(tab, file = NULL,
    title = NULL, # text to be written in bold above the table
    footer = NULL, # text to be writen in bold below the table
    highlight = NULL,  #vector of row indices of the table to highlight
    row.names = TRUE,  # includes row labels
    col.names = TRUE, # includes col labels
    high.color = 'yellow', # highlight color to use
    row.colors = c('lightgray', 'white'), # alternating colors to shade data rows
    header.colors = c('#4A4A4A', 'white'), # two element vector specifying background and text colors for header row, respectively,
    data.size = 15, # font size in px for data, title, and footer
    dt = TRUE,
    force = FALSE, ## force filename argument
    embed = FALSE,
    title.size = 15, footer.size = 20, header.size = round(1.1*data.size))
  {


    # require(hwriter)
    # require(gplots)

    if (is.data.table(tab))
      tab = as.data.frame(tab)

    if (!is.data.frame(tab))
      tab = as.data.frame(tab)

    if (is.null(rownames(tab)))
      row.names = F;

      if (!force)
      {
          if (!is.null(file))
          {
              if (!grepl('(^~)|(^\\/)', file))
                  file = paste('~/public_html/', file, sep = '')
          }
          else
          {
              if (nchar(Sys.getenv('HTAB.PATH'))>0)
                  file = Sys.getenv('HTAB.PATH')
              else
                  file = '~/public_html/htab.html'
          }
      }

      if (nchar(dirname(file))==0)
          file = paste0('./', file)

      if (!file.exists(dirname(file)))
          system(paste('mkdir -p', dirname(file)))

      file = paste(normalizePath(dirname(file)), basename(file), sep = '/')

     if (dt)
     {
         wij(DT::datatable(tab,
                           escape = FALSE,
                           filter = 'top',
                           options = list(
                               pageLength = 100),
                           rownames = row.names), file, embed = embed)
     }
     else
     {

         for (nm in names(tab))
             tab[[nm]] = as.character(tab[[nm]])
         tab[is.na(tab)] = '';
         tab = tab[1:nrow(tab), , drop = FALSE];  #not sure why this is necessary, but deflects occasional weird R bug

         if (any(lix <<- sapply(names(tab), function(x) is.list(tab[, x]))))
             for (i in which(lix))
                 tab[, i] = sapply(tab[, i], function(x) paste(x, collapse = ','))

         dir.create(dirname(normalizePath(dirname(file))), recursive=TRUE, showWarnings = FALSE)
         p = hwriter::openPage(file, link.css = 'hwriter.css')
         if (!is.null(title))
             hwriter::hwrite(title, p, style = sprintf('font-weight:bold; font-size:%spx; margin-top;50px', title.size), center = TRUE, div = TRUE, br = TRUE);

         row.bgcolor = as.list(as.character(gplots::col2hex(row.colors)[(1:nrow(tab))%%length(row.colors)+1]));
         names(row.bgcolor) = rownames(tab)
         if (!is.null(highlight))
             row.bgcolor[rownames(tab[highlight,, drop = FALSE])] = list(gplots::col2hex(high.color));

         row.bgcolor = c(gplots::col2hex(header.colors[1]), row.bgcolor)

                                        #    if (row.names)
         col.bgcolor = gplots::col2hex(header.colors[1])

         col.style = sprintf('font-weight:bold; font-size:%spx; color:%s; text-align:center', header.size, gplots::col2hex(header.colors[2]));

         row.style = rep(sprintf('font-size:%spx; text-align:center', data.size), nrow(tab))
         names(row.style) = rownames(tab)
         row.style = c(list(sprintf('font-weight:bold; font-size:%spx; color:%s; text-align:center', header.size, gplots::col2hex(header.colors[2]))), row.style)

         hwriter::hwrite(tab, p, row.style = row.style, col.style = col.style, col.bgcolor = col.bgcolor, row.names = row.names, col.names = col.names,
                         row.bgcolor = row.bgcolor, table.frame = 'void', table.style = 'margin-left: 30px; margin-top: 30px', br = TRUE)
         if (!is.null(footer))
             hwriter::hwrite(footer, p, style = sprintf('font-weight:bold; text-align:center; font-size:%spx; margin-top;50px', footer.size), center = TRUE, div = TRUE);
         hwriter::closePage(p)
     }
          }

      #' @name col.scale
      #' @title col.scale
      #'
      #' @description
      #' Assigns rgb colors to numeric data values in vector "x".. maps scalar values
      #' in val.range (default c(0,1)) to a linear color scale of between col.min (default white)
      #' and col.max (default black), each which are length 3 vectors or characters.  RGB values are scaled between 0 and 1.
      #'
      #' Values below and above val.min and val.max are mapped to col.max and col.max respectively
      #'
      #' @param x length n numeric or integer data to color
      #' @param val.range data range to assign to colors (= c(0,1))
      #' @param col.min character color to interpolate minimum value in val.range (='white')
      #' @param col.max character color interpolate maximum value in val.range (='black')
      #' @param na.col color to give to na.values (='white')
      #' @param invert logical flag whether to flip min and max (=FALSE)
      #' @author Marcin Imielinski
      #' @return length n vector of colors
      #' @export
      col.scale = function(x, val.range = c(0, 1), col.min = 'white', col.max = 'black', na.col = 'white',
                           invert = FALSE # if T flips rgb.min and rgb.max
                           )
      {
          if (!is.numeric(col.min))
              if (is.character(col.min))
                  col.min = col2rgb(col.min)/255
              else
                  error('Color should be either length 3 vector or character')

          if (!is.numeric(col.max))
              if (is.character(col.max))
                  col.max = col2rgb(col.max)/255
              else
                  error('Color should be either length 3 vector or character')

          col.min = as.numeric(col.min);
          col.max = as.numeric(col.max);

          x = (pmax(val.range[1], pmin(val.range[2], x))-val.range[1])/diff(val.range);
          col.min = pmax(0, pmin(1, col.min))
          col.max = pmax(0, pmin(1, col.max))

          if (invert)
          {
              tmp = col.max
              col.max = col.min
              col.min = tmp
          }

          nna = !is.na(x);

          out = rep(na.col, length(x))
          out[nna] = rgb((col.max[1]-col.min[1])*x[nna] + col.min[1],
          (col.max[2]-col.min[2])*x[nna] + col.min[2],
          (col.max[3]-col.min[3])*x[nna] + col.min[3])

          return(out)
      }


#' @name col.scale
#' @title col.scale
#'
#' @description
#' Assigns rgb colors to numeric data values in vector "x".. maps scalar values
#' in val.range (default c(0,1)) to a linear color scale of between col.min (default white)
#' and col.max (default black), each which are length 3 vectors or characters.  RGB values are scaled between 0 and 1.
#'
#' Values below and above val.min and val.max are mapped to col.max and col.max respectively
#'
#' @param x length n numeric or integer data to color
#' @param val.range data range to assign to colors (= c(0,1))
#' @param col.min character color to interpolate minimum value in val.range (='white')
#' @param col.max character color interpolate maximum value in val.range (='black')
#' @param na.col color to give to na.values (='white')
#' @param invert logical flag whether to flip min and max (=FALSE)
#' @author Marcin Imielinski
#' @return length n vector of colors
#' @export
col.scale = function(x, val.range = c(0, 1), col.min = 'white', col.max = 'black', na.col = 'white',
  invert = FALSE # if T flips rgb.min and rgb.max
  )
  {
    if (!is.numeric(col.min))
      if (is.character(col.min))
        col.min = col2rgb(col.min)/255
      else
        error('Color should be either length 3 vector or character')

    if (!is.numeric(col.max))
      if (is.character(col.max))
        col.max = col2rgb(col.max)/255
      else
        error('Color should be either length 3 vector or character')

    col.min = as.numeric(col.min);
    col.max = as.numeric(col.max);

    x = (pmax(val.range[1], pmin(val.range[2], x))-val.range[1])/diff(val.range);
    col.min = pmax(0, pmin(1, col.min))
    col.max = pmax(0, pmin(1, col.max))

    if (invert)
      {
        tmp = col.max
        col.max = col.min
        col.min = tmp
      }

    nna = !is.na(x);

    out = rep(na.col, length(x))
    out[nna] = rgb((col.max[1]-col.min[1])*x[nna] + col.min[1],
        (col.max[2]-col.min[2])*x[nna] + col.min[2],
        (col.max[3]-col.min[3])*x[nna] + col.min[3])

    return(out)
}

############################
#' @name capitalize
#' @title capitalize
#' @description
#' Capitalize first letter of each character element of vector "string"
#'
#' @param string character vector to capitalize
#' @param un logical flag whether to uncapitalize (=FALSE)
#' @return character vector of strings with capitalized values
#' @export
##############################
capitalize = function(string, un = FALSE)
{
    if (!un)
      {
        capped <- grep("^[^A-Z].*$", string, perl = TRUE)
        substr(string[capped], 1, 1) <- toupper(substr(string[capped],1, 1))
      }
    else
      {
        capped <- grep("^[A-Z].*$", string, perl = TRUE)
        substr(string[capped], 1, 1) <- tolower(substr(string[capped],1, 1))
      }

    return(string)
}

############################
#' @name fisher.pairwise
#' @title fisher.pairwise
#' @description
#' Performs fisher test on cols of matrix / df x vs cols of matrix / df y
#'
#' returns list with ncol(x) by ncol(y) matrices $p and $or denoting the p value and odds ratio of the result of the
#' fisher test on col i of x and col j of y
#'
#' If y is not provided, will correlate rows of x with themselves.
#'
#' @param x n x k1 data frame of categorical data on k1 variables
#' @param y n x k2 optional data frame of categorical data on k2 variables (= x)
#' @return list with field $p and $or correspodning to k1 x k2 matrices of p values and odds ratios for each pair of tests
#' @export
#' @author Marcin Imielinski
#############################
fisher.pairwise = function(x, y = x)
  {
    p = or = matrix(NA, nrow = ncol(x), ncol = ncol(y), dimnames = list(colnames(x), colnames(y)))

    if (nrow(x) != nrow(y))
      stop('x and y must have the same number of rows')

    logical.x = which(sapply(1:ncol(x), function(i) is.logical(x[,i])))
    logical.y = which(sapply(1:ncol(y), function(i) is.logical(y[,i])))

    if (length(logical.x)==0 | length(logical.y)==0)
      warning('No logical columns found')

    for (i in logical.x)
          for (j in logical.y)
            {
              O = table(x[,i], y[,j])
              if (min(dim(O))>1)
                {
                  res = fisher.test(O)
                  or[i, j] = res$estimate
                  p[i, j] = res$p.value
                }
            }

    out = list(p = p, or = or)
    return(out)
  }


################
#' @name strsplit2
#' @title strsplit2
#'
#' @description
#' Strsplit when there are two layers of separators (sep1, sep2) and one needs to extract
#' a collapsed vector of subitem j for all items i.
#'
#' Takes in a character vector and outputs a list of "separated" items
#'
#' @param x character vector
#' @param sep1 character specifying first level separator (=',')
#' @param sep2 character specifying second level separator (=' ')
#' @param j integer specifying which subitem to keep (=1)
#' @author Marcin Imielinski
#' @export
#' @return vector of values for subitem j
################
strsplit2 = function(x, sep1 = ",", sep2 = " ", j = 1)
  {
    return(lapply(strsplit(x, sep1), function(y) sapply(strsplit(y, sep2), function(z) z[[j]])))
  }


################
#' @name timestamp
#' @title timestamp
#'
#' @description
#' returns character time stamp
#' @author Marcin Imielinski
#' @export
################
timestamp = function()
  {
    return(gsub('[\\:\\-]', '', gsub('\\s', '_', Sys.time())))
  }


###############
#' @name img_link
#' @title img_link
#'
#' @description
#' Returns vector of html image links to files "file" with text "caption"
#'
#' if embed = T, then will make img link, and additional arguments may be supplied to image tag (eg height, width)
#' @param file vector of (relative) image file paths to link to
#' @param caption character vector of captions to add (= '')
#' @param embed logical flag whether to imbed images instead of returning links (=FALSE)
#' @param ... additional parameters to embed in tag (e.g. height and width)
#' @return character vector of links (<a> tags) or image tags (<img> or <embed) to dump into an html document
#' @author Marcin Imielinski
#' @export
################
img_link = function(file, caption = NULL, embed = F, ...) {
	if (is.null(caption)) {
		caption = ''
	}

	if (!embed) {
		return(paste('<a href = \"', file, '\">', caption, '</a>', sep = ''))
	} else {
		args = list(...);
		if (length(args)>0) {
			to.quote = is.numeric(unlist(args))+1
			more.args = paste((paste(', ', names(args), " = ", c('\"', '')[to.quote], unlist(args), c('\"', '')[to.quote], sep = "")), collapse="")
		} else {
			more.args = ''
		}

		parts = strsplit(basename(file), "\\.")[[1]]
		file.ext = parts[length(parts)]
		if (file.ext == "tif") { # handles tif images which won't display in chrome with a regular img src tag
			out = paste('<embed src = \"', file, '\", type = "image/tiff"', more.args, ' negative=yes>', sep = "")
		} else {
			out = paste('<img src = \"', file, '\", alt = \"', caption, '\"', more.args, '>', sep = "")
		}

		return(out)
	}
}

#############
#' @name img.html
#' @title img.html
#'
#' @description
#' takes img.paths and dumps out html with imgs +/- names
#'
#' can be dumped into a file for showing many images into a single page
#' alternative to img_link for "embedding images"
#'
#' @param paths vector of (relative) paths to embed in html
#' @param text optional text label to put above embedded images (default = names(paths))
#' @return character vector of img tags
#' @author Marcin Imielinski
#' @export
#############
img.html = function(paths, text = names(paths), height = 1024, width = 768, header = 1)
    {
        if (is.null(text))
            text = ''

        df = data.frame(paths = paths, text = text, height = height, width = width, header = header)
        header = ifelse(!is.na(df$text), paste0("<p> <h", df$header, "> ", df$text, " </h", df$header, "> <p> "), "")
        out = paste0(header,"<img src = '", df$paths, "' height = '", df$height, "' width = '", df$width, "'>", sep = '')
        return(out)
    }

#########
#' @name html_link
#' @title html_link
#'
#' @description
#' returns text with html link
#'
#' @param href character vector of paths to link
#' @param text text to display
#' @return character vector of html link text
#' @author Marcin Imielinski
#' @export
#########
html_link = function(href, text = NULL)
  {
    if (is.null(text))
      text = basename(href)

    return(mapply(function(x,y) html_tag('a', href = x,  text = y), href, text))
  }


#########
#' @name html_tag
#' @title html_tag
#'
#' @description
#' makes a open and close html tag with optional text inside and optional (named) vector of name=value pairs contained inside of opening tag
#'
#' @param tag character vector of tags (without brackets)
#' @param text text to put inside tags
#' @param collapse how to collapse tags (=newline)
#' @return character vector of html
#' @author Marcin Imielinski
#' @export
#########
html_tag = function(tag, text = NULL, collapse = ' ',  ...)
  {
    flags = unlist(list(...))

    if (!is.null(flags))
      {
        if (is.null(names(flags)))
          flag.str = ""
        else
          flag.str = paste(" ", paste(paste(names(flags), paste('"', flags, '"', sep = ""), sep = "=")), collapse = "")
      }
    else
      flag.str = ""

    return(paste(sprintf('<%s%s>', tag, flag.str), paste(text, collapse = collapse), sprintf('</%s>', tag), sep = collapse))
  }


################
#' @name set.comp
#' @title set.comp
#'
#' @description
#' Compares two sets and outputs data frame with "left", "middle", "right" members
#'
#'
#' @author Bryan Hernandez
#' @param s1 vector corresponding to "set 1"
#' @param s2 vector corresponding to "set 2"
#' @return list with fields $left, $middle, and $right corresponding to vectors that are in the left setdiff, intersection, right setdiff respectively
#' @export
################
set.comp = function(s1, s2)
  {
    universe = sort(union(s1, s2));
    out = data.frame(stringsAsFactors = F)
    tmp.comp = list(left = sort(setdiff(s1, s2)), middle = sort(intersect(s1, s2)), right = sort(setdiff(s2, s1)))
    lapply(names(tmp.comp), function(x) out[1:length(tmp.comp[[x]]),x] <<- tmp.comp[[x]])
    out[is.na(out)] = ''
    return(out)
  }



################################
#' @name dedup
#' @title dedup
#'
#' @description
#' relabels duplicates in a character vector with .1, .2, .3
#' (where "." can be replaced by any user specified suffix)
#'
#' @param x input vector to dedup
#' @param suffix suffix separator to use before adding integer for dups in x
#' @return length(x) vector of input + suffix separator + integer for dups and no suffix for "originals"
#' @author Marcin Imielinski
#' @export
################################
dedup = function(x, suffix = '.')
{
  dup = duplicated(x);
  udup = setdiff(unique(x[dup]), NA)
  udup.ix = lapply(udup, function(y) which(x==y))
  udup.suffices = lapply(udup.ix, function(y) c('', paste(suffix, 2:length(y), sep = '')))
  out = x;
  out[unlist(udup.ix)] = paste(out[unlist(udup.ix)], unlist(udup.suffices), sep = '');
  return(out)
}


#' @name nnnusvr
#' @title nnnusvr
#'
#' @description
#' Implementation of non-negative linear nu-support vector regression using CPLEX, given
#' length n response vector y and n x m data matrix X returns a length m
#' weight vector w linear SV regresssion result (or length m+1 result if
#' intercept = TRUE is specified).
#'
#' Additional feature nn (default = TRUE) allows constraining the weight vector to be non-negative
#' 
#'
#' @param X numeric matrix of n data points by m features
#' @param y numeric length(n) response vector
#' @param nu scalar numeric between 0 and 1 representing fraction of data to use as support vectors
#' @param C  scalar non-negative numeric representing the cost
#' @param nn logical scalar whether to constrain the output weights to be non-negative
#' @author Marcin Imielinski
#' @export
nnnusvr = function(X, y, nu = 0.85, C = 1, nn = TRUE, intercept = FALSE)
{
  if (length(nu)!=1)
    stop('nu must be length 1')

  if (length(C)!=1)
    stop('C must be length 1')
  
  if (!(nu>0 & nu<=1))
    stop('nu must be between 0 and 1')

  if (!(C>0))
    stop('C must be non-negative')
    
    ## keep track of variables
    variables =
      rbind(
        data.table(iid = 1:ncol(X), type = 'w', lb = ifelse(nn, 0, -Inf), ub = Inf),
        data.table(iid = 1:length(y), type = 'slack', lb = 0, ub = Inf),
        data.table(iid = 1:length(y), type = 'slack*', lb = 0, ub = Inf),
        data.table(iid = 1, type = 'coef0', lb = ifelse(intercept, -Inf, 0), ub = ifelse(intercept, Inf, 0)),
       data.table(iid = 1, type = 'eps', lb = 0, ub = Inf))
    variables[, name := paste(type, iid, sep = '_')][, id := 1:.N]

    ## keep track of constraints
    constraints =
      rbind(
        data.table(iid = 1:length(y), type = 'tube_ub', b = y, sense = 'L'),
        data.table(iid = 1:length(y), type = 'tube_lb', b = y, sense = 'G'))[, id := 1:.N]

    ## make a sparse Matrix of zeros
    Zero = function(d) sparseMatrix(1, 1, x = 0, dims = d)

    ## initialize A (constraints x variables) matrix to 0
    A = as.matrix(Zero(c(nrow(constraints), nrow(variables))))

    ## weights w
    A[constraints[type == 'tube_ub', id], variables[type == 'w', id]] = X
    A[constraints[type == 'tube_lb', id], variables[type == 'w', id]] = X

    ## intercept b
    A[constraints[type == 'tube_ub', id], variables[name == 'coef0_1', id]] = 1
    A[constraints[type == 'tube_lb', id], variables[name == 'coef0_1', id]] = 1

    ## epsilon
    A[constraints[type == 'tube_ub', id], variables[name == 'eps_1', id]] = -1
    A[constraints[type == 'tube_lb', id], variables[name == 'eps_1', id]] = 1

    ## slack constraintss
    A[cbind(constraints[type == 'tube_ub', id], variables[type == 'slack', id])] = -1
    A[cbind(constraints[type == 'tube_lb', id], variables[type == 'slack*', id])] = 1

    ## quadratic component of objective
    Q = Zero(c(nrow(variables), nrow(variables)))
    Q[cbind(variables[type == 'w', id], variables[type == 'w', id])] = 1

    ## linear component of objective 
    c = rep(0, nrow(variables))
    c[variables[name == 'eps_1', id]] = length(y)*C*nu #C*nu
    c[variables[type == 'slack', id]] = C # C/length(y)
    c[variables[type == 'slack*', id]]= C # C/length(y)
    variables[, c := c]
    
    sol = Rcplex(cvec = c,
                 Amat = A,
                 bvec = constraints$b,
                 Qmat = Q,
                 lb = variables$lb,
                 ub = variables$ub,
                 objsense = 'min',
                 sense = constraints$sense,
                 vtype = 'C')
   
    colnames(A) = variables$name
    variables[ , xopt := sol$xopt]
   
#    constraints[, constraint := paste(arrstring(A), ifelse(sense == 'G', '>=', '<='), constraints$b)]
#    constraints[, xopt := paste(arrstring(A, variables$xopt, signif = 10), ifelse(sense == 'G', '>=', '<='), constraints$b)]

    w = sol$xopt[variables[type == 'w', id]]
    names(w) = colnames(X)
    if (intercept)
      w = c(w, sol$xopt[variables[type == 'coef0', id]])
    return(w)
  }

################################
#' @name is.dup
#' @title is.dup
#'
#' @description
#' labels which vectors in x are part of a dup
#' returns logical TRUE if vector is part of a dup
#'
#' Note: this is a twist on "duplicated" which only returns TRUE if a given element is a duplicate (i.e. duplicated ()  is FALSE
#' for the original version for the duplicate, while is.dup() will be TRUE for that element)
#'
#' x can be vector or matrix
#'
#' @param x vector or matrix to check
#' @return logical vector of length(x) or nrow(x)
#' @author Marcin Imielinski
#' @export
################################
is.dup = function(x)
{
    if (is.matrix(x))
        x = as.data.frame(x)

    if (is.data.frame(x) | data.table::is.data.table(x))
        {
            tmp = x[[1]]
            if (ncol(x)>1)
                for (i in 2:ncol(x))
                    tmp = paste(tmp, x[[i]], sep = '@!$!$!@')
            x = tmp
        }

    d = duplicated(x)
    return(x %in% x[d])
}


########################
#' @name install.packages.bioc
#' @title install.packages.bioc
#'
#' @description
#' shortcut to install bioconductor packages
#'
#' @param pkg character vector of package names to install
#' @author Marcin Imielinski
#' @export
install.packages.bioc = function(pkg)
  {
##    source('http://bioconductor.org/biocLite.R')
    ##    sapply(pkg, biocLite)
    BiocManager::install(pkg)
  }

##########################
#' @name install.packages.github
#' @title install.packages.github
#'
#' @description
#' shortcut to install github packages
#'
#' @param pkg character vector of package names to install
#' @author Marcin Imielinski
#' @export
##########################
install.packages.github = function(pkg, username, branch)
  {
    devtools::install_github(repo = pkg, username = username, branch = branch)
}

####################
#' @name tabstring
#' @title tabstring
#'
#' @description
#' string representation of a named vector (ie the result of tab = table(x)
#
#' ie name1 (value1), name2 (value2), name3 (value3)
#'
#' @param tab "table" or any named(vector)
#' @param sep separator to use between table elements
#' @return character representation of table
#' @export
#' @author Marcin Imielinski
####################
tabstring = function(tab, sep = ', ', sep2 = '_', dt = FALSE)
    {
        if (length(dim(tab))<=1)
            if (dt)
                {
                    tab = data.table(key = names(tab), count = as.numeric(tab))
                    return(tab)
                }
            else
                {
                    return(paste(names(tab), '(', tab, ')', sep = '', collapse = sep))
                }

        else
            {
                # library(reshape)
                mtab = reshape2::melt(tab)
                nm = apply(mtab, 1, function(x) paste(x[-length(x)], collapse = sep2))
                if (dt)
                    {
                        tab = data.table(key = nm, count = mtab$value)
                        setkey(tab, key)
                        return(tab)
                    }
                else
                    {
                        tab = structure(mtab$value, names = nm)
                        return(paste(names(tab), '(', tab, ')', sep = '', collapse = sep))
                    }

            }
    }

####################
#' @name arrstring
#' @title arrstring
#'
#' @description
#' string representation of row array as linear combination of nonzero entries
#' of that row using column names as variables
#' 
#' @param A array
#' @param sep separator to use between table elements
#' @return character representation of table
#' @export
#' @author Marcin Imielinski
####################
arrstring = function(A, x = NULL, sep = ', ', sep2 = '_', signif = 3, dt = FALSE)
{
  if (is.null(dim(A)))
  {
    A = rbind(A)
  }

  if (is.null(colnames(A)))
  {
    colnames(A) = paste0('V', 1:ncol(A))
  }

  if (is.null(x))
  {
    x = colnames(A)
  }
  else
  {
    x = signif(x, signif)
  }

  str = apply(A, 1, function(y) paste(signif(y[y!=0],signif), x[y!=0], sep = '*', collapse = ' + '))

  return(str)    
}  



####################
#' @name dfstring
#' @title dfstring
#'
#' @description
#' "tuple" style chraacter representation of a table, key name1 = value1, name2 = value2
#' either as a single line or many lines
#' useful for quick eyeballing of tabular data
#'
#' @param df data.frame input
#' @param oneline logical flag whether to print on one line (=TRUE)
#' @param sep1 first level separator (=;) i.e. between rows
#' @param sep2 second level separator (=, ) i.e. between columns
#' @return character vector of string representation
#' @author Marcin Imielinski
#' @export
####################
dfstring = function(df, oneline = TRUE,  binary = FALSE, sep1 = '; ', sep2 = ', ')
    {
        if (!class(df)[1]=='data.frame')
            df = as.data.frame(df)

        if (binary)
            return(structure(apply(as.matrix(df), 1, function(x) paste(names(df)[which(as.logical(x))], collapse = sep1)), names = rownames(df)))


        df = as.list(df)

        if (!oneline)
            sep1 = NULL

        if (length(names(df))==0)
            return('')

        nm1 = nm2 = names(df)
        nm1[1] = paste('(', nm1[1], sep = '')
        ## paste inception YIKES
        cmd = paste('paste(paste(', paste('paste("', nm1, '=", df$"', names(df), '", sep = "")', sep = '', collapse = ','),
            ', sep = sep2), ")", sep = "", collapse = sep1)', sep = '')
        return(eval(parse(text=cmd)))
    }



#############################
#' @name levapply
#' @title levapply
#'
#' @description
#' Applies FUN locally to levels of x and returns vector of length()
#' (eg can do a "local" order within levels)
#'
#' @param x input vector of data
#' @param by length(x) vector of categorical labels
#' @param FUN function that takes a length k vector and outputs a length k vector, used for processing each "level" of by
#' @return length(x) vector of outputs, the results of applying FUN to each "by" defined level of x
#' @export
#' @author Marcin Imielinski
#############################
levapply = function(x, by, FUN = 'order')
  {
    if (!is.list(by))
      by = list(by)

    f = factor(do.call('paste', c(list(sep = '|'), by)))
    ixl = split(1:length(x), f);
    ixv = lapply(ixl, function(y) x[y])
    res = structure(unlist(lapply(ixv, FUN)), names = unlist(ixl))
    out = rep(NA, length(x))
    out[as.numeric(names(res))] = res;
    return(out)
  }


## ####################
## #' @name cytoscape
## #' @title cytoscape
## #'
## #' @description
## #' shortcut to connect to local cytoscape instance running on LOCAL.COMPUTER (unix environment variable) via RPC call
## #'
## #' graph must be igraph, or adjacency matrixx
## #'
## #' @param graph igraph or adjacency matrix
## #' @param sessionName character name of cytoscape session to open on local computer
## #' @param host character name of host on which cytoscape is running (=Sys.getenv('LOCAL.COMPUTER'))
## #' @param port integer port on which cytoscape is running (=9000)
## #' @param display logical flag whether to display graph locally (=TRUE)
## #' @param layout character specifying layout to display as (=degree-circle)
## #' @param verbose logical vector whether to use verbose output
## #' @param ... additional arguments to new.CytoscapeWindow
## #' @export
## #' @author Marcin Imielinski
## #' #@importFrom igraph V E E<- get.edgelist list.edge.attributes get.edge.attribute get.vertex.attribute
## cytoscape = function(graph = NULL, sessionName = 'M-ski', host = Sys.getenv('LOCAL.COMPUTER'), port = 9000, display = T, layout = 'degree-circle', verbose = T, ...)
##   {
##     # require(RCytoscape)
##     # require(igraph)

##     if (is(graph, 'matrix'))
##       graph = igraph::graph.adjacency(graph, weighted = 'weight');

##     if (is(graph, 'igraph'))
##       {
##         if (!is.null(E(graph)$weight))
##           E(graph)$weight = 1

##         if (!is.null(E(graph)$arrow.shape))
##           E(graph)$arrow.shape = 'ARROW'

##         graph.nel = igraph2graph(graph)
##       }

##     cw = RCytoscape::new.CytoscapeWindow(sessionName, host = host, rpcPort = port, graph = graph.nel,  ...)

##     if (display)
##       {
##       RCytoscape::displayGraph(cw)
##       RCytoscape::setDefaultBackgroundColor(cw, gplots::col2hex('white'))

##         eG = paste(igraph::get.edgelist(graph)[,1], get.edgelist(graph)[,2], sep = '~')
##         ceG = RCytoscape::cy2.edge.names(cw@graph)

##         if (verbose)
##           cat('Setting line styles\n')

##         if ('line.style' %in% igraph::list.edge.attributes(graph))
##           {
##             uls = setdiff(E(graph)$line.style, NA)
##             RCytoscape::setEdgeLineStyleRule(cw, 'line.style', uls, uls)
##           }

##         if (verbose)
##           cat('Setting arrow shape\n')

##         if (igraph::is.directed(graph))
##           if ('arrow.shape' %in% igraph::list.edge.attributes(graph))
##             RCytoscape::setEdgeTargetArrowRule(cw, 'arrow.shape', unique(E(graph)$arrow.shape), unique(E(graph)$arrow.shape))

##         if (verbose)
##           cat('Setting edge color\n')

##         if ('col' %in% igraph::list.edge.attributes(graph))
##           {
##             uc = setdiff(unique(E(graph)$col), NA);
##             RCytoscape::setEdgeColorRule(cw, 'col', uc, uc, mode = 'lookup')
##           }

##         if (verbose)
##           cat('Setting edge width\n')

##         if ('width' %in% igraph::list.edge.attributes(graph))
##           {
##             uw = setdiff(igraph::E(graph)$width, NA)
##             RCytoscape::setEdgeLineWidthRule(cw, 'width', as.character(uw), uw)
##           }

##         if (verbose)
##           cat('Setting node size\n')

##         if ('size' %in% igraph::list.vertex.attributes(graph))
##           {
##             us = setdiff(unique(igraph::V(graph)$size), NA)
##             RCytoscape::setNodeSizeRule(cw, 'size', us, us, mode = 'lookup')
##           }

##         if (verbose)
##           cat('Setting node color\n')

##         if ('col' %in% igraph::list.vertex.attributes(graph))
##           {
##             uc = setdiff(unique(igraph::V(graph)$col), NA)
##             RCytoscape::setNodeColorRule(cw, 'col', uc, uc, mode = 'lookup')
##           }

##         if (verbose)
##           cat('Setting node labels\n')

##         if ('label' %in% igraph::list.vertex.attributes(graph))
##           RCytoscape::setNodeLabelRule(cw, 'label')

##         if (verbose)
##           cat('Setting node shapes\n')

##         if ('shape' %in% igraph::list.vertex.attributes(graph))
##           {
##             us = setdiff(unique(V(graph)$shape), NA)
##             RCytoscape::setNodeShapeRule(cw, 'shape', us, us, default = 'ELLIPSE')
##           }

##         if (verbose)
##           cat('Setting node width\n')

##         if ('border.width' %in% igraph::list.vertex.attributes(graph))
##           {
##             ubw = setdiff(unique(V(graph)$border.width), NA)
##             RCytoscape::setNodeBorderWidthRule(cw, 'border.width', ubw, ubw)
##           }

##         if (all(c('x', 'y') %in% igraph::list.vertex.attributes(graph)))
##           {
##             good.ix = !is.na(V(graph)$x) & !is.na(V(graph)$y)
##             if (any(good.ix))
##               RCytoscape::setNodePosition(cw, V(graph)$name[good.ix], V(graph)$x[good.ix], V(graph)$y[good.ix])
##           }
##         else
##           RCytoscape::layoutNetwork(cw, layout)


##         RCytoscape::redraw(cw)
##       }

##     return(cw)
##   }


## ####################
## #' @name igraph2graph
## #' @title igraph2graph
## #'
## #' @description
## #' Converts igraph object into janky graphNEL object (for visualization in cytoscape)
## #' and populates all edge features both via the edgeL and as NodeAttributes for visualization
## #' in cytoscape
## #'
## #' #@importFrom graph edgeData<- nodeData<-
## #' @param g igraph object
## #' @author Marcin Imielinski
## #' @export
## #' @return graph object
## igraph2graph = function(g)
##   {
##     # require(igraph)
##     # require(graph)
##     # require(RCytoscape)

##     if (class(V) != 'function' | class(E) != 'function')
##       stop('Namespace conflict - either V() or E() no longer mapping to igraph functions')

##     if (!is.null(V(g)$name))
##       node.labels = V(g)$name
##     else
##       node.labels = as.character(V(g));

##     edge.df = structure(as.data.frame(get.edgelist(g), stringsAsFactors = F), names = c('vertices', 'edges'))
##     if (length(igraph::list.edge.attributes(g))>0)
##       {
##         tmp = do.call('cbind', lapply(igraph::list.edge.attributes(g), function(x) get.edge.attribute(g, x)))
##         colnames(tmp) = igraph::list.edge.attributes(g)
##         edge.df = cbind(edge.df, as.data.frame(tmp, stringsAsFactors = F))
##       }
##     if (!is.null(edge.df$weights))
##       edge.df$weights = as.numeric(edge.df$weights)
##     edge.df[is.na(edge.df)] = "NA"
##     edge.df[,1] = as.character(edge.df[,1])
##     edge.df[,2] = as.character(edge.df[,2])

##     vertex.df = data.frame(vertices = node.labels, stringsAsFactors = F)
##     if (length(igraph::list.vertex.attributes(g))>0)
##       {
##         tmp = do.call('cbind', lapply(igraph::list.vertex.attributes(g), function(x) get.vertex.attribute(g, x)))
##         colnames(tmp) = igraph::list.vertex.attributes(g)
##         vertex.df = cbind(vertex.df, as.data.frame(tmp, stringsAsFactors = F))
##       }

##     ## have to reciprocate edges in undirected otherwise graphNEL will barf
##     if (!igraph::is.directed(g))
##       {
##         edge.df.rev = edge.df;
##         tmp.col = edge.df[,2]
##         edge.df.rev[,2] = edge.df.rev[,1]
##         edge.df.rev[,1] = tmp.col;
##         edge.df = rbind(edge.df, edge.df.rev)
##       }

##     edgeL = lapply(split(edge.df, edge.df$vertices), as.list)[node.labels]
##     names(edgeL) = node.labels;

##     ## retarded GraphNEL object format necessitates these gymnastics
##     null.vert = sapply(edgeL, is.null)
##     blank.edge.item = list(structure(rep(list(c()), length(igraph::list.edge.attributes(g))+1), names = c('edges', igraph::list.edge.attributes(g))))
##     edgeL[null.vert] = blank.edge.item
##     edgemode = c('undirected', 'directed')[1 + igraph::is.directed(g)]

##     out.g = new('graphNEL', node.labels, edgeL, edgemode = edgemode)

##     ## populate edge and node attribute for RCytoscape to access
##     if (ncol(edge.df)>2)
##       {
##         attr.cols = names(edge.df)[3:ncol(edge.df)]
##         for (attr in attr.cols)
##           {
##             if (is.numeric(edge.df[, attr]))
##               out.g = RCytoscape::initEdgeAttribute(out.g, attr, 'numeric', NA)
##             else if (is.integer(edge.df[, attr]))
##               out.g = RCytoscape::initEdgeAttribute(out.g, attr, 'integer', NA)
##             else
##               {
##                 cast.numeric = suppressWarnings(as.numeric(edge.df[, attr]))
##                 cast.character = edge.df[, attr]

##                 if (any(is.na(cast.numeric) & cast.character != "NA"))
##                   {
##                     edge.df[, attr] = as.character(edge.df[, attr])
##                     out.g = RCytoscape::initEdgeAttribute(out.g, attr, 'char', '')
##                   }
##                 else
##                   {
##                     cast.numeric[is.na(cast.numeric)] = 0
##                     edge.df[, attr] = cast.numeric
##                     out.g = RCytoscape::initEdgeAttribute(out.g, attr, 'numeric', '')
##                   }
##               }
##             edgeData(out.g, edge.df[,1], edge.df[,2], attr) = edge.df[,attr]
##           }
##       }

##     if (ncol(vertex.df)>1)
##       {
##         attr.cols = names(vertex.df)[2:ncol(vertex.df)]
##         for (attr in attr.cols)
##           {
##             if (is.numeric(vertex.df[, attr]))
##               out.g = RCytoscape::initNodeAttribute(out.g, attr, 'numeric', NA)
##             else if (is.integer(vertex.df[, attr]))
##               out.g = RCytoscape::initNodeAttribute(out.g, attr, 'integer', NA)
##             else
##               {
##                 cast.numeric = suppressWarnings(as.numeric(vertex.df[, attr]))
##                 cast.character = suppressWarnings(as.character(vertex.df[, attr]))
##                 if (any(setdiff(is.na(cast.numeric) & cast.character != "NA", NA)))
##                   {
##                     vertex.df[, attr] = as.character(vertex.df[, attr])
##                     out.g = RCytoscape::initNodeAttribute(out.g, attr, 'char', '')
##                   }
##                 else
##                   {
##                     cast.numeric[is.na(cast.numeric)] = 0
##                     vertex.df[, attr] = cast.numeric
##                     out.g = RCytoscape::initNodeAttribute(out.g, attr, 'char', '')
##                   }
##               }
##             nodeData(out.g, vertex.df[,1], attr) = vertex.df[,attr]
##           }
##       }

##     return(out.g)
##   }

## #' @name cyto2igraph
## #' @title cyto2igraph
## #'
## #' @description
## #' Pulls graph in CytoscapeWindow "cw" and returns as igraph object
## #'
## #' @param cw CytoscapeWindow object to grab from (see Rcytoscape)
## #' @author Marcin Imielinski
## #' @return igraph object
## #' @export
## cyto2igraph = function(cw)
## {
##   node.attr = RCytoscape::getAllNodeAttributes(cw)
##   edge.attr = RCytoscape::getAllEdgeAttributes(cw)
##   directed = graph::edgemode(cw@graph) == 'directed'

##   edge.df = cbind(from = edge.attr$source, to = edge.attr$target,
##     edge.attr[, setdiff(names(edge.attr), c('source', 'target'))])

##   if ('name' %in% node.attr)
##     node.df = node.attr[, c('name', setdiff(names(node.attr), 'name'))]
##   else
##     node.df = cbind(name = rownames(node.attr), node.attr)

##   return(igraph::graph.data.frame(edge.df, directed = directed, vertices = node.df))
## }

#' @name brewer.master
#' @title brewer.master
#'
#' @description
#' Makes a lot of brewer colors using an "inexhaustible" brewer palette ie will not complain if number of colors requested is too high.
#'
#' Yes - this technically violates the "grammar of graphics", but meant for quick and dirty use.
#'
#' @param n TODO
#' @param palette character specifyign pallette to start with (options are: Blues, BuGn, BuPu, GnBu, Greens Greys, Oranges, OrRd, PuBu, PuBuGn, PuRd, Purples, RdPu, Reds, YlFn, YlFnBu, YlOrBr, YlOrRd, BrBg, PiYG, PRGn, PuOr, RdBu, RdGy, RdYlBu, RdYlGn, Spectral, Accent, Dark2, Paired, Pastel1, Pastel2, Set2, Set3)
#' @return length(n) character vector of colors
#' @author Marcin Imielinski
#' @export
brewer.master = function(n, palette = NULL, wes = FALSE,  list = FALSE)
{
    if (wes)
    {
      palettes = c("Royal2"=5, "Chevalier1"=4, "Darjeeling1"=5, "IsleofDogs1"=6, "Darjeeling2"=5, "Moonrise1"=4, "BottleRocket1"=7, "Rushmore"=5, "Moonrise3"=5, "Cavalcanti1"=5, "Rushmore1"=5, "FantasticFox1"=5, "BottleRocket2"=5, "Royal1"=4, "IsleofDogs2"=5, "Moonrise2"=4, "GrandBudapest1"=4, "GrandBudapest2"=4, "Zissou1"=5)
    }
    else
    {
    palettes = list(
      sequential = c('Blues'=9,'BuGn'=9, 'BuPu'=9, 'GnBu'=9, 'Greens'=9, 'Greys'=9, 'Oranges'=9, 'OrRd'=9, 'PuBu'=9, 'PuBuGn'=9, 'PuRd'=9, 'Purples'=9, 'RdPu'=9, 'Reds'=9, 'YlGn'=9, 'YlGnBu'=9, 'YlOrBr'=9, 'YlOrRd'=9),
      diverging = c('BrBG'=11, 'PiYG'=11, 'PRGn'=11, 'PuOr'=11, 'RdBu'=11, 'RdGy'=11, 'RdYlBu'=11, 'RdYlGn'=11, 'Spectral'=11),
          qualitative = c('Accent'=8, 'Dark2'=8, 'Paired'=12, 'Pastel1'=8, 'Pastel2'=8, 'Set1'=9, 'Set2'=8, 'Set3'=12)
        );
      }

  palettes = unlist(palettes);
  if (list)
    return(palettes)


  if (is.null(palette))
    palette = names(palettes)[1]

  nms = NULL
    if (is.character(n) | is.factor(n))
    {
        nms = unique(n)
        n = length(nms)
    }
  
    names(palettes) = gsub('\\w+\\.', '', names(palettes))

    if (palette %in% names(palettes))
      i = match(palette, names(palettes))
    else
      i = ((max(c(1, suppressWarnings(as.integer(palette))), na.rm = T)-1) %% length(palettes))+1

    col = c();
    col.remain = n;

    while (col.remain > 0)
    {
      if (col.remain > palettes[i])
      {
        next.n = palettes[i]
        col.remain = col.remain-next.n;
      }
      else
      {
        next.n = col.remain
        col.remain = 0;
      }

      if (!wes)
        {
          col = c(col, RColorBrewer::brewer.pal(max(next.n, 3), names(palettes[i])))
        }
      else
      {
        col = c(col, wesanderson::wes_palettes[[names(palettes[i])]])
      }

      i = ((i) %% length(palettes))+1
    }

    col = col[1:n]
    names(col) = nms
    return(col)
}

#' @name charToDec
#' @title charToDec
#'
#' @description
#' converts character vector to byte vector in decimal representation
#' @param c character vector
#' @return length(c) integer vector of byte representation of c
#' @author Marcin Imielinski
#' @export
charToDec = function(c)
  {
    return(as(charToRaw(c), 'integer'))
  }

#' @name which.char
#' @title which.char
#'
#' @description
#' finds the index of the character in subject (length 1 character vector) matching
#' nchar = 1 single character query
#' eg which.char('a', 'cat') = 2
#'
#' if query has more than one char (or has length>1) then will return indices matching <any one> of the characters in any
#' element of query
#' @param subject length 1 character vector
#' @param query length 1, nchar 1 character
#' @return indices in subject that query appears
#' @author Marcin Imielinski
#' @export
which.char = function(subject, query)
  {
    if (length(query)>1)
      query = paste(query, collapse = '')

    which(charToRaw(subject[1]) %in% charToRaw(query))
  }

#' @name modix
#' @title modix
#'
#' @description
#' Takes integer input ix and projects on to 1-based modulus over base l
#'
#' ie modix(1, 5) -> 1, modix(5, 5) -> 5, modix(6, 5) -> 1
#'
#' @param ix input indices to apply module
#' @param l base of ix
#' @return ((ix-1) mod l) - 1
#' @author Marcin Imielinski
#' @export
modix = function(ix, l)
  {
    return(((ix-1) %% l)+1)
  }

#' @name elcycles
#' @title elcycles
#'
#' @description
#' enumerates all elementary cycles in a graph via igraph library
#'
#' A is either an adjacency matrix or igraph object
#'
#' @param A adjacency matrix
#' @return
#' list with fields:
#' $cycles = list of vertices in elementary cycles
#' $cycles.eix = list of edges in elementary cycles, where edges are numbered according to the 1D index of adj matrix A
#' @export
#' @author Marcin Imielinski
elcycles = function(A)
  {
    if (inherits(A, 'igraph'))
      A = as.matrix(igraph::graph.adjacency(A));

    A = abs(sign(A)) * matrix(1:length(A), nrow = nrow(A))

    # list of cycles (ie lists of node indices, seeded with self cycles
    out = lapply(which(diag(A)!=0), function(x) x)
    cl = igraph::clusters(igraph::graph.adjacency(A!=0), 'strong')

    while (length(cl.left <- which(cl$csize>1))>0)
      {
        # get all cycles
        seeds = match(cl.left, cl$membership)
        cycles.eix = cycles = list();

        for (s in seeds)
          {
            tail = which(A[s, ]!=0); # this is vector containing tail item of these.cycles
            these.cycles = lapply(tail, function(x) x)
            these.cycles.eix = lapply(1:length(tail), function(x) list())
            visited = array(FALSE, dim = c(length(these.cycles), ncol(A)))   ## visited is cycles x nodes matrix keeps track of visited nodes in each path
            done = rep(FALSE, length(these.cycles)) ## cycles are done if their tail node comes back to "s"
            while (any(!done)) ## cycle through !done cycles
              {
                print(these.cycles)
                j = which(!done)[1]; # pick next !done path
                i = tail[j]  # i is the last element of that path
                visited[j, i] = TRUE;

                if (i != s)
                  {
                    children = which(A[i, ] != 0 & !visited[j, ])
                    children.eix = A[i, children]

                    if (length(children)>0)
                      {
                        these.cycles = c(these.cycles[-j], lapply(children, function(x) c(these.cycles[[j]], x)))
                        these.cycles.eix = c(these.cycles.eix[-j], lapply(children.eix, function(x) c(these.cycles.eix[[j]], x)))
                        visited = rbind(visited[-j, , drop = FALSE], visited[rep(j, length(children)), , drop = FALSE])
                        tail = c(tail[-j], children);
                        done = c(done[-j], done[rep(j, length(children))]);
                      }
                    else ## remove this path since it is childless and does not end in s
                      {
                        these.cycles = these.cycles[-j]
                        these.cycles.eix = these.cycles.eix[-j]
                        visited = visited[-j, , drop = FALSE]
                        tail = tail[-j];
                        done = done[-j];
                      }
                  }
                else
                  done[j] = TRUE ## yes we have found a cycle
              }

            # collect cycles
            cycles = c(cycles, these.cycles);
            cycles.eix = c(cycles.eix, these.cycles.eix)

            # recompute strongly connected components (zeroing out seeds in A matrix);
            A[seeds, ] = 0
            A[, seeds] = 0
            cl = igraph::clusters(igraph::graph.adjacency(A), 'strong')
          }
      }
  }



#########################
#' @name rmix
#' @title rmix
#'
#' @description
#'
#' sample N points from a mixture of k densities of a single functional form (eg norm, beta, multinomial)
#' where n is either an integer vector of length k denoting how many samples to be drawn from each density
#' (in which case N = sum(n)) or n is a scalar, in which case n points are drawn from each density and N = n*k.
#'
#' p = params data frame whose named columns correspond to arguments to rdens (eg $n, $shape1, $shape2 for rbeta or $n, $mean, $sd for rnorm)
#' rdens = function encoding random number generator for given density, that takes as input named columns of params
#' n = either an nrow(p) integer vector or scalar denoting how many samples to draw from each density
#'
#' n can also be just be a column of p
#'
#' useful for plotting "smears" of points
#'
#' Output is the rbind-ed output of individual rdens calls
#'
#' @param p k x p data frame of k parameter sets of rdens density functions, each column is a parameter value, each row is a parameter setting for a mixture component
#' @param rdens R density specific random number generator function object (eg rnorm)
#' @param n length k or legnth 1 integer specifying number of samples to draw from each mixture component
#' @export
#' @author Marcin Imielinski
########################
rmix = function(p, rdens, n = NULL)
  {
    if (!is.null(n))
      p$n = n

    tmp = do.call('mapply', c(list(FUN = rdens, SIMPLIFY = FALSE), as.list(p)))

    if (length(tmp)>0)
      if (!is.null(dim(tmp[[1]])))
        return(do.call('rbind', tmp))
      else
        return(do.call('c', tmp))
    else
      return(tmp)
  }

#' @name dmix
#' @title dmix
#'
#' @description
#' generates data frame of density points in a provided range for a provided mix of k densities of a singlen family
#' useful for plugging into downstream plotting (eg ggplot 2)
#'
#' "..." variables depend on density function, arguments should be provided as they would to the
#' corresponding R function (ie with respect to vectorization)
#'
#' if collapse = TRUE then the density will be summed according to the mixing parameter yielding a single density
#' (ie a fuzzy histogram) summarizing the mixing distribution
#'
#' @param dens character specifying R density function, the possibilities include (with ... arguments shown alongside the density names)
#' dnorm: mean, sd
#' dbinom: size, prob
#' dmultinom: size, prob
#' dgamma: shape, rate
#' dbeta: shape1, shape2
#' @param ... additional density specific arguments each with vectorized values of length k, where k is the number of desired mixture componetns, see dens arugment)
#' @param alpha length(k) numeric vector specifying mixing probability
#' @param xlim length 2 vector specifying plot bounds (=NULL)
#' @param n integer number of points to draw distribution over (=500)
#' @param plot logical flag specifying whether to draw the plot (=FALSE)
#' @param fill logical flag specifying whether to fill the colored plots (=TRUE)
#' @param collapse collapse logical flag whether to collapse the mixture components into a single mixture)
#' @return if plot == TRUE then ggplot2 object of plot, otherwise data.frame of data points with fields $id specifying thee mixture id, $x = data value,
#' @author Marcin Imielinski
#' @export
dmix = function(dens = 'dnorm', xlim = NULL, n = 500, alpha = NULL, plot = F, fill = T, collapse = F,  ...)
  {
    if (is.null(xlim))
      {
        if (dens == 'dnorm')
          xlim = c(min(list(...)$mean-3*list(...)$sd, na.rm = T), max(list(...)$mean+3*list(...)$sd, na.rm = T))
        else if (dens == 'dbinom')
          xlim = c(0, max(list(...)$size))
        else if (dens == "dbeta")
          xlim = c(0, 1)
        else if (dens == 'dmultinom')
          xlim = c(0, max(list(...)$size))
        else if (dens == "gamma")
          {
            mean = list(...)$alpha/list(...)$beta
            sd = sqrt(list(...)$alpha/list(...)$beta^2)
            xlim = c(mean-3*sd, mean+3*sd)
          }
      }

    args = do.call('data.frame', list(...));

    if (is.null(alpha))
      alpha = rep(1, nrow(args))/nrow(args)

    x = seq(xlim[1], xlim[2], length = n);

    id = prob = NULL ## NOTE fix
    out = data.frame(id = rep(1:nrow(args), each = length(x)),  x = rep(x, nrow(args)), prob = unlist(lapply(1:nrow(args), function(i)
         d = alpha[i]*do.call(dens, c(list(x=x), as.list(args[i, ]))))))

    if (collapse)
      {
          out = as.data.table(out)[, list(prob = sum(prob), id = 1), by = x]
#          out = aggregate(prob ~ x, data = out, FUN = sum)
#        out$id = 1;
      }

    if (plot)
      {
        if (fill)
          return(ggplot2::ggplot(out, aes(x = x, y = prob, group = id, fill = id)) + geom_ribbon(alpha = 0.3, color = 'black', aes(ymin = 0, ymax = prob)) + scale_fill_gradient(low = 'red'))
        else
          return(ggplot2::ggplot(out, aes(x = x, y = prob, group = id, color = id)) + geom_line(alpha = 0.3) + scale_color_gradient(low = 'red'))
      }
    else
      return(out)
  }

#################################
# svec
#
# makes "sparsely" defined numeric vector of length n
# using name= value arguments
#
# eg svec(10, '5'=c(1,2,3,4)) makes vector of length 10 with 1,2,3,4 having value 5.
# Note: numeric keys have to be enclosed in quotes
#
# conflicts (ie values hitting the same index) are resolved with FUN (eg sum = adds)
# (similar to matlab accumarray)
#
# dval is "default" val for unspecified indices
#################################
svec = function(n = 0, dval = 0, op = '+', ...)
  {
    args = list(...);
    n = pmax(n, sapply(args, function(x) max(x, na.rm = T)))
    out = blank = rep(dval, n);
    for (a in names(args))
      {
        tmp = blank;
        tmp[args[[a]]] = as.numeric(a);
        out = do.call(op, list(out, tmp))
      }
    return(out)
  }

##################################
#' @name nz
#' @title nz
#'
#' @description
#' outputs the nonzero entries of a vector or array
#'
#' @param x length(x)
#' @param zero integer specifying what to use as the "zero" value in the input (=0)
#' @author Marcin Imielinski
#' @return data.frame of row id col id value pairs
#' @export
##################################
nz = function(x, zero = 0, full = FALSE, matrix = TRUE)
  {
    if (length(x)==0)
      return(data.frame())

    if (is.null(nrow(x)))
      x = as.matrix(x)

    if (is.matrix(x))
    {
      ix = which(x!=zero, arr.ind = T)
    }
    else if (is(x, 'Matrix'))
      {
        ix = Matrix::which(x!=zero, arr.ind = T)
      }
    else
    {
      stop('x should be matrix or Matrix')
    }

    if (nrow(ix)==0)
      return(data.frame())

    out = NULL;

    if (!is.null(rownames(x)))
      out = data.frame(rowname = rownames(x)[ix[,1]], stringsAsFactors = F)

    if (is.null(out))
      out = data.frame(row = ix[,1])
    else
      out = cbind(out, data.frame(row = ix[,1]))

    if (!is.null(colnames(x)))
      out = cbind(out, data.frame(colname = colnames(x)[ix[,2]], stringsAsFactors = F))

    out = cbind(out, data.frame(col = ix[,2], val = x[ix], stringsAsFactors = F))
    rownames(out) = NULL

    if (matrix)
        {
            tmp = out;
            if ('rowname' %in% colnames(tmp))
                rid = factor(tmp[, 'rowname'])
            else
                rid = factor(tmp[, 'row'])

            if ('colname' %in% colnames(tmp))
                cid = factor(tmp[, 'colname'])
            else
                cid = factor(tmp[, 'col'])

            ##out = sparseMatrix::sparseMatrix(as.integer(rid), as.integer(cid), x = tmp[, 'val'], dimnames = list(levels(rid), levels(cid)))
            ## sparseMatrix not avaiable on R-3.2
            out = sparseMatrix(as.integer(rid), as.integer(cid), x = tmp[, 'val'], dimnames = list(levels(rid), levels(cid)))

            if (full)
                out = as.matrix(out)
        }

    return(out)
  }


###############################################
#' @name sparse_subset
#' @title sparse_subset
#'
#' @description
#' given k1 x n matrix A and k2 x n matrix B
#' returns k1 x k2 matrix C whose entries ij = 1 if the set of nonzero components of row i of A is
#' a (+/- strict) subset of the nonzero components of row j of B
#'
#' @param A k1 x n matrix
#' @param B k2 x n matrix
#' @param strict logical flag whether to return strict subset (=FALSE)
#' @param chunksize integer size of rows to process from each matrix at a single iteration (=100)
#' @param quiet logical flag (=FALSE)
#' @return k1 x k2 matrix C whose entries ij = 1 if the set of nonzero components of row i of A is
#' a (+/- strict) subset of the nonzero components of row j of B
#' @export
#' @author Marcin Imielinski
###############################################
sparse_subset = function(A, B, strict = FALSE, chunksize = 100, quiet = FALSE)
  {
    nz = colSums(A!=0, 1)>0
    if (is.null(dim(A)) | is.null(dim(B)))
      return(NULL)

    C = sparseMatrix(i = c(), j = c(), dims = c(nrow(A), nrow(B)))

    for (i in seq(1, nrow(A), chunksize))
      {
        ixA = i:min(nrow(A), i+chunksize-1)
        for (j in seq(1, nrow(B), chunksize))
          {
            ixB = j:min(nrow(B), j+chunksize-1)

            if (length(ixA)>0 & length(ixB)>0 & !quiet)
              cat(sprintf('\t interval A %s to %s (%d) \t interval B %d to %d (%d)\n', ixA[1], ixA[length(ixA)], nrow(A), ixB[1], ixB[length(ixB)], nrow(B)))
            if (strict)
              C[ixA, ixB] = (sign((A[ixA, , drop = FALSE]!=0)) %*% sign(t(B[ixB, , drop = FALSE]!=0))) * (sign((A[ixA, , drop = FALSE]==0)) %*% sign(t(B[ixB, , drop = FALSE]!=0))>0)
            else
              C[ixA, ixB] = (sign(A[ixA, nz, drop = FALSE]!=0) %*% sign(t(B[ixB, nz, drop = FALSE]==0)))==0
          }
      }

    return(C)
  }


######################################################
#' @name morder
#' @title morder
#'
#' @description
#' matrix order wrt columns ..
#' ie ordering rows matrix based on left to right ordering of columns (if MARGIN = 1)
#' OR  ordering columns of  matrix based on top to bottom ordering of rows (if MARGIN = 2)
#'
#' @param A matrix of values
#' @param orient integer orientation, if 1 will do row-wise ordering, otherwise column ordering (=1)
#' @return input matrix with rows and columns ordered
#' @export
#' @author Marcin Imielinski
######################################################
morder = function(A, orient = 1)
  {
    if (orient==1)
      return(do.call('order', lapply(1:ncol(A), function(x) A[,x])))
    else
      return(do.call('order', lapply(1:nrow(A), function(x) A[x,])))
  }


######################################################
#' @name mmatch
#' @title mmatch
#'
#' @description
#' match rows of matrix A to matrix B
#'
#' @param A query matrix k1 x n
#' @param B subject matrix k2 x n
#' @param dir 1
#' @return length k1 vector specifying first row of B matching row i of A
#' @export
#' @author Marcin Imielinski
######################################################
mmatch = function(A, B, dir = 1, default.value = 0)
{
  nzix = which(A!=default.value, arr.ind = TRUE)
  Adt = as.data.table(nzix)[, v := A[nzix]]
  if (dir == 2)
    setnames(Adt, c('row', 'col'), c('col', 'row'))
  sA = Adt[, paste(col, v, collapse = ' '), by = row]
  setkey(sA, row)

  nzix = which(B!=default.value, arr.ind = TRUE)
  Bdt = as.data.table(nzix)[, v := B[nzix]]
  if (dir == 2)
    setnames(Bdt, c('row', 'col'), c('col', 'row'))
  sB = Bdt[, paste(col, v, collapse = ' '), by = row]
  setkey(sB, V1)

  ix = sB[.(sA[.(1:nrow(A)), ]$V1), unname(row)]

  return(ix)
}

##########################################################
#' @name bisort
#' @title bisort
#'
#' @description
#' "bisorts" matrix according to rows and columns (and optionally removes empty rows, ie with no nonzero)
#'
#' @param A matrix to sort
#' @param drop logical flag whether to drop empty rows (=FALSE)
#' @export
#' @author Marcin Imielinski
##########################################################
bisort = function(A, drop = F)
  {
    if (drop)
      A = A[which(rowSums(A!=0)>0), , drop = FALSE]

    A = A[, morder(t(A)), drop = FALSE];
    A = A[morder(A), , drop = FALSE];
    return(A)
  }

##############################################################
#' @name setxor
#' @title setxor
#'
#' @param A vector specifying set A
#' @param B vector specifying set B
#' @export
#' @author Marcin Imielinski
#' @return elements in A or B that are not in the intersection of A and B
##############################################################
setxor = function(A, B)
  {
    return(setdiff(union(A,B), intersect(A,B)))
  }


##############################################################
#' @name sub2ind
#' @title sub2ind
#'
#' @description
#' MATLAB style sub2ind function in R physical essence.  Provides the one dim matrix index
#' of row-column locations in matrix
#'
#' (RIP matlab)
#'
#' @param dim dimension of matrix to return index for
#' @param r integer vector of row index to look up
#' @param c length(r) integer vector of column index to look up
#' @param byrow whether to calculate indices by row or column (= FALSE)
#' @return length(r) vector of 1D indices into matrix with dim "dim"
#' @author Marcin Imielinski
#' @export
##############################################################
sub2ind = function(dim, r, c, byrow = F) if (byrow) (r-1)*dim[2] + c else (c-1)*dim[1] + r

##############################################################
#' @name ind2sub
#' @title ind2sub
#'
#' @description
#' MATLAB style ind2sub function in R physical essence.  Provides the 2D row / column index for a
#' given 1D query
#'
#' @param dim dimensions of matrix to query
#' @param ind 1D index
#' @return length(ind) x 2 matrix of row and column index pairs corresponding to input ind in dim "dim" matrix
#' @param byrow whether to calculate indices by row or column (= FALSE)
#' @author Marcin Imielinski
#' @export
##############################################################
ind2sub= function(dim, ind, byrow = F)
  {
    if (byrow)
      cbind(floor((ind-1) / dim[2])+1, ((ind-1) %% dim[2]+1))
    else
      cbind(((ind-1) %% dim[1])+1, floor((ind-1) / dim[1])+1)
  }

#############################################################
#' @name munlist
#' @title munlist
#'
#' @description
#' unlists a list of vectors, matrices, data frames into a n x k matrix
#' whose first column specifies the list item index of the entry
#' and second column specifies the sublist item index of the entry
#' and the remaining columns specifies the value(s) of the vector
#' or matrices.
#'
#' force.cbind = T will force concatenation via 'cbind'
#' force.rbind = T will force concatenation via 'rbind'
#'
#' @param x list of vectors, matrices, or data frames
#' @param force.rbind logical flag to force concatenation via rbind (=FALSE), otherwise will guess
#' @param force.cbind logical flag to force concatenation via cbind (=FALSE), otherwise will guess
#' @param force.list logical flag to force concatenation via unlist (=FALSE), otherwise will guess
#' @return data.frame of concatenated input data with additional fields $ix and $iix specifying the list item and within-list index from which the given row originated from
#' @author Marcin Imielinski
#' @export
#############################################################
munlist = function(x, force.rbind = F, force.cbind = F, force.list = F)
  {
    if (!any(c(force.list, force.cbind, force.rbind)))
      {
        if (any(sapply(x, function(y) is.null(dim(y)))))
          force.list = T
        if (length(unique(sapply(x, function(y) dim(y)[2]))) == 1)
          force.rbind = T
        if ((length(unique(sapply(x, function(y) dim(y)[1]))) == 1))
          force.cbind = T
      }
    else
      force.list = T

    if (force.list)
      return(cbind(ix = unlist(lapply(1:length(x), function(y) rep(y, length(x[[y]])))),
                   iix = unlist(lapply(1:length(x), function(y) if (length(x[[y]])>0) 1:length(x[[y]]) else NULL)),
                   unlist(x)))
    else if (force.rbind)
      return(cbind(ix = unlist(lapply(1:length(x), function(y) rep(y, nrow(x[[y]])))),
                   iix = unlist(lapply(1:length(x), function(y) if (nrow(x[[y]])>0) 1:nrow(x[[y]]) else NULL)),
                   do.call('rbind', x)))
    else if (force.cbind)
      return(t(rbind(ix = unlist(lapply(1:length(x), function(y) rep(y, ncol(x[[y]])))),
                     iix = unlist(lapply(1:length(x), function(y) if (ncol(x[[y]])>0) 1:ncol(x[[y]]) else NULL)),
                   do.call('cbind', x))))
  }


#############################################################
#' @name dunlist
#' @title dunlist
#'
#' @description
#' unlists a list of vectors, matrices, data.tables into a data.table indexed by the list id
#' $listid
#'
#' does fill = TRUE in case the data.tables inside the list do not have compatible column names 
#' 
#' @param x list of vectors, matrices, or data frames
#' @return data.frame of concatenated input data with additional fields $ix and $iix specifying the list item and within-list index from which the given row originated from
#' @author Marcin Imielinski9
#' @export
#############################################################
dunlist = function (x) 
{
    listid = rep(1:length(x), elementNROWS(x))
    if (!is.null(names(x))) 
        listid = names(x)[listid]
    xu = unlist(x, use.names = FALSE)
    if (is.null(xu)) {
        return(as.data.table(list(listid = c(), V1 = c())))
    }
    if (!(inherits(xu, "data.frame")) | inherits(xu, "data.table")) 
        xu = data.table(V1 = xu)
    out = cbind(data.table(listid = listid), xu)
    setkey(out, listid)
    return(out)
}


############################################################
#' @name readRDA
#' @title readRDA
#'
#' @description
#' loads Rdata environment into a list variable and returns
#' (to mirror RDS functionality)
#'
#' @param fn file name of .rda or .RData file
#' @return object containing all the elements of the environment stored in fn
#' @author Marcin Imielinski
#' @export
############################################################
readRDA = function(fn)
  {
    my.env  = new.env()
    load(fn, my.env);
    return(as.list(my.env))
  }
############################################################
gr.flip = function(...)
  {
      return(gr.flipstrand(...))
  }


#' @name dplot
#' @title dplot
#'
#' @description
#' Plots dotplot of grouped data
#'
#' @param y numeric vector of data
#' @param group length(y) vector of catageories
#' @param ylab y axis label (='')
#' @param xlab x axis label (='')
#' @param log logical flag whether to plot y axis in log format (=FALSE)
#' @param dotsize integer dot size to plot with, as function of 0.02 category width plot real estate (= NULL)
#' @param binwidth numeric binwidth of histogram in units of data quantiles (= NULL)
#' @param title character title of plot (='')
#' @param ylim y limits of plot (= NULL)
#' @param text.size text size of legend (= NULL)
#' @author Marcin Imielinski
#' @export
dplot = function(y, group, ylab = '', xlab = '', log = F, dotsize = NULL, binwidth = 0.02, title = NULL, ylim = NULL, text.size = NULL)
  {

      df = data.frame(y = y, group = as.character(group), stringsAsFactors = F)

      binwidth = as.numeric((quantile(y, c(0.99)) - quantile(y, c(0.01))))*binwidth
      maxstack = max(hist(y, diff(range(y, na.rm = TRUE))/binwidth, plot = FALSE)$counts)

      if (is.null(dotsize)) ## control sizing if ntot specified based on max stack size (which is function of binwidth)
          dotsize = pmin(1, 50/maxstack)


    if (is.null(dotsize))
        g = ggplot(df, aes(x = group, y = y)) + theme_bw() + theme(text = element_text(size = text.size)) + geom_dotplot(binaxis = 'y', method = 'dotdensity', stackdir = 'center', position = 'identity', binwidth = binwidth)
    else
        g = ggplot(df, aes(x = group, y = y)) + theme_bw() + theme(text = element_text(size = text.size)) + geom_dotplot(binaxis = 'y', method = 'dotdensity', stackdir = 'center', position = 'identity', dotsize = dotsize, binwidth = binwidth)

    if (!is.null(ylab))
      g = g + labs(y = ylab)

    if (!is.null(xlab))
      g = g + labs(x = xlab)

    if (!is.null(title))
      g = g + ggtitle(title)

    if (log)
      {
        if (!is.null(ylim))
          g = g + scale_y_log10(limits = ylim)
        else
          g = g + scale_y_log10()
      }
    else if (!is.null(ylim))
      g = g + scale_y_continuous(limits = ylim)

    print(g)
    'voila'
}


#' @name dirr
#' @title dirr
#'
#' @description
#' a variant of dir that gsubs pattern from normal output of dir to name output vector
#'
#' eg dirr(path, '.txt') will return dir output with .txt removed
#' eg dirr(path, '.txt', '.rds' ) will return dir output with .txt subbed with .rds
#'
#' @param x character of path to run dir on
#' @param pattern character pattern to limit files to and to replace with rep
#' @param rep character pattern to replace filenames with
#' @param full whether to return full path
#' @param ... additional arguments to dir
#' @return named vector of file paths, named by file names in dir gsub-stripped with pattern
#' @author Marcin Imielinski
#' @export
dirr = function(x, pattern = NULL, rep = '', full = TRUE,  ...)
  {
      out = dir(x, pattern, full.names = full, ...)
      if (!is.null(pattern))
          names(out) = gsub(pattern, rep, basename(out))
      else
          names(out) = basename(out)
    return(out)
  }


## takes k  n_i x m matrices with n_1, ..., n_k rows and outputs
## a (n_1 * n_2 .. * n_k) x m x k matrix of all k cartesian combinations of
## of the rows of these matrices
##
## first matrix can have 3 dimensions, i.e. be n_1 x m x k_0, in which case
## the additional combinations will be added to the end (i.e. the final
## matrix will havve (n_1 * n_2 * ... * n_k) x  (k_9 + k -1 ) x m combos
.matcart = function(...)
  {
    mats = list(...)
    if (length(mats)==0)
      return(NULL)
    out = mats[[1]]
    if (length(dim(out))==2)
      out = array(out, dim = c(dim(out), 1))
    if (length(mats)==1)
      return(out)
    if (length(mats)>1)
      for (i in 2:length(mats))
        {
          y = mats[[i]]
          ix = cbind(rep(1:nrow(out), nrow(y)), rep(1:nrow(y), each = nrow(out)))
          out = array(c(out[ix[, 1],,], y[ix[,2], ]), dim = c(nrow(ix), dim(out)[2], dim(out)[3]+1))
        }
    return(out)
                                        #        return(aperm(out, c(1, 3, 2)))
  }

###################
#' @name padding
#' @title padding
#'
#' @description
#' pads an (integer) vector with k places below and above its lowest and highest value
#'
#' (by default, clips at 0)
#'
#' useful for querying around specific entires of vector, matrix, data.frame, GRanges ewtc
#'
#' @param x integer vector to pad
#' @param k window around each entry to pad
#' @param clip logical flag whether to clip elements below 0 (=TRUE)
#' @return "padded" integer vector of unique entires with entries in k window around each input included
#' @author Marcin Imielinski
#' @export
###################
padding = function(x, k, clip = T)
  {
    out = unique(as.vector(rbind(
      apply(cbind(x), 1, function(y) (y-k):(y-1)),
      x,
      apply(cbind(x), 1, function(y) (y+1):(y+k))
      )))

    if (clip)
      out = out[out>0]
    return(out)
  }


###################################
#' @name ppng
#' @title ppng
#'
#' @description
#' sends quick plot to ~public_html/plot.png.  If PPNG.DIR env variable defined, then will send to that directory (i.e. instead of public_html)
#' Useful for doing quick standard plots to a static location which one views through tabs in Chrome or other web browser.
#'
#' @param expr Plotting expression eg plot(runif(1000), runif(1000))
#' @param filename filename under ~/public_html/ or Sys.getenv('PPNG.DIR') to dump plots to (='plot.png')
#' @param height integer pixel height of plot (=1000)
#' @param width integer pixel width of plot (=1000)
#' @param dim length 2 integer vector, if expr contains multiple plot calls then will output to matrix of plots with specified dim (=NULL)
#' @param cex expansion factor of plot from "default size" either length 1 scalar or length 2 vector specifying height and width expansion (=1)
#' @param title title to add to plot (='')
#' @param cex.title character expansion factor to title)
#' @param ... additional arguments to png()
#' @author Marcin Imielinski
#' @export
###################################
ppng = function(expr, filename = 'plot.png', height = 1000, width = 1000, dim = NULL, cex = 1, title = NULL, cex.pointsize = min(cex), cex.title = 1,  ...)
  {

    if (length(cex) == 1)
      cex = rep(cex, 2)
    height = cex[1]*height
    width = cex[2]*width

    DEFAULT.OUTDIR = Sys.getenv('PPNG.DIR')
    if (nchar(DEFAULT.OUTDIR)==0)
        DEFAULT.OUTDIR = normalizePath('~/public_html/')

    if (!grepl('^[~/]', filename))
        filename = paste(DEFAULT.OUTDIR, filename, sep = '/')

    if (!file.exists(dirname(filename)))
        system(paste('mkdir -p', dirname(filename)))

    cat('rendering to', filename, '\n')
    png(filename, height = height, width = width, pointsize = 24*cex.pointsize, ...) ## R default for pointsize is 12...

    if (!is.null(dim))
    {
            if (length(dim)==1)
                dim = rep(dim, 2)
            dim = dim[1:2]
            layout(matrix(1:prod(dim), nrow = dim[1], ncol = dim[2], byrow = TRUE))
        }

    eval(expr)
    if (!is.null(title))
        title(title, cex.main = cex.title*max(cex))
    dev.off()
  }


###################################
#' @name ppdf
#' @title ppdf
#'
#' @description
#' sends quick plot to ~public_html/plot.pdf.  If PPDF.DIR env variable defined, then will send to that directory (i.e. instead of public_html)
#' Useful for doing quick standard plots to a static location which one views through tabs in Chrome or other web browser.
#'
#' @param expr Plotting expression eg plot(runif(1000), runif(1000))
#' @param filename filename under ~/public_html/ or Sys.getenv('PPDF.DIR') to dump plots to (='plot.pdf')
#' @param height integer pixel height of plot (=1000)
#' @param width integer pixel width of plot (=1000)
#' @param dim length 2 integer vector, if expr contains multiple plot calls then will output to matrix of plots with specified dim (=NULL)
#' @param cex expansion factor of plot from "default size" either length 1 scalar or length 2 vector specifying height and width expansion (=1)
#' @param title title to add to plot (='')
#' @param cex.title character expansion factor to title)
#' @param ... additional arguments to pdf()
#' @author Marcin Imielinski
#' @export
###################################
ppdf = function(expr, filename = 'plot.pdf', height = 10, width = 10, cex = 1, title = NULL, byrow = TRUE, dim = NULL, cex.title = 1, ...)
  {
    if (length(cex) == 1)
      cex = rep(cex, 2)
    height = cex[1]*height
    width = cex[2]*width


    DEFAULT.OUTDIR = Sys.getenv('PPDF.DIR')
    if (nchar(DEFAULT.OUTDIR)==0)
        DEFAULT.OUTDIR = normalizePath('~/public_html/')

    if (!grepl('^[~/]', filename))
        filename = paste(DEFAULT.OUTDIR, filename, sep = '/')

    if (!file.exists(dirname(filename)))
        system(paste('mkdir -p', dirname(filename)))

    cat('rendering to', filename, '\n')
    pdf(filename, height = height, width = width, ...)

    if (!is.null(dim))
        {
            if (length(dim)==1)
                dim = rep(dim, 2)
            dim = dim[1:2]
            layout(matrix(1:prod(dim), nrow = dim[1], ncol = dim[2], byrow = byrow))
        }

    eval(expr)

    if (!is.null(title))
        title(title, cex.main = cex.title*max(cex))
    dev.off()
  }



#' @name wij
#' @title wij
#'
#' @description
#'
#' Evaluates output of htmlwidget generating expression (e.g. via highcharter) and send to filename in predefined WIDGET.DIR
#' by default plot.html
#'
#' @export
wij = function(expr, filename = 'plot.html', zoom = NULL, cex = 1, force = FALSE, quiet = FALSE, embed = FALSE)
    {
        if (length(cex)==1)
            cex = rep(cex,2)

        if (!force)
            {
                DEFAULT.OUTDIR = Sys.getenv('WIDGET.DIR')
                if (nchar(DEFAULT.OUTDIR)==0)
                    DEFAULT.OUTDIR = normalizePath('~/public_html/')

                if (!grepl('^[~/]', filename))
                    filename = paste(DEFAULT.OUTDIR, filename, sep = '/')
            }

        if (nchar(dirname(filename))==0)
          filename = paste0('./', filename)

        if (!file.exists(dirname(filename)))
            system(paste('mkdir -p', dirname(filename)))

        filename = paste(normalizePath(dirname(filename)), basename(filename), sep = '/')

        widg = eval(expr)

        toWidget <- function(x) {
            htmlwidgets::createWidget(
                name = "plotly",
                x = plotly_build(x),
                width = cex[1]*1000,
                height = cex[2]*1000,
                htmlwidgets::sizingPolicy(
                    padding = 5,
                    browser.fill = TRUE
                    )
                )
        }

        if (!inherits(widg, 'htmlwidget'))
            if (inherits(widg, 'plotly_hash'))
                widg = toWidget(widg)
            else
                stop('expression does not produce valid htmlwidget object')

        if (!is.null(zoom))
            {
                if (is.logical(zoom))
                    zoom = 'x'
                widg = widg %>% hc_chart(zoomType = zoom)
            }

        if (embed)
            return(widg)

        if (quiet == FALSE)
            message('rendering to ', filename)
        htmlwidgets::saveWidget(widg, paste(filename), selfcontained = FALSE)
    }



#' @name wijj
#' @title wijj
#'
#' @description
#'
#' Embeds widget in jupyter notebook
#'
#' @export
wijj = function (x, width = NULL, height = NULL, file = paste0("plotlyJupyterHTML/",
    digest::digest(Sys.time()), ".html"))
{
    if (system.file(package = "IRdisplay") == "") {
        warning("You need the IRdisplay package to use this function: \n",
            "devtools::install_github(c('IRkernel/repr', 'IRKernel/IRdisplay'))")
        return(x)
    }
    l <- plotly_build(x)
    src <- if (is.null(l$url)) {
        dir <- dirname(file)
        if (!dir.exists(dir))
            dir.create(dir, recursive = TRUE)
        owd <- setwd(dir)
        on.exit(setwd(owd), add = TRUE)
        htmlwidgets::saveWidget(as.widget(l), file = basename(file), selfcontained = FALSE)
        file
    }
    else {
        paste0(l$url, ".embed")
    }

    .fun = function (x, y)
    {
        if (length(x) > 0 || is_blank(x))
            x
        else y
    }

    is_blank = function (x)
    {
        inherits(x, "element_blank") && inherits(x, "element")
    }
   #     iframe <- plotly_iframe(src, width %||% l$width, height %||%   l$height)
    iframe <- plotly:::plotly_iframe(src, .fun(width, l$width), .fun(height, l$height))

    get("display_html", envir = asNamespace("IRdisplay"))(iframe)
}



#' @name sortable
#' @title sortable
#' @description
#'
#' dumps sortable list for manual sorting
#' into list.html (in public_html by default)
#'
#' @export
sortable = function(x, filename = 'list.html', title = NULL)
    {

        DEFAULT.OUTDIR = Sys.getenv('WIDGET.DIR')
        if (nchar(DEFAULT.OUTDIR)==0)
            DEFAULT.OUTDIR = normalizePath('~/public_html/')

        if (!grepl('^[~/]', filename))
            filename = paste(DEFAULT.OUTDIR, filename, sep = '/')

        if (!file.exists(dirname(filename)))
            system(paste('mkdir -p', dirname(filename)))

        cat('dropping list to', filename, '\n')
        head1 = '<!doctype html>
  <html lang="en">
  <head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>jQuery UI Sortable - Default functionality</title>
  <link rel="stylesheet"
	href="//code.jquery.com/ui/1.12.0/themes/base/jquery-ui.css">
  <link rel="stylesheet" href="/resources/demos/style.css">
  <style>
  #sortable { list-style-type: none; margin: 0; padding: 0; width:
  60%; }
  #sortable li { margin: 0 3px 3px 3px; padding: 0.4em; padding-left:
  1.5em; font-size: 1.4em; height: 18px; }
  #sortable li span { position: absolute; margin-left: -1.3em; }
  </style>
  <script src="https://code.jquery.com/jquery-1.12.4.js"></script>
  <script src="https://code.jquery.com/ui/1.12.0/jquery-ui.js"></script>
  <script>
  $( function() {
    $( "#sortable" ).sortable();
    $( "#sortable" ).disableSelection();
  } );
  </script>
</head>
<body>'

        head2 = c('<h1>', title, '</h1>', '<ul id="sortable">')

            if (length(x)>0)
          middle = sapply(x, function(y) sprintf('  <li class="ui-state-default"><span class="ui-icon
  ui-icon-arrowthick-2-n-s"></span>%s</li>', y))
      else
          middle = c()

     tail = '</ul>
</body>
</html>
'
     writeLines(c(head1, head2, middle, tail), filename)
    }


###################################
#' @name plop
#' @title plop
#'
#' @description
#' grabs file and plops into public_html (or Sys.getenv('PLOP.DIR') if defined)
#'
#' prefix will be added to left of file name (can include firectory)
#'
#' if fn is list then prefix is expanded to unlisted fn
#'
#' Useful for inspecting a specific subset of analysis files eg when debugging.
#'
#' @param fn character vector of filenames to "plop" into ~public_html/prefix (Sys.getenv('PLOP.DIR') is used as alternative if defined)
#' @param prefix character prefix to add to filenames after plopping (can include subdirectories which can be inspected)
#' @author Marcin Imielinski
#' @export
###################################
plop = function(fn, prefix = NULL, force = NULL)
  {
      if (is.list(fn))
          if (!is.data.frame(fn))
              if (all(sapply(fn, class)=='character'))
                  {
                      l = sapply(fn, length)
                      if (!any(l>0))
                          return(NULL)
                      prefix = rep(prefix[l>0], l[l>0])
                      fn = unlist(fn)
                  }

      if (!is.character(fn))
          {
              new.fn = paste('~/public_html/', prefix, gsub('\\W+', '_', deparse(substitute(fn)), perl = TRUE), '.rds', sep = '')
              if (!file.exists(dirname(new.fn)))
                  system(paste('mkdir -p', dirname(new.fn)))

              saveRDS(fn, new.fn)

              return(new.fn)
          }

      new.fn = paste(prefix, basename(fn), sep = '')

      DEFAULT.OUTDIR = Sys.getenv('PLOP.DIR')
      if (nchar(DEFAULT.OUTDIR)==0)
          DEFAULT.OUTDIR = normalizePath('~/public_html/')

      if (!file.exists(DEFAULT.OUTDIR))
          system(paste('mkdir -p', dirname(filename)))

      if (!is.null(force))
          {
              if (length(force)==1)
                  force = rep(force, length(new.fn))
              new.fn = force
          }

      new.fn = paste(DEFAULT.OUTDIR, new.fn, sep = '/')

      if (any(ix <- !file.exists(dirname(new.fn))))
          sapply(new.fn[ix], function(x)
              system(paste('mkdir -p', dirname(x))))

      mapply(function(x, y) system(paste('cp', x, y)), fn, new.fn)


      return(new.fn)
  }



###################################
#' @name splot
#' @title splot
#' @description
#' convenient formatted scatter plot with additional features as defaults, useful for fast interactive data inspection / exploration of large
#' datasets (eg 1000s of points):
#' - autoamtic setting of solid dots (pch = 19)
#' - transparent colors for over plotting
#' - automatic setting of x and y limits parametrized by "p.outlier"
#' - quick setting of jiggle / jitter on plot
#' - automatic fitting and plotting of regression line (fit = FALSE)
#'
#' @param x numeric vector of x data
#' @param y numeric vector y data
#' @param cex character inspection
#' @param poutlier numeric value between 0 and 1 specifying quantile threshold of  outliers to remove (=0.01)
#' @param col character vector color (=alpha('black', 0.3))
#' @param xlim length 2 numeric vector specifying x limits (=quantile(x, na.rm = T, prob = c(poutlier[1], 1-poutlier[length(poutlier)])))
#' @param ylim length 2 numeric vector specifying y limits (=quantile(y, na.rm = T, prob = c(poutlier[1], 1-poutlier[length(poutlier)])))
#' @param log standard plot log string
#' @param jiggle numeric value between 0 and 1 specifying what percentage of plot area to jiggle each point (useful for overplotting) (= NULL)
#' @param fit logical flag whether to fit a linear regression line to the data (=FALSE)
#' @param col.fit character specifying color of linear regression fit (='blue')
#' @param cex.fit character specifying size of text associated with linear regerssion line (=1)
#' @param square logical flag whether to make square plot
#' @param pch pch
#' @param ... ...
#' @author Marcin Imielinski
#' @export
###################################
splot = function(x, y, cex = 0.4, poutlier = 0.01, col = alpha('black', 0.3),
                 intercept = TRUE, 
    xlim = quantile(x, na.rm = T, prob = c(poutlier[1], 1-poutlier[length(poutlier)])),
    ylim = quantile(y, na.rm = T, prob = c(poutlier[1], 1-poutlier[length(poutlier)])),
    label = NULL,
    cex.label = 1,
    adj.label = c(1, 0.5),
    col.label = 'black',
    log = '',
    jiggle = NULL, ## number between 0 and 1 as percentage of plot
    fit = FALSE,
    col.fit = 'blue',
    cex.fit = 1,
    square = FALSE,
    pch = 19, ...)
    {

        is.inf = is.infinite(x) | is.infinite(y)
        is.inf[is.na(x) | is.na(y)] = TRUE
        if (any(is.inf))
            {
                warning(paste('Removing', sum(is.inf), 'infinite values'))
                x = x[!is.inf]
                y = y[!is.inf]
            }

        if (length(x)==0)
            {
                warning('Empty plot produced')
                plot(0, type = "n")
                return()
            }

        if (grepl('x', log))
            xlim = c(pmax(xlim[1], min(x[which(x>0)])), xlim[2])

        if (grepl('y', log))
            ylim = c(pmax(ylim[1], min(y[which(y>0)])), ylim[2])

        if (square)
            {
                xylim = range(c(xlim, ylim))
                xlim = ylim = xylim;
            }

        if (!is.null(jiggle))
            {
                jiggle = pmax(0, pmin(0.1, jiggle))
                if (grepl('x', log))
                    x = exp(log(x + rnorm(length(x))*jiggle*diff(range(log(xlim)))))
                else
                    x = x + rnorm(length(x))*jiggle*diff(range(xlim))

                if (grepl('y', log))
                    y = exp(log(y + rnorm(length(y))*jiggle*diff(range(log(ylim)))))
                else
                    y = y + rnorm(length(y))*jiggle*diff(range(ylim))
            }

        plot(x, y, cex = cex, col = col, pch = pch, xlim = xlim, ylim = ylim, log =log,  ...)

        if (!is.null(label))
            text(x, y, label, cex = cex.label, adj = adj.label)

        if (fit)
            {
              dat = data.frame(x, y)
              ix = rowSums(is.infinite(as.matrix(dat)), na.rm = TRUE)>0 | rowSums(is.na(dat))
              dat = as.data.table(dat[!ix, ])[ x>=xlim[1] & x<=xlim[2] & y>=ylim[1] & y<=ylim[2], ]

              if (intercept)
              {
                m = lm(y ~ x, dat)
                a = coef(m)[1]
                b = coef(m)[2]
              }
              else
              {
                m = lm(y ~ x-1, dat)
                a = 0
                b = coef(m)[1]
              }

              abline(m, lwd = 3, lty = 2, col = col.fit)
              
              eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                               list(a = format(a, digits = 2),
                                    b = format(b, digits = 2),
                                    r2 = format(summary(m)$r.squared, digits = 3)))
              
              
              if (b>0)
              {
                adj = c(1, 0.5)
                xpos = par('usr')[1] + diff(par('usr')[1:2])*0.8
                ypos = par('usr')[3] + diff(par('usr')[3:4])*0.2
              }
              else
              {
                adj = c(0, 0.5)
                xpos = par('usr')[1] + diff(par('usr')[1:2])*0.2
                ypos = par('usr')[3] + diff(par('usr')[3:4])*0.2
              }
              
              text(xpos, ypos, eq, col = col.fit, cex = cex.fit, adj = adj)
            }
    }



#' @name phist
#' @title phist
#'
#' @description
#' Quick plotlyhistogram
#'
#' @author Marcin Imielinski
#' @import plotly
#' @export
phist = function(expr, data = data.frame(), ...)
    {
        plot_ly(data = data, x = eval(expr), type = 'histogram', ...)
    }


#' @name pscatter
#' @title pscatter
#'
#' @description
#' Quick plotly scatterplot
#'
#' @author Marcin Imielinski
#' @import plotly
#' @export
pscatter = function(x, y, text = '', color = NULL, size = NULL, mode = 'markers', type = 'scatter')
{
  if (!is.character(text))
  {
    text = as.matrix(text)
    if (is.null(colnames(text)))
      colnames(text) = paste0('V', 1:ncol(text))
    text = apply(text, 1, function(x) paste(colnames(text), '=', x, collapse = ', '))    
  }

  data = data.table(x = x, y = y, text = text)

  if (!is.null(color))
    if (is(color, 'character'))
      data$color = color

  if (!is.null(size))
    if (is(color, 'character') | is(color, 'numeric') | is(color, 'integer'))
      data$size = size

  t = paste0("plot_ly(data = data, mode = mode, type = type, ",
            paste0(names(data), '=~', names(data), collapse = ', '), ')')
  eval(parse(text = t))
}


#' @name bubblemap
#' @title bubblemap
#'
#' @description
#' Quick bubble heatmap from matrix mat, use cex to tweak bubble size. 
#' @export
bubblemap = function(mat, col = 'darkgreen', cluster = TRUE, cex = 1, cex.text = 1, zlim = cex*c(0, 10), col.text = 'white', show.legend = FALSE)
{
  if (is.null(rownames(mat)))
    rownames(mat) = 1:nrow(mat)

  if (is.null(colnames(mat)))
    colnames(mat) = 1:ncol(mat)

  rowind = 1:nrow(mat)
  colind = 1:ncol(mat)
  if (cluster)
    {
      rowind = hclust(dist(mat))$order
      colind = hclust(dist(t(mat)))$order
    }
  uv1 = rownames(mat)[rowind]
  uv2 = colnames(mat)[colind]

  res = as.data.table(melt(mat))
  res[, Var1 := factor(Var1, uv1)]
  res[, Var2 := factor(Var2, uv2)]

  gg = ggplot(res, aes(Var1, Var2)) +
    geom_point(aes(size = cex*value), alpha=0.8, color=col, show.legend=show.legend) +
    geom_text(aes(label = signif(value,2), size = cex*cex.text*value/10), color=col.text) +
    #geom_text(aes(label = signif(value,2), size = cex*cex.text/10), color=col.text) +
    scale_size(range = zlim) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.background = element_rect(fill = 'white', colour = 'white')
          )
  print(gg)
}

#' @name bplot
#' @title bplot
#' @description
#'
#' Takes vector y of binary / categorical data (eg genotype) plots a bar plot of the
#' fraction +/- confidence intervals (from prop.test) of each unique value in y. 
#' 
#' Optional argument by and facet provide additonal grouping / faceting variables that will
#' result in the fractions being computed within those groups and the results plotted vertically
#' (for by = ) and horizontally (for facet = ). 
#' 
#'
#' @param y categorical vector across N samples whose fraction will be plotted on the y axis 
#' @param by optional grouping vector, same length as y, to group on prior to plotting fractions
#' @param facet optional second grouping vector, same length as y, to group on
#' @param col single color or named vector to plot or NA (then brewer is used)
#' @param wes if col = NA wes = TRUE will cause wesanderson colors to be used, otherwise standard brewer is used (TRUE)
#' @param keep.base whether to keep the base level when plotting, will be FALSE unless y is logical (ie suppress plotting FALSE value)
#' @param xlab x label
#' @param ylab y label
#' @param counts print counts alongside the by / facet labels (TRUE)
#' @param title title
#' @param return.stats return stats with confidence intervals and proportions in each tranche
#' @return gg plot object or stats
#' @author Marcin Imielinski
#' @export 
bplot = function(y, by = NULL, facet = NULL, col = NA, ylim = NA, keep.base = NULL, xlab = '', counts = TRUE, wes = TRUE, ylab = '', title = '', print = TRUE, return.stats = FALSE)
{

  if (is.null(keep.base))
    keep.base = !is.logical(y)

  if (!inherits(y, 'factor')) ## factor
  {
    if (is.character(y))
      y = factor(y, names(rev(sort(table(y)))))
    else
      y = factor(y)
  }

  dat = data.table(y = y, by = by, facet = facet)

  if (is.null(dat$by))
    dat$by = factor('dummy')
  else if (counts)
    dat$by = paste0('(', dat[, .N, keyby = by][.(dat$by), N], ') ', dat$by)

  if (is.null(dat$facet))
    dat$facet = factor('dummy')
  else if (counts)
    dat$facet = paste0(dat$facet, ' (', dat[, .N, keyby = facet][.(dat$facet), N], ')')
    
  base = levels(y)[1]

  if (is.na(col))
    col = skitools::brewer.master(levels(y), wes = wes)
  else if (is.null(names(col))) ## colors provided as named or unnamed vector     
    col = data.table(uy = levels(y), col = col)[, structure(col, names = uy)] ## replicate if necessary

  dat[, tot := .N, by = .(by, facet)] ## tally totals in each facet

  if (!keep.base)
    dat = dat[y != base, ]

  stats = dat[, .(ndom = .N, frac =.N/tot[1], tot = tot[1]), by = .(y, by, facet)]
  prop.stats = stats[, prop.test(frac*tot, tot) %>% dflm, by = .(y, by, facet)][, .(y, by, facet, frac.lower = ci.lower, frac.upper = ci.upper)] 
  stats = stats %>% merge(prop.stats, by = c("y", "by", "facet"), allow.cartesian = TRUE)

  if (!inherits(stats$by, 'factor'))
    stats$by = factor(stats$by, stats[, sum(frac), by = by][rev(order(V1)), by])
  
  if (!inherits(stats$facet, 'factor'))
    stats$facet = factor(stats$facet, stats[, sum(frac), by = facet][rev(order(V1)), facet])

  if (is.null(by))
    stats[, x := y]
  else
  {
    stats[, x := by]
    stats[, facet1 := y]
  }

  if (!is.null(facet))
    stats[, facet2 := facet]

  if (is.na(ylim))
    ylim = c(0, pmin(1, pmax(max(stats$frac.upper*1.1))))

  favtheme = theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())

  g = ggplot(stats, aes(x = x, y = frac, fill = y)) +
    geom_bar(position = position_dodge(0.7), stat = 'identity') +
    geom_errorbar(aes(ymin = frac.lower, ymax = frac.upper, width = 0.3), position = position_dodge(0.7)) + 
    scale_fill_manual(values = col) + 
    favtheme + 
    ylim(ylim) + ylab(ylab) +
    xlab(xlab) +
    ggtitle(title)

  if (!is.null(stats$facet1) & !is.null(by))
    {
      if (!is.null(stats$facet2))
        g = g + facet_grid(facet1 ~ facet2, scales = 'fixed')
      else
        g = g + facet_grid(facet1 ~ ., scales = 'fixed')
    }
  else
  {
       if (!is.null(stats$facet2))
        g = g + facet_grid(. ~ facet2, scales = 'fixed')
  }

  if (!print)
    return(g)

  if (return.stats)
    return(stats)

  print(g)
}



#' @name vplot
#' @title vplot
#'
#' @description
#' Quick violin plot
#'
#' @param y numeric data vector
#' @param group length(y) vector of categories
#' @param facet1 optional length(y) vector of row categories to facet on (=NULL)
#' @param facet2 optional length(y) vector of column categories to facet on (=NULL)
#' @param transpose logical vector whether flip row / column orientation of facets (=FALSE)
#' @param mapping mapping of groups to colors (=NULL)
#' @param scale scale parameter to geom_vplot (=width), can also be "area" and "count"
#' @param log logical flag whether to log transform (=FALSE)
#' @param count logical flag whether to include counts in ylabels (=TRUE)
#' @param xlab character xlabel (=NULL)
#' @param ylab character ylabel (=NULL)
#' @param reorder logical flag whether to reorder groups by order.fun (which by default is mean)
#' @param reorder.fun function on which to order from left to right (default mean)
#' @param minsup minimum support to include in a group (=NA)
#' @param scatter logical flag whether to include scatter of points (=FALSE)
#' @param alpha numeric vector between 0 and 1 to specify alpha transparency of points if scatter is TRUE (0.3)
#' @param title character specifying plot title (=NULL)
#' @param facet_scales character specifying scales arg in ggplot2::facet_grid()
#' @author Marcin Imielinski
#' @import ggplot2
#' @export
vplot = function(y, group = 'x', facet1 = NULL, facet2 = NULL, transpose = FALSE, flip = FALSE,  mapping = NULL,
                 stat = "ydensity",
    position = "dodge",
    trim = TRUE, sample = NA, scale = "width", log = FALSE, log1p = FALSE, count = TRUE, xlab = NULL, ylim = NULL, ylab = NULL, minsup = NA,
    scatter = FALSE,
    wes = 'Royal1',
    col = NULL, 
    method = 'count',
    sina = FALSE,
    sina.scale = FALSE,
    text = NULL,
    reorder = FALSE,
    reorder.fun = mean,
    cex = 1,
    cex.axis = 1,
    cex.title = 1, 
    cex.scatter = cex,
    col.scatter = alpha('black', 0.5),
    alpha = 0.3, title = NULL, legend.ncol = NULL, legend.nrow = NULL, vfilter = TRUE, vplot = TRUE, dot = FALSE, stackratio = 1, binwidth = 0.1, plotly = FALSE, print = TRUE,facet_scales = "fixed")
    {
        # require(ggplot2)
      if (!is.factor(group))
          group = as.factor(group)
      dat = data.table(y = suppressWarnings(as.numeric(y)), group)

      if (reorder)
      {
        newlev = dat[, reorder.fun(y, na.rm = TRUE), by = group][order(V1), group]
        dat[, group := factor(group, levels = newlev)]
      }

      if (!is.na(sample))
        if (sample>0)
        {
          if (sample<1)
            dat = dat[sample(nrow(dat), round(sample*nrow(sample))), ]
          else
            dat = dat[sample(nrow(dat), round(sample)), ]
        }

        if (is.null(facet1))
            {
                facet1 = facet2
                facet2 = NULL
            }

        if (!is.null(facet1))
            if (!is.factor(facet1))
                facet1 = factor(facet1, unique(facet1))


        if (!is.null(facet2))
            if (!is.factor(facet2))
                facet2 = factor(facet2, unique(facet2))

      suppressWarnings(dat[, facet1 := facet1])
      suppressWarnings(dat[, facet2 := facet2])

      dat = dat[rowSums(is.na(dat))==0, ]
      
      ## remove 0 variance groups
      dat$vgroup = paste(dat$group, dat$facet1, dat$facet2)

      vgroup = NULL ## NOTE fix
      good = as.data.table(dat)[, list(var = var(y)), keyby = vgroup][var>0, vgroup]
      dat = dat[, vfilter := dat$vgroup %in% as.character(good)]
      
      if (!is.na(minsup))
      {
        num = NULL ## NOTE fix
        good = as.data.table(dat)[, list(num = length(y)), keyby = vgroup][num>minsup, vgroup]
        dat = dat[(dat$vgroup %in% as.character(good)), ]
      }
      
      if (nrow(dat)==0)
        stop('No groups exist with >0 variance')
      
      if (count)
      {
        tmp = table(dat$group)
        ix = match(levels(dat$group), names(tmp))
        newnames =  paste(names(tmp)[ix], '\n(', tmp[ix], ')', sep = '')
        if (!is.null(col))
          names(col)[match(levels(dat$group), names(col))] = newnames
        levels(dat$group) = newnames
      }

      if (is.null(mapping))
        mapping = aes(fill=group)

      dat = dat[vfilter!=0, ]
      g = ggplot(dat, aes(y = y, x = group)) + theme_bw(base_size = 15*cex.axis) %+replace% theme(plot.title = element_text(size = 11*cex.axis*cex.title), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())

      if (vplot)
        g = g + geom_violin(mapping = mapping, stat = stat, position = position, trim = trim, scale = scale)
      
      scatter = sina | scatter
      if (scatter)
      {
        if (sina)
        {
          g = g + ggforce::geom_sina(scale = sina.scale, method = method, colour = col.scatter, size = 2*cex.scatter, alpha = alpha)
        }
        else if (dot)
          {
            if (is.null(text))
              g = g + geom_dotplot(data = dat, mapping = aes(x = group, y = y, fill = group), binaxis = 'y', position = 'identity', col = NA, alpha = alpha, method = 'dotdensity', dotsize = cex.scatter, stackratio = stackratio, binwidth = binwidth, stackdir = 'center')
            else
              g = g + geom_dotplot(data = dat, mapping = aes(x = group, y = y, fill = group, text = text), binaxis = 'y', position = 'identity', col = NA, alpha = alpha, method = 'dotdensity', dotsize = cex.scatter, stackratio = stackratio, binwidth = binwidth, stackdir = 'center')
          }
          else
          {
            if (is.null(text))
            {
              if (is.null(col.scatter))
                g = g + geom_jitter(data = dat, mapping = aes(fill = group), shape = 21, size = cex.scatter, alpha = alpha, position = position_jitter(height = 0))
              else
                g = g + geom_jitter(data = dat, fill = alpha(col.scatter, alpha), shape = 21, position = position_jitter(height = 0))
              
            }
            else
            {
              if (is.null(col.scatter))
                g = g + geom_jitter(data = dat, mapping = aes(fill = group, text = text), shape = 21, size = cex.scatter, alpha = alpha, position = position_jitter(height = 0))
              else
                g = g + geom_jitter(data = dat, mapping = aes(text = text), fill = alpha(col.scatter, alpha), shape = 21, position = position_jitter(height = 0))
            }
          }

        }
      
      

        if (log)
            {
                if (!is.null(ylim))
                    if (length(ylim)!=2)
                        ylim = NULL

                if (is.null(ylim))
                    g = g + scale_y_log10()
                else
                    g = g+ scale_y_log10(limits = ylim)
            }
        else if (log1p)
      {
        if (!is.null(ylim))
          if (length(ylim)!=2)
            ylim = NULL
        
        if (is.null(ylim))
          g = g + scale_y_continuous(trans = "log1p")
        else
          g = g+ scale_y_log10(limits = ylim)
      }
        else
            {
                if (!is.null(ylim))
                    if (length(ylim)==1)
                        g = g+ ylim(ylim[1])
                    else if (length(ylim)==2)
                        g = g+ ylim(ylim[1], ylim[2])
            }

        if (!is.null(xlab))
            g = g+ xlab(xlab)

        if (!is.null(ylab))
            g = g+ ylab(ylab)

        if (!is.null(title))
            g = g + ggtitle(title)

        if (!is.null(legend.ncol))
            g = g + guides(fill = guide_legend(ncol = legend.ncol, byrow = TRUE))

        if (!is.null(legend.nrow))
            g = g + guides(fill = guide_legend(nrow = legend.nrow, byrow = TRUE))


      if (!is.null(col))
      {
        g = g + scale_fill_manual(values = col)
        g = g + scale_color_manual(values = col)
      }
      else if (!is.null(wes))
      {
        g = g + scale_fill_manual(values = skitools::brewer.master(length(unique(dat$group)), wes, wes = TRUE))
#        g = g + scale_fill_manual(values = wesanderson::wes_palette(wes))
      }
     

      if (flip)
        g = g + coord_flip()
      
      if (!is.null(dat$facet1))
      {
        if (!is.null(dat$facet2))
        {
          if (transpose)
            g = g + facet_grid(facet2 ~ facet1, scales = facet_scales)
          else
            g = g + facet_grid(facet1 ~ facet2, scales = facet_scales)
        }
        else
        {
          if (transpose)
            g = g + facet_grid(. ~ facet1, scales = facet_scales)
          else
            g = g + facet_grid(facet1 ~ ., scales = facet_scales)
        }
      }

      if (plotly)
        return(ggplotly(g))
      
      if (print)
        print(g)
      else
        g
      
    }


#' @name clock
#' @title clock
#'
#' times expression
#'
#' @param expr R code to eval while suppressing all errors
#' @author Marcin Imielinski
#' @export
clock = function(expr)
  {
    now = Sys.time()
    eval(expr)
    return(Sys.time()-now)
  }


########
#' @name muffle
#' @title muffle
#'
#' Runs code returning NULL is there is any error
#'
#' @param code R code to eval while suppressing all errors
#' @param ... additional tryCatch arguments
#' @return output of evaluated R code or NULL if error
#' @author Marcin Imielinski
#' @export
#########
muffle = function(code, ...)
    {
        return(tryCatch(code, error = function(e) NULL, ...))
    }


##########
# .cc
#
# grabs data.frame columns whose names match the given regex character vector
##########
.cc = function(df, regex, ignore.case = FALSE)
    {
        gr = NULL
        if (is(df, 'GRanges'))
            {
                gr = df[, c()]
                df = values(df)
            }

        if (data.table::is.data.table(df))
            key = key(df)
        else
            key = NULL

        cols = names(df)
        if (is.character(regex))
            {
                cols = cols[rowSums(do.call(cbind, lapply(regex, grepl, cols, ignore.case = ignore.case)))!=0]
                if (!is.null(key))
                    cols = c(key, cols)

                if (data.table::is.data.table(df))
                    df = df[, cols, with = F]
                else
                    df = df[, cols, drop = F]
            }
        else
            {
                cols = cols[regex]

                if (length(cols)==1)
                    return(df[[cols]])

                if (!is.null(key))
                    cols = c(key, cols)

                if (length(cols)==1)
                    if (data.table::is.data.table(df))
                        df = as.data.frame(df)

                if (data.table::is.data.table(df))
                    df[, cols, with = FALSE]
                else
                    df = df[, cols]
            }

        if (!is.null(gr))
            {
                values(gr) = df
                df = gr
            }

        if (data.table::is.data.table(df))
            data.table::setkeyv(df, key)

        return(df)
    }


#' @name col.slice
#' @title col.slice
#' @aliases %!%
#' @description
#' Hacked operator for subsetting columns of data.frames, DataFrames, data.tables, GRanges using
#' a vector of regexps
#'
#' df %!% c('string.*to.*match', 'another.string.to.match')
#'
#' @param df data.frame
#' @param regex string to match or number in which case that column is returned (same behavior for data.table)
#' @return slices of data.frame matching regex
#' @rdname col.slice
#' @exportMethod %!%
#' @export
#' @author Marcin Imielinski
#' @import GenomicRanges
#' @importFrom data.table data.table fread setnames as.data.table
#' @importFrom S4Vectors DataFrame
setGeneric('%!%', function(df, ...) standardGeneric('%!%'))

setMethod("%!%", signature(df = "GRanges"), function(df, ...) {
    return(.cc(df, ...))
})
setMethod("%!%", signature(df = "data.table"), function(df, ...) {
    return(.cc(df, ...))
})
setMethod("%!%", signature(df = "data.frame"), function(df, ...) {
    return(.cc(df, ...))
})

setMethod("%!%", signature(df = "DataFrame"), function(df, ...) {
    return(.cc(df, ...))
})




#' @name row.slice
#' @title row.slice
#' @aliases %~%
#'
#' @description
#' Hacked operator for subsetting rows of data.frames, DataFrames, data.tables, GRanges
#'
#' df %~% 1:5 retrieves first five rows
#'
#' df %~% NA retrieves 5 random rows
#'
#' @param df data.frame
#' @param regex string to match or number in which case that column is returned (same behavior for data.table)
#' @return slices of data.frame matching regex
#' @rdname row.slice
#' @exportMethod %~%
#' @export
#' @author Marcin Imielinski
#' @importFrom data.table setkey := key
setGeneric('%~%', function(df, ...) standardGeneric('%~%'))
setMethod("%~%", signature(df = "GRanges"), function(df, x = NULL) {
    if (is.null(x))
        x = NA
    if (all(is.na(x)))
        x = sample(1:length(df), min(length(df), 5))
    return(df[x, ])
})
setMethod("%~%", signature(df = "data.table"), function(df, x = NULL) {
    if (is.null(x))
        x = NA
    if (all(is.na(x)))
        x = sample(1:nrow(df), min(nrow(df), 5))
    if (is.character(x))
        x = grep(x, df[[gplots::key(df)]], ignore.case = TRUE)
    return(df[x, ])
})
setMethod("%~%", signature(df = "data.frame"), function(df, x = NULL) {
    if (is.null(x))
        x = NA
    if (all(is.na(x)))
        x = sample(1:nrow(df), min(nrow(df), 5))
    return(df[x, ])
})
setMethod("%~%", signature(df = "DataFrame"), function(df, x = NULL) {
    if (is.null(x))
        x = NA
    if (all(is.na(x)))
        x = sample(1:nrow(df), min(nrow(df), 5))
    return(df[x, ])
})
setMethod("%~%", signature(df = "vector"), function(df, x = NULL) {
    if (is.null(x))
        x = NA
    if (all(is.na(x)))
        x = sample(1:length(df), min(length(df), 5))
    return(df[x])
})


#' @name mtable
#' @title mtable
#' @description
#' tabulates unique rows values for matrix / data frame
#'
#' @param mat input matrix
#' @return unique rows of mat, with additional column $count on the left
#' @author Marcin Imielinski
#' @export
mtable = function(mat)
    {
        val = apply(as.matrix(mat), 1, paste, collapse = ',')
        tab = table(val)
        ix = structure(match(names(tab), val), names = names(tab))
        out = cbind(count=as.numeric(tab), as.matrix(mat)[ix, , drop = FALSE])
        rownames(out) = NULL
        return(out)
    }


#' @name ucount
#' @title ucount
#'
#' @description
#' returns vector of same length as input with number of counts of each value
#' in the whole list
#'
#' @param x vector
#' @return length(x) vector with number of instances of each item in x
#' @author Marcin Imielinski
#' @export
ucount = function(x)
    {
        id = id2 = count = NULL ## NOTE fix
        tmp = data.table(id = 1:length(x), id2 = x)[, count := length(id), keyby = id2]
        setkey(tmp, id)
        return(tmp[list(1:length(x)), count])
    }

#' @name more
#' @title more
#'
#' @description
#' "more" +/- grep vector of files
#'
#' @param x vector of iles
#' @param grep string to grep in files (=NULL)
#' @author Marcin Imielinski
#' @export

more = function(x, grep = NULL)
{
    if (is.null(grep))
        x = paste('more', paste(x, collapse = ' '))
    else
      x = paste0('grep -H "', grep, '" ', paste(x, collapse = ' '), ' | more')
    system(x)
}


#' @name sgrep
#' @title System grep
#'
#' @description
#' Greps for expression and returns list with hits, indexed by filename
#' 
#' @param paths character of file paths
#' @param expr regular expression to match 
#' @return data.table of all non
#' @author Marcin Imielinski
sgrep = function(paths, expr) {
  MAXCHUNK = 500
  if (length(paths)>MAXCHUNK)
  {
    numsplits = ceiling(length(paths)/MAXCHUNK)
    message(numsplits)
    pathl = split(paths, ceiling(runif(length(paths))*numsplits))
    out = do.call('c', lapply(pathl, sgrep, expr))
  } else
  {
    p = readLines(pipe(paste0('xargs grep -H "', expr, '" <<< ', paste(paths, collapse = ' '))))
    if (length(p)==0)
      return(lapply(paths, function(x) c()))
    nms = sapply(strsplit(p,'\\:'), '[', 1)
    out = split(p, nms)
  }

  return(out[paths])
}


#' @name sgrepl
#' @title System grep logical
#'
#' @description
#' Greps for expression and returns logical vector if match exists in file
#' 
#' @param paths character of file paths
#' @param expr regular expression to match 
#' @return data.table of all non
#' @author Marcin Imielinski
sgrepl = function(paths, expr) {
  return(paths %in% names(sgrep(paths, expr)))
}


#' @name tailf
#' @title tailf
#'
#' @description
#' "tail -f" +/- grep vector of files
#'
#' @param x vector of iles
#' @param grep string to grep in files (=NULL)
#' @author Marcin Imielinski
#' @export
tailf = function(x, n = NULL, grep = NULL)
{
    if (is.null(n))
        x = paste('tail -f', paste(x, collapse = ' '))
    else
        x = paste('tail -n', n, paste(x, collapse = ' '))

    if (!is.null(grep))
        x = paste(x, '| grep -H', grep, ' | more')

    system(x)
}


#' @name headf
#' @title headf
#'
#' @description
#' "head" +/- grep vector of files
#'
#' @param x vector of iles
#' @param grep string to grep in files (=NULL)
#' @author Marcin Imielinski
#' @export
headf = function(x, n = 5, grep = NULL)
{

    if (is.null(n))
        x = paste('head -f', paste(x, collapse = ' '))
    else
        x = paste('head -n', n, paste(x, collapse = ' '))

    if (!is.null(grep))
        x = paste(x, '| grep -H', grep, ' | more')

    system(x)
}



##  __       __   ______   ________          __                          __
## |  \     /  \ /      \ |        \        |  \                        |  \
## | $$\   /  $$|  $$$$$$\| $$$$$$$$       _| $$_     ______    ______  | $$  _______
## | $$$\ /  $$$| $$__| $$| $$__          |   $$ \   /      \  /      \ | $$ /       \
## | $$$$\  $$$$| $$    $$| $$  \          \$$$$$$  |  $$$$$$\|  $$$$$$\| $$|  $$$$$$$
## | $$\$$ $$ $$| $$$$$$$$| $$$$$           | $$ __ | $$  | $$| $$  | $$| $$ \$$    \
## | $$ \$$$| $$| $$  | $$| $$              | $$|  \| $$__/ $$| $$__/ $$| $$ _\$$$$$$\
## | $$  \$ | $$| $$  | $$| $$               \$$  $$ \$$    $$ \$$    $$| $$|       $$
##  \$$      \$$ \$$   \$$ \$$                \$$$$   \$$$$$$   \$$$$$$  \$$ \$$$$$$$


##############
#' @name maf_disp
#' @title maf_disp
#'
#' @description
#' Quick display of rows of data.frame holding contents of Oncotated maf file
#'
#' @param maf data.frame with Oncotated maf columns
#' @param flavor character specifying 'functional' or 'validation' flavors, which correspond to special column slices of maf data.frame (= NULL)
#' @param sorted logical flag whether to return output sorted by gene, variant classificaiton, uniprot site, then patient (=FALSE)
#' @param show.pat logical flag whether to show patient (=TRUE)
#' @param extra_col character vector of additional columns to include (=NULL)
#' @param gene character vector of Hugo_Symbol to subset on (=NULL)
#' @param pat character vector of Tumor_Sample_Barcodes to subset on (=NULL)
#' @return character vector or sliced data.frame
#' @author Marcin Imielinski
#' @export
##############
maf_disp = function(maf, flavor = NULL, # if null then outputs a string for each event, if == "functional" then outputs special rows of annotated maf files
  sorted = F, show.pat = TRUE, extra_cols = NULL, gene = NULL, pat = NULL)
  {
    # get patient name (different naming conventions in different maf files)

    if (show.pat)
      if (is.null(maf$patient_name))
        {
          if (!is.null(maf$patient))
            maf$patient_name = maf$patient
          else if (!is.null(maf$name))
            maf$patient_name = maf$name
          else if (!is.null(maf$Tumor_Sample_Barcode))
            maf$patient_name = gsub('\\-Tumor', '', maf$Tumor_Sample_Barcode)
        }

    if (!is.null(gene))
      maf = maf[maf$Hugo_Symbol %in% gene, ]

    if (!is.null(pat))
      maf = maf[maf$patient_name %in% pat, ]

    if (is.null(flavor)) # vanilla
      {
        if (show.pat)
          out = paste(maf$patient_name, maf$Hugo_Symbol, maf$Protein_Change, sep = " ")
        else
          out = paste(maf$Hugo_Symbol, maf$Protein_Change, sep = " ");

        if (!is.null(extra_cols))
          out = paste(out, apply(maf[, extra_cols, drop = FALSE], 1, function(x) paste(x, collapse = " ")))
        return(out)
      }
    else if (!is.na(charmatch(flavor, "functional")))
     {
        maf$coord = paste('chr', maf$Chromosome, ":", maf$Start_position, sep = "");

        if (sorted)
          maf = maf[order(maf$Hugo_Symbol, maf$Variant_Classification, as.numeric(maf$UniProt_AApos), maf$patient), ];

        out = maf[, c('Hugo_Symbol', 'patient_name', 'coord', 'Variant_Classification',
          grep('Protein_Change', names(maf), value = TRUE),
          grep('^i\\_t*.*((count)|(tumor))', names(maf), value = TRUE),
          grep('PPH2*.*((Class)|(Prob))', names(maf), value = TRUE), extra_cols)];
        return(out)
      }
    else if (!is.na(charmatch(flavor, "validation")))
     {
        maf$coord = paste('chr', maf$Chromosome, ":", maf$Start_position, sep = "");

        if (sorted)
          maf = maf[order(maf$Hugo_Symbol, maf$Variant_Classification, as.numeric(maf$UniProt_AApos), maf$patient), ];

        out = maf[, c('Hugo_Symbol', 'patient_name', 'coord', 'Chromosome', 'Start_position', 'End_position',
          grep('Protein_Change', names(maf), value = TRUE),  extra_cols)];
        return(out)
      }
  }

##################
#' @name reorder_maf
#' @title reorder_maf
#'
#' @description
#' Re-orders maf to comply with the TCGA MAF specifications (v2.2), tacking on all "non-standard" columns after the initial 32
#'
#' @param maf data.frame of MAF
#' @param include.other logical flag whether to include non-standard maf columns after the standard ones (=TRUE)
#' @return data.frame representing MAF in standard column order
#' @export
#' @author Marcin Imielinski
##################
reorder_maf = function(maf, include.other = TRUE)
{
  maf.cols = c('Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Chromosome', 'Start_position', 'End_position', 'Strand', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS', 'dbSNP_Val_Status', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele2', 'Tumor_Validation_Allele1', 'Tumor_Validation_Allele2', 'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele2', 'Verification_Status', 'Validation_Status', 'Mutation_Status', 'Sequencing_Phase', 'Sequence_Source', 'Validation_Method', 'Score', 'BAM_file', 'Sequencer')

  if (any(setdiff(maf.cols, names(maf))))
    stop(sprintf('MAF is missing the following required columns:\n%s', paste(names(maf), collapse = ", ")))

  other.cols = setdiff(names(maf), maf.cols)

  if (include.other)
    maf = cbind(maf[, maf.cols], maf[, other.cols])

  return(maf)
}


#################
#' @name maflite
#' @title maflite
#'
#' @description
#' take maf data frame and returns columns corresponding to "maflite" format
#' https://confluence.broadinstitute.org/display/CGA/Onctoator#Oncotator-Mafliteformat
#'
#' @param maf maf.data.frame
#' @return data.frame in maf.lite format
#' @export
#' @author Marcin Imielinski
#################
maflite = function(maf)
  {
    cols = c('NCBI_Build', 'Chromosome', 'Start_position', 'End_position', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode');

    if (inherits(maf, 'GRanges'))
      {
        names(maf) = NULL;
        maf = as.data.frame(maf)
        maf$Start_position = maf$start
        maf$End_position = maf$end
        maf$Chromosome = maf$seqnames
      }

    if (length(setdiff(cols, names(maf)))>0)
      stop(sprintf('maf is missing these cols: %s', paste(setdiff(cols, names(maf)), collapse = ', ')))

    return(maf[, cols])
  }



#################
#################
# Variant Territories
#################
#################

#################
#' @name maf_coding
#' @title maf_coding
#'
#' @description
#' Scans "Variant_Classification" field in maf and outputs TRUE if variant is in protein coding region (includes synonymous)
#'
#' @param maf maf data.frame
#' @return logical vector specifying whether row satisfies the criterion
#' @author Marcin Imielinski
#' @export
##################
maf_coding = function(maf, inclusive=T)
{
	if(inclusive) {
		return(maf$Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Silent", "Nonstop_Mutation", "Translation_Start_Site", "Splice_Site", "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame"))
	} else {
		return(!(maf$Variant_Classification %in% c("3'UTR", "5'UTR", "Intron", "3'Flank", "5'Flank", "IGR")))
	}
}

#################
#' @name maf_exonic
#' @title maf_exonic
#'
#' @description
#' Scans "Variant_Classification" field in maf and outputs TRUE if variant is exonic
#'
#' @param maf maf data.frame
#' @return logical vector specifying whether maf row satisfies the criterion
#' @author Marcin Imielinski
#' @export
#'
##################
maf_exonic = function(maf, inclusive=T)
{
  if(inclusive) {
    return(maf$Variant_Classification %in% c("3'UTR", "5'UTR") | maf_coding(maf, inclusive=T) )
  } else {
    return(!(maf$Variant_Classification %in% c("Intron", "3'Flank", "5'Flank", "IGR")))
  }
}

#################
#' @name maf_genic
#' @title maf_genic
#'
#' @description
#' Scans "Variant_Classification" field in maf and outputs TRUE if variant is genic
#'
#' @param maf maf data.frame
#' @return logical vector specifying whether maf row satisfies the criterion
#' @author Marcin Imielinski
#' @export
##################
maf_genic = function(maf, inclusive=T)
{
  if(inclusive) {
    return(maf$Variant_Classification %in% c("Intron", "3'Flank", "5'Flank") | maf_exonic(maf, inclusive=T) )
  } else {
    return(!(maf$Variant_Classification %in% c("IGR")))
  }
}


#' @name maf_sub
#' @title maf_sub
#'
#' @description
#' Scans "Variant_Classification" field in maf and outputs TRUE if variant is a substitution
#'
#' @param maf maf data.frame
#' @return logical vector specifying whether maf row satisfies the criterion
#' @author Marcin Imielinski
#' @export
maf_sub = function(maf, inclusive=T)
  {
    if(inclusive) {
      return(maf$Variant_Type %in% c("SNP" , 'DNP', 'TNP', 'ONP'))
    } else {
      return(!(maf$Variant_Type %in% c("INS", "DEL", "Consolidated")))
    }
  }

#' @name maf_indel
#' @title maf_indel
#'
#' @description
#' Scans "Variant_Classification" field in maf and outputs TRUE if variant is a indel
#'
#' @param maf maf data.frame
#' @return logical vector specifying whether maf row satisfies the criterion
#' @author Marcin Imielinski
#' @export
maf_indel = function(maf, inclusive=T)
  {
    if(inclusive) {
      return(maf$Variant_Type %in% c("INS" , 'DEL'))
    } else {
      return(!(maf$Variant_Type %in% c("SNP" , 'DNP', 'TNP', 'ONP', "Consolidated")))
    }
  }


#' @name maf_syn
#' @title maf_syn
#'
#' @description
#' Scans "Variant_Classification" field in maf and outputs TRUE if variant is synonymous
#'
#' @param maf maf data.frame
#' @return logical vector specifying whether maf row satisfies the criterion
#' @author Marcin Imielinski
#' @export
maf_syn = function(maf, inclusive=T)
  {
    if (!is.character(maf$Protein_Change))
      maf$Protein_Change = as.character(maf$Protein_Change)

    if(inclusive) {
      return(nchar(maf$Protein_Change) == 0 | is.na(maf$Protein_Change) | maf$Variant_Classification == "Silent")
    } else {
      return(!maf_nonsyn(maf, inclusive=T))
    }
  }



#' @name maf_nonyn
#' @title maf_nonyn
#'
#' @description
#' Scans "Variant_Classification" field in maf and outputs TRUE if variant is non-synonymous
#'
#' @param maf maf data.frame
#' @return logical vector specifying whether maf row satisfies the criterion
#' @author Marcin Imielinski
#' @export
maf_nonsyn = function(maf, inclusive=T)
  {
    if (!is.character(maf$Protein_Change))
      maf$Protein_Change = as.character(maf$Protein_Change)

    sil_prot_change = function(pc) {
		AA = strsplit(gsub("^p\\.", "", pc), "\\d+")
		sapply(1:length(AA), function(i) AA[[i]][1] == AA[[i]][2] )
	}

    if(inclusive) {
      return(maf$Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame") |
             (maf$Variant_Classification %in% c("Translation_Start_Site", "Splice_Site") & nchar(maf$Protein_Change)>0) )
    } else {
      return(!maf_syn(maf, inclusive=T))
    }
  }



#' @name maf_ns
#' @title maf_ns
#' @description
#' Scans "Variant_Classification" field in maf and outputs TRUE if variant is non-synonymous
#'
#' @param maf maf data.frame
#' @return logical vector specifying whether maf row satisfies the criterion
#' @author Marcin Imielinski
#' @export
maf_ns = function(maf, inclusive=T) {
	return(maf_nonsyn(maf, inclusive))
}


#' @name maf_truncating
#' @title maf_truncating
#'
#' @description
#' Scans "Variant_Classification" field in maf and outputs TRUE if variant is truncating
#'
#' @param maf maf data.frame
#' @return logical vector specifying whether maf row satisfies the criterion
#' @author Marcin Imielinski
#' @export
maf_truncating = function(maf, inclusive=T)
  {
    if(inclusive) {
      return(maf$Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Nonstop_Mutation", "De_novo_Start_OutOfFrame", "Splice_Site", "Translation_Start_Site"))
    } else {
      return(!(maf$Variant_Classification %in% c("In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Silent", "3'UTR", "3'Flank", "5'UTR", "5'Flank", "IGR", "Intron", "RNA", "Targeted_Region", "De_novo_Start_InFrame")))
    }
  }


####################
# maf_classify
#
# Refines variant classificaiton, plugs into maf_colors
#
###################


#' @name maf_classify
#' @title maf_classify
#'
#' @description
#' Re-classifies oncotated variants
#'
#' @param maf maf data.frame
#' @return variant categories
#' @author Marcin Imielinski
#' @export
maf_classify = function(maf)
  {
    maf$variant.label = "";
    maf$variant.label[maf$Variant_Classification %in% c('Translation_Start_Site', 'De_novo_Start_OutOfFrame', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Start_Codon_Del', 'Start_Codon_Ins', 'Stop_Codon_DNP')] = 'Truncating'
    maf$variant.label[maf$Variant_Classification %in% c("3'UTR", "5'UTR", "5'Flank", "3'Flank")] = 'UTR'
    maf$variant.label[maf$Variant_Classification %in% c("Intron")] = 'Intronic'
    maf$variant.label[maf$Variant_Classification %in% c("Splice_Site")] = 'Splice_Site'
    maf$variant.label[maf$Variant_Classification %in% c('In_Frame_Del', 'In_Frame_Ins')] = 'In_Frame_Indel'
    maf$variant.label[maf$Variant_Classification %in% c('Silent')] = 'Silent'
    if (!is.null(maf$PPH2_Class))
      {
        maf$variant.label[maf$PPH2_Class == 'deleterious'] = 'PPH2_Deleterious';
        maf$variant.label[maf$PPH2_Class == 'neutral'] = 'PPH2_Neutral';
      }
    maf$variant.label[nchar(maf$variant.label)==0] = "Other_Missense";
    return(maf$variant.label)
  }


##########################
# maf_to_simple
#
# Dumps maf files to "simple" format for input into oncotator, adds dummy ref_allele and tum_allele1 cols if not provided
#
##########################


#' @name maf_to_simple
#' @title maf_to_simple
#'
#' @description
#' Dumps maf files to "simple" format for input into oncotator, adds dummy ref_allele and tum_allele1 cols if not provided
#'
#' @param maf maf data.frame
#' @param filename output file
#' @param genome An BSgenome object (was "build genome build (='hg19')")
#' @author Marcin Imielinski
#' @export
maf_to_simple = function(maf, filename, genome)
  {
    if (is.null(maf$ref_allele))
      {
        maf$ref_allele =  Biostrings::getSeq(genome, paste('chr', gsub('chr', '', maf$chr),  sep = ""), maf$start, maf$end)
      }
    if (is.null(maf$tum_allele1))
      maf$tum_allele1 = '-';

     write.table(maf[, c('chr', 'start', 'end', 'ref_allele', 'tum_allele1')], filename, quote = F, sep = "\t", row.names = F)
  }

#' @name maf2vcf
#' @title maf2vcf
#'
#' @description
#' Dumps maf data frame to VCF file "fn"
#'
#' @param maf maf data.frame
#' @param fn output file
#' @author Marcin Imielinski
#' @export
maf2vcf = function(maf, fn)
  {
    ## Prep
    ##

    # currently keeping INS, DEL, SNP, and variants with different end positions separate .. (technically should collapse)
    # however given the current maf setup that will require playing w the reference allele and padding certain bases etc.
    # so skipping for now
    maf$variant_id = paste(maf$Chromosome, maf$Start_position, maf$End_position, maf$Variant_Type, sep = "_");


    # identify unique alt alleles and ref allele
    alt.alleles = aggregate(formula = Tumor_Seq_Allele1 ~ variant_id, data = maf, function(x) paste(unique(x), collapse = ","))

    # code alt alleles
    allele.lists = strsplit(alt.alleles[,2], ',');
    names(allele.lists) = alt.alleles[,1];
    multallele.var = names(which(sapply(allele.lists, length)>1));

    maf$allele_code = 1;
    mult_ix = which(maf$variant_id %in% multallele.var)

    if (length(mult_ix)>0)
      maf$allele_code[mult_ix] = sapply(mult_ix, function(x) match(maf$Tumor_Seq_Allele1[x], allele.lists[[maf$variant_id[x]]]))

    if (is.logical(maf$Tumor_Seq_Allele1) | is.logical(maf$Reference_Allele))
      stop('Either tumor or reference allele is encoded in the maf as a logical value!  How can this be? (Check your maf, genius)')

    # make pre vcf tables that will be later catted
    alt_count = xtabs(t_alt_count ~ variant_id + patient, data = maf, sparse = F)
    ref_count = xtabs(t_ref_count ~ variant_id + patient, data = maf, sparse = F)
    alt_gt = xtabs(allele_code  ~ variant_id + patient, data = maf, sparse = F)

    vcf.table = as.data.frame(matrix(paste(paste('0/', alt_gt, sep = ""), alt_count, ref_count, sep = ':'),
      nrow = nrow(alt_count), ncol = ncol(alt_count), dimnames = dimnames(alt_count)));

    out = maf[!duplicated(maf$variant_id), c('Chromosome', 'Start_position', 'variant_id', 'Reference_Allele'), drop = F]
    names(out)[1:2] = c('chr', 'pos');
    out$alt.allele = alt.alleles[match(out$variant_id, alt.alleles$variant_id), 'Tumor_Seq_Allele1']

    # Dummy cols
    out$Qual = 100;
    out$Filter = "PASS"
    out$Info = ".";
    out$Format = "GT:AC:RC";

    out = cbind(out, vcf.table[out$variant_id, , drop = F]);
    out$chr = as.numeric(gsub('Y', 24, gsub('X', 23, out$chr)))

    #sort
    out = out[order(out$chr, out$pos), ]
#    out$chr = paste('chr', out$chr, sep = "");

    ## Write VCF file
    ##
    f = file(fn, "w");
    vcf.source = "maf2vcf"
    hg.build = maf$NCBI_Build[1]

    # vcf meta data
    writeLines(sprintf('##fileformat=VCFv4.0\n##fileDate=%s\n##source=%s\n##reference=%s', gsub('-', '', Sys.Date()), vcf.source, hg.build), f);
    writeLines('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', f);
    writeLines('##FORMAT=<ID=AC,Number=2,Type=Integer,Description="Alt allele Count">', f);
    writeLines('##FORMAT=<ID=RC,Number=2,Type=Integer,Description="Ref allele Count">', f);

    # vcf header
    writeLines(paste(c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', colnames(vcf.table)), collapse = "\t"), f);
    write.table(out, f, sep = "\t", quote = F, row.names = F, col.names = F);

    close(f);

    return(out);
  }




#########################
# mutpairsd
#
#
##########################

#' @name mutpairsd
#' @title mutpairsd
#'
#' @description
#' Takes maf data.frame and outputs table that lists how many pairs there are <= distance d in amino acid space
#'
#' @param maf maf data.frame
#' @param d distance threshold in amino acid space
#' @author Marcin Imielinski
#' @return clusters of mutations with how many pairs of variants supporting
#' @export
mutpairsd = function(maf, d = 0)
  {
    maf$UniProt_AApos = suppressWarnings(as.numeric(maf$UniProt_AApos));
    tmp = aggregate(formula = UniProt_AApos ~ Hugo_Symbol,  data = maf, function(x) sapply(d, function(y) length(which(dist(x)<=y))));
    out = cbind(tmp[,1, drop = F], as.data.frame(tmp[,2]));
    names(out)[1+c(1:length(d))] = paste('npairs', d, sep = "");
    return(out);
  }


##########################
# mutclusters
#
#
# if max.cluster = TRUE returns maximum size cluster in gene where either all (method == complete) or at least one (method single)
# mutation pair is within distance d
#
# if max.cluster = F then returns number of clusters of mutations of count greater than k within a distance d per gene
#
# eg d = 0, k = 1, will give the number of unique sites with more than 1 perfectly recurrent mutation per gene
#
# Clustering is by default using single-linkage agglomerative clustering, but any method that is input to hclust can be used
##########################

#' @name mutclusters
#' @title mutclusters
#'
#' @description
#' Returns genes with a degree of mutaiton clustering (e.g. ranked by how many k>k.thresh clusters with d<d.thresh pairwise distance, or
#' the largest cluster with those characteristics)
#'
#' if max.cluster = TRUE returns maximum size cluster in gene where either all (method == complete) or at least one (method single)
#' mutation pair is within distance d
#'
#' if max.cluster = F then returns number of clusters of mutations of count greater than k within a distance d per gene
#'
#' eg d = 0, k = 1, will give the number of unique sites with more than 1 perfectly recurrent mutation per gene
#'
#' Clustering is by default using single-linkage agglomerative clustering, but any method that is input to hclust can be used
#'
#' @param maf maf data.frame
#' @param d max distance threshold in amino acid space
#' @param k minimum number of mutations in returned clusters
#' @param method character specifying "single" or "complete" linkage clustering of mutations
#' @param max.cluster logical flag whether to return the gene with the largest cluster (if TRUE) or the most number of clusters (if FALSE) (=TRUE)
#' @author Marcin Imielinski
#' @return genes ranked by numbers of cluster or max.cluster size
#' @export
mutclusters = function(maf, d = 0, k = 1, method = "single", max.cluster = TRUE)
  {
    CHUNKSIZE = 1000;
    maf$UniProt_AApos = gsub('_\\d+', '', maf$UniProt_AApos); # get start position of indels
    maf$UniProt_AApos = suppressWarnings(as.numeric(maf$UniProt_AApos));
    gene.count = table(maf[!is.na(maf$UniProt_AApos), ]$Hugo_Symbol);
    tocluster = names(gene.count[gene.count>k]); # remove singleton genes
    out.list = list();

    cat('Processing genes ')
    for (i in seq(1, length(tocluster), CHUNKSIZE)) # chunk up the computation to conserve memory
      {
        cat(sprintf('%d/%d ', i-1, length(tocluster)));
        chunk = i:min(length(tocluster), i+CHUNKSIZE-1);
        maf.tmp = maf[maf$Hugo_Symbol %in% tocluster[chunk], ];

        if (max.cluster)
          out.list = c(out.list, list(aggregate(formula = UniProt_AApos ~ Hugo_Symbol,  data = maf.tmp,
            function(x) max(table(cutree(hclust(dist(x), method = method), h = d))))))
        else
          out.list = c(out.list, list(aggregate(formula = UniProt_AApos ~ Hugo_Symbol,  data = maf.tmp,
            function(x) length(which(table(cutree(hclust(dist(x), method = method), h = d))>k)))))
      }
    cat('\n Done\n')
    out = do.call('rbind', out.list);

    if (max.cluster)
      names(out)[2] = 'max.cluster.size'
    else
      names(out)[2] = 'nclusters';

    out = out[out[,2]>0, ]
    out = out[order(-out[,2]), ]

    return(out);
  }


##################################
#' @name gatk_callvariants
#' @title gatk_callvariants
#'
#' @description
#' Outputs variant calling pipeline on a set of whole genome (or whole exome) bams to shell script and prints
#' instructions on how to execute from command line.
#'
#' @param bams bams to run variant calling on
#' @param outdir out directory to output to
#' @param run name to give run (=timestamp)
#' @param memlimit memory limit in GB
#' @param queue LSF queue to run on (='week')
#' @param skipclean logical flag whether to skipclean (=TRUE)
#' @param chunk integer chunk of variants to run each section on (=250e3)
#' @param runtype character run type (='wg')
#' @param dry logical flag whether to do dry run (=FALSE)
#' @param dcov dcov parameter to GATK (=200)
#' @param bsub whether to run on LSF vs local (=FALSE)
#' @author Marcin Imielinski
#' @export
##################################
gatk_callvariants = function(bams, outdir = './', run = gsub('[^\\w+]', '_', as.character(Sys.time()), perl = T), queue = 'hour', memlimit = 2, sep = '/', skipclean = T, chunk = 250e3, runtype = 'wg', dry = FALSE, dcov = 200, bsub = FALSE, vqsr.memlimit = memlimit)
{
  # require(tools)
  system(paste('mkdir -p', outdir))
  outdir = tools::file_path_as_absolute(outdir)
  jobsdir = paste(outdir, 'jobs', sep = '/');
  tmpdir = paste(outdir, 'tmp', '', sep = '/')
  runfile = paste(outdir, 'run', sep = '/')
  cleaner.tmpdir = paste(outdir, 'cleaner.tmp', '', sep = '/')
  bamlist = paste(outdir, paste(run, '.bam.list', sep = ''), sep = '/')
  writeLines(bams, bamlist)
  system(paste('mkdir -p', outdir, jobsdir, tmpdir, cleaner.tmpdir))
  jobreport.path = paste(outdir, 'jobreport.txt', sep = '/')
    logfile = paste(outdir, 'log', sep = '/');

  runtype = intersect(c('wg', paste('chr', 1:22, sep = ''), 'chrX', 'chrY'), runtype)

  if (length(runtype)==0)
    stop('runtype must be either wg or chr1-22, chrX, chrY')

  if (skipclean)
    skipclean = '-skipClean'
  else
    skipclean = ''

  if (!dry)
    dry = '-run'
  else
    dry = ''

  if (bsub)
    bsub = '-bsub'
  else
    bsub = ''

  cmd = sprintf('#!/bin/bash\ntouch %s\ncp %s %s.$RANDOM\njava -jar %s -S %s -gvsg -jobReport %s -runName %s -jobProject %s -jobQueue %s -memLimit %s -resMemLimit %s -tempDir %s -jobSGDir %s -log %s -I %s %s -runType %s -chunk %s -pipeMem %s -callMem %s -cleanerTmpDir %s %s -dcov %s %s -vqsrMem %s' ,
          logfile, logfile, logfile, QUEUE.JAR, WGP.PATH, jobreport.path, run, run, queue, memlimit, memlimit, tmpdir, jobsdir, logfile, bamlist, skipclean, runtype, as.integer(chunk), memlimit, memlimit, cleaner.tmpdir, dry, as.character(dcov), bsub, vqsr.memlimit)

  writeLines(cmd, runfile)
  system(paste('chmod +x', runfile))

  cat(sprintf('\nCreated directory %s and wrote bash run script:\n\n %s \n\nto file %s.  \nTo run Queue pipeline type \"source %s\" in the shell\n\n', outdir, cmd, runfile, runfile))
#  return(cmd)
}

#' @name gatk_oncotate
#' @title gatk_oncotate
#'
#' @description
#' Makes shell script to oncotate variants outputted from GATK UG run in directory "gatk.dir", outputs instructions how
#' to run the script.
#'
#' @param gatk.dir output directory containing GATK UG output
#' @param jname job name to run jobs with
#' @param mem max memory in GB (=3)
#' @param queue queue to run on
#' @author Marcin Imielinski
#' @export
gatk_oncotate = function(gatk.dir, jname = 'gatk.oncotate', mem = 3, queue = 'week', ...)
{
  if (!file.exists(gatk.dir))
    stop('this gatk dir does not exist')

  oncotator.dir = sprintf('%s/Oncotator', normalizePath(gatk.dir))
  system(sprintf('mkdir -p %s', oncotator.dir))
  system(sprintf('cp %s/chrs/*/*recal*vcf %s', gatk.dir, oncotator.dir))

  vcf.paths = dir(oncotator.dir, 'vcf', full.names = TRUE)
  cmd = paste(ONCOTATOR.PATH, '-v', vcf.paths, '--output-format MAF', paste(vcf.paths, '.maf', sep = ''), 'hg19')
  bcmd = bsub_cmd(cmd, jname, mem = mem, queue = queue, ...)

  cmd.path = sprintf('%s/cmd.sh', oncotator.dir)
  bcmd.path = sprintf('%s/bcmd.sh', oncotator.dir)
  writeLines(cmd, cmd.path)
  writeLines(bcmd, bcmd.path)

  cat(sprintf('\nCreated oncotator directory %s within gatk dir %s and put cmd.sh and bcmd.sh for local and LSF run of oncotator.  \nTo run these scripts just type "source %s" or "source %s" in shell \n\n', oncotator.dir, gatk.dir, cmd.path, bcmd.path))
}



#' @name gatk_haplotypecaller
#' @title gatk_haplotypecaller
#'
#' @description
#' calls haplotype caller on a set of input bams and a given set of targets,
#' outputting to target.dir
#'
#' intervals are given as GRanges
#'
#' @param outdir out directory to output to
#' @param bams bams to run variant calling on
#' @param intervals GRanges intervals to run haplotype caller on
#' @param dbsnp dbSNP path (=Sys.getenv('GATK.DBSNP'))
#' @param hg genome fasta location (=Sys.getenv('GATK.FASTA'))
#' @param genome name of genome build (='hg19')
#' @param outroot prefix to give output files (='out')
#' @param stand_call_conf confidence for calls (=30)
#' @param stand_emit_conf confience for emission (=30)
#' @param min_mapq minimum mapping quality (=20)
#' @param other.arg other arguments to give to haplotype caller (='')
#' @param run logical flag whether to run immediately or just return character vector of command (= FALSE)
#' @param verbose logical flag (=TRUE)
#' @param write_bam logical flag whether to write the bam (=FALSE)
#' @param oncotate logical flag whether to oncotate output into MAF files (=TRUE)
#' @author Marcin Imielinski
#' @return character vector of command(s) (only if run = FALSE), otherwise just runs command with system call.
#' @export
#'
gatk_haplotypecaller = function(outdir, bams, intervals = NULL, dbsnp = Sys.getenv('GATK.DBSNP'), hg = Sys.getenv('GATK.FASTA'), genome = 'hg19', outroot = 'out', stand_call_conf = 30, stand_emit_conf = 2, min_mapq = 20, other.args = '', run = F, verbose = T, write_bam = F, oncotate = T)
  {
    system(paste('mkdir -p', outdir))
    outdir = normalizePath(outdir)

    out.file = paste(outdir, '/', outroot, '.vcf', sep = '')

    cmd = paste('java -jar', GATK.JAR, '-T HaplotypeCaller', '--dbsnp', dbsnp, '-R', hg,
      paste('-I', normalizePath(bams), collapse = ' '),
      '-stand_call_conf', stand_call_conf, '-stand_emit_conf', stand_emit_conf, '--min_mapping_quality_score', min_mapq,
      '-o', out.file)

    if (!is.null(intervals))
      if (inherits(intervals, 'GRanges'))
        if (length(intervals)>0)
          {
            interval.file = paste(outdir, '/', outroot, '.intervals', sep = '')
            gr2gatk(intervals, interval.file)
            cmd = paste(cmd, '-L', interval.file)
          }

    if (write_bam)
      {
        out.bam.file = paste(outdir, '/', outroot, '.bam', sep = '')
        cmd = paste(cmd, '-bamOutput', out.bam.file)
      }

    if (oncotate)
      {
        maf.file = paste(outdir, '/', outroot, '.maf', sep = '')
        cmd = c(cmd, paste(ONCOTATOR.PATH, '-v  --output-format=MAF', out.file, maf.file, genome))
      }

    cmd.file = paste(outdir, 'run.sh', sep = '/')
    writeLines(c('#!/bin/bash', cmd), cmd.file)

    cat('Created GATK haplotype caller run file', cmd.file, 'which can be run from command line\n')

    if (run)
      system(cmd)
    else
      return(cmd)
}

#' @name pindel
#' @title pindel
#'
#' @description
#' calls pindel on a set of input bams and a given set of targets,
#' outputting to target.dir
#'
#' intervals are given as GRanges
#'
#' @param outdir out directory to output to
#' @param bams vector of input bams
#' @param intervals GRanges of to run on (=NULL)
#' @param isizes integer insert.size to use
#' @param hg genome fasta location (=Sys.getenv('GATK.HG19'))
#' @param genome genome build (='hg19')
#' @param outroot prefix to give output files (='out')
#' @param run logical flag whether to run immediately or just return character vector of command (= FALSE)
#' @param verbose logical flag
#' @param write_bam logical flag whether to write the bam (=FALSE)
#' @param oncotate logical flag whether to oncotate output into MAF files (=TRUE)
#' @param threads number of threads to use
#' @param window.size integer window size to use (=10)
#' @param other.args other args to add (='')
#' @author Marcin Imielinski
#' @return character vector of command(s) (only if run = FALSE), otherwise just runs command with system call.
#' @export
pindel = function(outdir, bams, intervals = NULL, isizes = NULL, hg = Sys.getenv('GATK.FASTA'), genome = 'hg19', outroot = 'out',
  run = F, verbose = T, write_bam = F, oncotate = F, threads = 1, window.size = 10, other.args = '')
  {
    system(paste('mkdir -p', outdir))
    outdir = normalizePath(outdir)
    DEFAULT.INSERT.SIZE = 250

    out.prefix =  paste(outdir, '/', outroot, sep = '')
    config.file = paste(outdir, '/', outroot, '.config', sep = '')

    if (is.null(isizes))
      {
        isizes = DEFAULT.INSERT.SIZE
        cat('Applying default insert size of', DEFAULT.INSERT.SIZE, '\n')
      }
    bams = normalizePath(bams)
    if (is.null(names(bams)))
      names(bams) = basename(bams)

    write.tab(data.frame(bams, isizes, names(bams)), config.file, col.names = F)

    flag = F
    if (!is.null(intervals))
       if (inherits(intervals, 'GRanges'))
        if (length(intervals)>0)
          {
            intervals = gUtils::gr.string(intervals, mb = F)
            flag = T
          }

    if (!flag)
      intervals = 'ALL'

    cmd = paste(PINDEL.PATH, '-f', hg, '-i', config.file, '-o', out.prefix, '-c', intervals, '-T', threads, '-w', window.size)

    if (oncotate)
      {
        maf.file = paste(outdir, '/', outroot, '.maf', sep = '')
        cmd = c(cmd, paste(ONCOTATOR.PATH, '-v  --output-format=MAF', out.file, maf.file, genome))
      }

    cmd.file = paste(outdir, 'run.sh', sep = '/')
    writeLines(c('#!/bin/bash', cmd), cmd.file)

    cat('Created pindel run file', cmd.file, 'which can be run from command line')

    if (run)
      system(cmd)
    else
      return(cmd)
}

#' Label discordant read pairs
#'
#' Labels read pairs discordant based on whether they have (1) ++ or -- strand orientation (2) "-" strand read start
#' is not greater than dmin or less than dmin ahead of  "+" strand read on same chromosome
#'
#' @note need to merge with gr.isdisc
#' @name discordant.pairs
#' @export
discordant.pairs = function(pairs, inter.only = F, ## will only include interchromosomal pairs
  dmin = 50, dmax = 500)
{
  pairs.gr = grl.unlist(pairs)
  chr = as.character(seqnames(pairs.gr));
  str = as.character(strand(pairs.gr));
  st = start(pairs.gr);
  en = end(pairs.gr);

  chr.l = split(chr, pairs.gr$grl.ix)
  str.l = split(str, pairs.gr$grl.ix)
  st.l = split(st, pairs.gr$grl.ix)
  en.l = split(en, pairs.gr$grl.ix)

  tmp.out = sapply(chr.l, function(x) length(unique(x)))>1 ## any item w more than one chrom is discordant

  if (any(!tmp.out))
    tmp.out[!tmp.out] = sapply(str.l[!tmp.out], function(x) length(unique(x)))==1 ## any item w only a single strand is discordant

  if (any(!tmp.out))
    {
      d = sapply(en.l[!tmp.out], max) - sapply(st.l[!tmp.out], min) ## any item w only a single strand is discordant
      tmp.out[!tmp.out] = d<dmin | d>dmax
    }

  out = rep(NA, length(pairs))
  out[as.numeric(names(chr.l))] = tmp.out

  return(out)
}


##########################
# match.seg.id
#
# Just like match.seg for segs1, segs2 which have an $id field (in additinon to $chr, $pos1, $pos2)
# and returns a numeric vector of length nrow(segs1) that gives the position of the (first) segment in segs2 that intersects segment i of segs1
#
# Unlike match.seg, match.seg.id will respect individual $id relationships so that segs can only intersect if they
# share an $id and also intersect on the genome.
#
# Also works on segs with nomenclature handles by "standardize_segs" function.
#
##########################
match.seg.id = function(seg1, seg2, verbose = F)
{
  seg1 = standardize_segs(seg1)
  seg2 = standardize_segs(seg2)

  out = rep(NA, nrow(seg1));
  uid = unique(seg1$id);
  for (id in uid)
  {
    if (verbose==T) print(sprintf('Finished patient %s', id))

    ix1 = which(seg1$id == id);
    ix2 = which(seg2$id == id);
    out[ix1] = ix2[match.seg(seg1[ix1, ], seg2[ix2,])]
  }

  return(out)
}


#' @name allele_multiplicity
#' @title allele_multiplicity
#'
#' @description
#' Given individuals file (with either Purity value or a bsp_participant_id or Tumor_scan_name columns)
#' and maf file annotated with absolut copy number (maf$cn.tot) and column $patient for firehose id.
#'
#' You can provide maf without cn.tot annotations latter case will pull absolute segs
#' using bsp_participant_id in ind and annotate, but may take a few minutes to pull.
#' @author Marcin Imielinski
#' @export
allele_multiplicity = function(ind, maf, abs.seg = NULL, verbose = TRUE)
{
  rownames(ind) = ind$Individual_Id;

  if (is.null(ind$Purity))
    {
      if (verbose) print('Fetching purity values ...')
      ind = get_pp_ind(ind)
    }

  if (is.null(maf$cn.tot))
    {
      if (is.null(abs.seg))
        {
          if (verbose) print('Grabbing absolute segs')
          abs.seg = get_seg(ind, absolut = TRUE);
        }
      if (verbose) print('Matching maf to absolute segs')
      seg_ix = match.seg.id(maf, abs.seg, verbose = TRUE);
      maf[, c('cn.low', 'cn.high', 'cn.tot')] = abs.seg[maf$seg_ix, c('abs_cn.low', 'abs_cn.high', 'tot.cn')];
    }

  if (verbose) print('Computing posteriors for allele multiplicity')

  CHUNKSIZE = 1e4; # chunk up computation in case mafs get insanely large
  maf$totreads = maf$t_alt_count + maf$t_ref_count;
  maf$mult = NA;
  maf$mult.c = NA;
  maf$clonal.p = NA;
  maf$subclonal.p = NA;
  maf$homozygous.c = NA;

  for (i in seq(1, nrow(maf), CHUNKSIZE))
    {
      if (verbose == T) print(sprintf('Doing chunk %d', i));

      this.chunk = i:min(nrow(maf), i+CHUNKSIZE-1)

      if (any(!is.na(maf$cn.tot[this.chunk])))  # otherwise no non-NA cn.tot in this chunk of maf
        {
          max.cn = max(maf$cn.tot[this.chunk], na.rm = TRUE)
          purity = ind[maf$patient[this.chunk], 'Purity']
          multmat = array(NA, dim = c(length(this.chunk), max(maf$cn.tot[this.chunk], na.rm = TRUE))); #stores likelihoods for multiplicity
          multmat2 = array(NA, dim = c(length(this.chunk), max(maf$cn.tot[this.chunk], na.rm = TRUE))); #stores likelihoods for multiplicity
          for (mult in 1:max.cn)
            {
              ix1 = which(maf$cn.tot[this.chunk]>=mult)
              multmat[ix1, mult] = dbinom(maf$t_alt_count[this.chunk[ix1]],
                       maf$totreads[this.chunk[ix1]], purity[ix1]*mult/maf$cn.tot[this.chunk[ix1]])
              multmat2[ix1, mult] = pbinom(maf$t_alt_count[this.chunk[ix1]],
                       maf$totreads[this.chunk[ix1]], purity[ix1]*mult/maf$cn.tot[this.chunk[ix1]])
            }

          multmat = multmat/rowSums(multmat, na.rm = TRUE); # normalize rows to obtain posterior probabilities (s.t. uniform prior)
          ix = which(!is.na(multmat[,1]))
          mult.bestguess = apply(multmat[ix,], 1, which.max);
          maf$mult[this.chunk[ix]] = mult.bestguess
          maf$mult.c[this.chunk[ix]] = sapply(1:length(ix), function(x) multmat[ix[x],mult.bestguess[x]]);
          maf$clonal.p[this.chunk[ix]] = sapply(1:length(ix), function(x) multmat2[ix[x],mult.bestguess[x]])
          maf$clonal.p[this.chunk[ix]] = pmin(maf$clonal.p[this.chunk[ix]], 1-maf$clonal.p[this.chunk[ix]])/0.5;
          maf$subclonal.p[this.chunk[ix]] = 1-multmat2[ix,1]; # subclonal is probability of being on left side of mult 1 binomial
          maf$homozygous.c[this.chunk[ix]] = sapply(ix, function(x) multmat2[x,maf$cn.tot[this.chunk[x]]]); # probability of being inside mult cn.tot binomial
          maf$homozygous.c[this.chunk[ix]] = pmin(maf$homozygous.p[this.chunk[ix]], 1-maf$homozygous.p[this.chunk[ix]])/0.5;
        }
    }
  return(maf)
}

##################
#' @name plot_multiplicity
#' @title plot_multiplicity
#'
#' @description
#' Plots allele fractions and colors in an individual  based on multiplicity call stratified by total copy number value
#' "individual" is a list or data frame with field $Individual_Id, $Purity
#'
#' maf is a maf file annotated with multiplicity (ie the output of allele_multiplicity) with additional fields
#' $cn.tot, $mult, $mult.p
#' @author Marcin Imielinski
#' @export
#################
plot_multiplicity = function(individual, maf, plot.reads = F)
{
  ucn = sort(unique(maf$cn.tot[maf$patient == individual$Individual_Id]))
  ucn = ucn[ucn<10]
  layout(t(matrix(1:((length(ucn)+2)*3), nrow = 3, ncol = length(ucn)+2)), widths = c(1, 4, 2))
  col = suppressWarnings(RColorBrewer::brewer.pal(max(k, 3), 'Accent'));
  par(mar=c(0,0,0,0), las=1)
  plot.new()
  plot.new()
  text(0.5, 0.5, individual$Individual_Id, adj = c(0.5, 0.5))
  plot.new()
  for (k in ucn)
    {
      par(mar=c(0,0,0,0), las=1)
      plot.new()
      text(0.5, 0.5, paste('CN', k))
      ix = which(maf$cn.tot == k & maf$patient == individual$Individual_Id);
      par(mar=c(0.5,1,0.5,1), las=1)
      if (plot.reads)
        plot(maf$t_alt_count[ix]/maf$totreads[ix], maf$totreads[ix], col = col[maf$mult[ix]], pch = 19, xlim = c(0,1), axes = FALSE)
      else
        plot(maf$t_alt_count[ix]/maf$totreads[ix], maf$mult.c[ix], col = col[maf$mult[ix]], pch = 19, xlim = c(0,1), ylim = c(0,1), axes = FALSE);

      if (!plot.reads)
        {
          yticks = seq(0, 1,  0.1)
          ylab = 'Probability';
        }
      else
        {
          m = max(maf$totreads[ix])
          yticks = seq(0, m, round(m/4/10)*10)
          ylab = 'Reads';
        }

      axis(2, at = yticks, labels = TRUE, tick = TRUE, las = 2, ylab = ylab)

      xticks = seq(0, 1, 0.2)
      if (k==ucn[length(ucn)])
        axis(1, at = xticks, labels = TRUE, tick = TRUE, las = 2, ylab = ylab)
      else
        axis(1, at = xticks, labels = FALSE, tick = TRUE, las = 2, ylab = ylab)

      if (k>1)
        abline(v = (1:(k-1))/k*individual$Purity, lty = 3)

      abline(v = individual['Purity'], lty = 1)
      plot.new()
      par(mar=c(0,1,0,1), las=1)
      gplots::legend(legend = paste(1:k), ncol = ceiling(k/4), col = col, x = 'left', y = 'center', pch = 19)
    }
  plot.new()
  plot.new()
  text(0.5, 0.5, 'Allele Fraction', adj = c(0.5, 1))
}


##################
#' @name mut_genecluster
#' @title mut_genecluster
#'
#' @description
#' Greedy divisive clustering of genes based on mutation rates along a provided order (eg order of gene expression)
#'
#' Outputs a list of gene clusters (each list a character vector of gene symbols)
#' @author Marcin Imielinski
#' @export
#################
mut_genecluster = function(genes, maf, cov, p.thresh = 0.05, min.cluster.size = 2, bonferonni = TRUE)
  {
    mut = table(maf$Hugo_Symbol)
    genes = intersect(intersect(genes, rownames(mut)), rownames(cov))

    if (length(dim(cov))==3) #collapse along 3rd dimension if provided
      cov = t(apply(cov, 1, rowSums))

    mut = mut[genes];
    cov = cov[genes, ,drop = FALSE]

    # stack will implement the recursion
    stack = data.frame(begin = 1, end = length(genes), mutrate = sum(mut)/sum(cov), p = NA)
    out.stack = NULL;
    while (nrow(stack)>0)
      {
        pop = stack[1,];
        stack = stack[-1, ];
        tmp = data.frame(p = rep(NA, pop$end-pop$begin), eff = rep(NA, pop$end-pop$begin))

        cov1 = cov[pop$begin]
        mut1 = mut[pop$begin]
        cov2 = sum(cov[(pop$begin+1):pop$end]);
        mut2 = sum(mut[(pop$begin+1):pop$end]);
        for (i in 1:(nrow(tmp)-1))
          {
            O = array(c(cov1-mut1, mut1, cov2-mut2, mut2), dim = c(2, 2));
            tmp[i, 'p'] = suppressWarnings(chisq.test(O)$p.value)
            tmp[i, 'eff'] = abs(log((mut1/cov1)/(mut2/cov2)))
            tmp[i, 'mutrate1'] = mut1/cov1
            tmp[i, 'mutrate2'] = mut2/cov2
            cov1 = cov1 + cov[pop$begin + i]
            mut1 = mut1 + mut[pop$begin + i]
            cov2 = cov2 - cov[pop$begin + i]
            mut2 = mut2 - mut[pop$begin + i]
            if (i %% 1000 == 0)
              print(sprintf('Interval %d %d Iteration %d', pop$begin, pop$end, i))
          }

        best = order(tmp$p, -tmp$eff)[1]
        if (bonferonni)
          P.THRESH = p.thresh/nrow(tmp)
        else
          P.THRESH = p.thresh;

        if (tmp$p[best] < P.THRESH & best>min.cluster.size & (pop$end - pop$begin - best>=min.cluster.size)) # check if we pass bonferonni correction and clusters are big enough
            stack = rbind(stack, data.frame(begin = c(pop$begin, pop$begin+best),
              end = c(pop$begin+best-1, pop$end),
              mutrate = c(tmp[best, 'mutrate1'],  tmp[best, 'mutrate2']),
              p = tmp[c(best, best), 'p']))
        else
          out.stack = rbind(out.stack, data.frame(begin = pop$begin, end = pop$end, mutrate = pop$mutrate, p = pop$p))

        if (length(dev.list())<5)
        {
          dev.new()
          layout(matrix(c(1,2), nrow = 2, ncol = 1))
          plot(-log10(tmp$p), main = sprintf("Interval %d %d Best = %d P = %0.2e eff = %0.2e", pop$begin, pop$end, pop$begin+best-1, tmp$p[best], tmp$eff[best]), pch = 19)
          abline(v = best, lty = 3, lwd = 2, col = 'red')
          plot(tmp$eff, pch = 19)
          abline(v = best, lty = 3, lwd = 2, col = 'red')
        }

        print(out.stack)
      }

    out.stack$genes = lapply(1:nrow(out.stack), function(i) genes[out.stack$begin[i]:out.stack$end[i]])

    out.stack = out.stack[order(out.stack$begin), ]
    return(out.stack)
  }


##################
#' @name mutrate_window
#' @title mutrate_window
#'
#' @description
#' Computes mutation rates along k gene "windows" along an ordered list "genes" of genes.
#' @author Marcin Imielinski
#' @export
#################
mutrate_window = function(genes, maf, cov, window = 100)
  {
    genes = intersect(genes,  rownames(cov))

    mut = rep(0, length(genes))
    names(mut) = genes;

    tmp = table(maf$Hugo_Symbol);
    mut[names(tmp)] = as.numeric(tmp);

    if (length(dim(cov))==3) #collapse along 3rd dimension if provided
      cov = t(apply(cov, 1, rowSums))

    mut = mut[genes];
    cov = rowSums(cov[genes, ,drop = FALSE])

    cov.window = sum(cov[1:window-1])
    mut.window = sum(mut[1:window-1])

    if (length(genes)>(window-1))
      {
        out = rep(NA, length(genes)-window+1)
        for (i in 1:(length(genes)-window))
          {
            cov.window = cov.window - cov[i]
            mut.window = mut.window - mut[i]

            cov.window = cov.window + cov[i+window]
            mut.window = mut.window + mut[i+window]

            out[i] = mut.window / cov.window;
          }
      }

    names(out) = genes[1:length(out)];
    return(out)
  }

#################
#' @name quickSig
#' @title  quickSig
#'
#' @description
#' Quick implementation of Mike / Gaddy binomial / poisson model.  Requires category based coverage output (*.per_gene.coverage.txt file)
#' computed during mutsig preprocess step.
#' (run mutsig_preprocess setting the following additional flags:
#'   P.output_per_gene_coverage = true;
#'   P.output_per_gene_mutation_counts = true;
#'   P.simplified_gene_sample_coverage_table = false;
#'   P.simplified_gene_sample_mutation_counts_table = false;
#' )
#'
#' Computes context-category specific mutation rates either across whole cohort or within strata of gene-patient categories.
#'
#' Outputs a significance table with the following columns for each category (if analyze.categories = TRUE) and/or category "tot" which is across all categories
#' o.k = observed mutations of category k
#' e.k = expected mutations of category k given background model
#' eff.k = log(o.k / e.k) for category k
#' p.k = p value of deviation from expectation under poisson model
#' q.k = q value of deviation
#'
#' @author Marcin Imielinski
#' @export
#################
quickSig = function(maf, # this is the maf file made by mutsig preprocess *** needs mutation categories $categ field that maps to dim 3 in cov array
  cov, # 3D array imported prior to running using load_genecov or string corresponding to file path to mutsig_preprocess per_gene.coverage.txt file
  patients = NULL, # either list of patient names to limit analysis to, or a (labeled) integer vector of patient categories within which to stratify the analysis
  genes = NULL, # either list of patient names to limit analysis to, or a (labeled) vector of integer gene categories within which to stratify the analysis
  analyze.categories = F, # will also analyze each subcategory for deviation from expectation
  remove.silent = TRUE, # remove mutations with variant classification "silent"
  limit.cat = NULL, # integer vector corresponding to the number of categories to limit to
  two.tailed = TRUE
  )
  {

    ##########
    # prep maf if in mutsig style
    #########
    if (is.null(maf$Hugo_Symbol))
      maf$Hugo_Symbol = maf$gene

    if (is.null(maf$Variant_Classification))
      maf$Variant_Classification = maf$classification

    ######
    # Argument processing
    ######

    # load cov from file if string provided
    if (is.character(cov))
      cov = load_genecov(cov, collapse.patients = F)

    # parse patient names /  classes
    if (is.null(patients))
      {
        patients = rep(1, ncol(cov))
        names(patients) = colnames(cov)
      }

    if (is.null(names(patients))) # if nameless patient vector provided, treat values of provided vector as patient names and create a single patient stratum
      {
        names(patients) = patients;
        patients[1:length[patients]] = rep(1, length(patients))
      }

    if (length(setdiff(names(patients), colnames(cov)))>0)
      stop(sprintf('Variable "patients" contains %s uncovered patients, ie patients that are not in cov matrix', length(setdiff(names(patients), colnames(cov)))))

    if (is.null(maf$patient_name))
      maf$patient_name = maf$patient;

    if (is.null(maf$categ))
      stop('MAF file does not have "categ" field!!')

    if (remove.silent)
        maf = maf[maf$Variant_Classification != "Silent", ]

    if (is.null(genes))
      genes = rownames(cov)
    else if (is.character(genes))
      {
        thrown.out = setdiff(genes, rownames(cov))
        if (length(thrown.out)>0) warning(sprintf('%s genes in list thrown out', length(thrown.out)))
        genes = intersect(genes, rownames(cov));
      }

    ######
    # Calculate sig
    ######

    cov = cov[,,-dim(cov)[3]]; # last category in cov is "total" category, which we remove
    maf = maf[maf$Hugo_Symbol %in% genes & maf$patient_name %in% names(patients), ];
    tmp = table(maf$Hugo_Symbol, maf$patient_name, maf$categ);
##    muts = cov*0;  # compute muts from maf file
    muts[rownames(tmp), colnames(tmp), as.numeric(dimnames(tmp)[[3]])] = tmp;
    cov = cov[genes, names(patients), ];
    muts = muts[genes, names(patients), ];

    if (!is.null(limit.cat))
      {
        limit.cat = intersect(1:dim(cov)[3], limit.cat);
        if (length(limit.cat)>0)
          {
            cov = cov[,,limit.cat, drop = FALSE];
            muts = muts[,,limit.cat, drop = FALSE];
          }
        else
          stop('Improper categories given for limit.cat!');
      }

    pat.cat.mat = sign(as.matrix(table(1:length(patients), patients))); # binary matrix of patient categories that we can do nifty multiplifications with

    # rates is dimension patient categories x mutation categories
    rates = array(0, dim = c(ncol(pat.cat.mat), dim(cov)[3]))
#    rates = apply(muts, 3, sum)/apply(cov, 3, sum)
    denom = t(apply(cov, 1, colSums))

    ## rates is patient categories x context categories
    rates = sapply(1:dim(muts)[3], function(k) colSums(muts[,,k]%*%pat.cat.mat, na.rm = TRUE) / colSums(cov[,,k]%*%pat.cat.mat, na.rm = TRUE))
    if (is.null(dim(rates))) # unflatten if no rates
        rates = matrix(rates, nrow = 1)
#    denom = apply(cov[,,, drop = FALSE], 3, function(x) x %*% pat.cat.mat);  # sum coverage across patient categories

    # denom correspond to denominator (coverage)
    # in each gene x patient category x mutation category
    denom = array(0, dim = c(nrow(cov), ncol(pat.cat.mat), dim(cov)[3]), dimnames = list(rownames(cov)))
    for (i in 1:dim(cov)[[3]])
      denom[,,i] = cov[,,i]%*%pat.cat.mat;

    # lambda is expected number of mutations of category k in each gene, ie the poisson lambda
    lambda = array(0, dim = dim(denom))

    # output table
    sig = data.frame(gene = genes, stringsAsFactors = F);

    Ps = array(NA, dim = c(nrow(cov), dim(cov)[3]))

    for (k in 1:dim(cov)[3])
      {
        tmp.p = array(NA, dim = c(nrow(cov), ncol(pat.cat.mat)))
        for (j in 1:ncol(pat.cat.mat))
          {
            lambda[,j,k] = denom[,j, k] %*% rates[j, k, drop = FALSE]; # lambda is the patient / context specific expected value

            # observed mut of context category k in patient category j
            obsmut = rowSums(muts[,pat.cat.mat[,j, drop = FALSE]!=0,k, drop = FALSE])

            ## compute p value for deviation of all genes in this category from poisson expectation
            tmp.p[, j] = ppois(obsmut, lambda[,j,k, drop = FALSE], log.p = F, lower.tail = TRUE);
            if (two.tailed)
              {
                tmp.p.2 = ppois(obsmut, lambda[,j,k], log.p = F, lower.tail = FALSE);
                tmp.p[,j] = pmin(tmp.p[,j], tmp.p.2)*2; # two tailed p value
              }
#            tmp.p[lambda[,j,k, drop = FALSE]==0, j]=1;

          }
                                        #            sig[,paste('c.', k, sep = "")] = rowSums(denom[,, k, drop = FALSE]);
        sig[,paste('c.', k, sep = "")] = rowSums(denom[,, k, drop = FALSE])
        sig[,paste('n.', k, sep = "")] = rowSums(muts[,,k,drop = FALSE]);
        sig[,paste('e.', k, sep = "")] = rowSums(lambda[,,k, drop = FALSE]);
        Ps[, k] = sig[,paste('p.', k, sep = "")] = fisher.combined(tmp.p);
        sig[,paste('eff.', k, sep = "")] = log2(sig[,paste('n.', k, sep = "")]/sig[,paste('e.', k, sep = "")]); # eff > 0 = positive selection, eff<0 = negative selection
        sig[,paste('q.', k, sep = "")] = p.adjust(Ps[,k], 'BH');
      }

    p = 1-fisher.combined(Ps); # do fisher combined of category p values to get overall p value

    sig[,'c.tot'] = rowSums(denom);
    sig[,'n.tot'] = rowSums(muts);
    sig[,'e.tot'] = rowSums(lambda);
    sig[,'p.tot'] = p;
    sig[,'eff.tot'] = log(obsmut/lambda); # eff > 0 = positive selection, eff<0 = negative selection

    ix = which(!is.na(p))
    qval = p.adjust(p[ix], 'BH');
    sig[,'q.tot']  =  NA;

    sig = sig[order(sig$p.tot), ];
    rownames(sig) = sig$gene;
    sig = sig[,-1];
    return(sig)
  }

#' @name pmGSEA
#' @title poor mans GSEA
#' @description
#'
#' pmGSEA "poor man's GSEA " ***
#'
#' Given a gene.set (character vector) or gene.sets (list of character vectors)
#' and given a named vector of significance values or table of significant genes (sig.table)
#' (if table then significance column is $p or first column) identifies gene sets that have significant
#' negative deviation of a "signed K-S" statistic vs uniform distribution  (ie have p values significantly
#' clustering towards zero) ie are significantly enriched in genes showing positive selection.
#'
#' if positive.selection = F, will identify sets with significantly positive deviation of a "signed K-S" statistic (ie have p values significantly clustering towards 1)
#' these are sets showig significant negative selection.
#'
#' All p-values are computed against a distribution of signed K-S statistic obtained through permutation using random gene sets of the same size chosen from sig.table
#'
#' Will adaptively perform permutations between minperms and maxperms using following rule of thumb: if there are <PERM.THRESH permutations with
#' greater than (lower.tail = F) or less than (lower.tail = T) score than observed score, then will compute additional perms
#'
#' *** actually not much poorer than the original GSEA, basically a reimplementation of Mootha et al Nat Gen 2002
#'
#' @param gene.sets a named list of character vectors, each list item is a gene set, i.e. a character vector of genes
#' @param sig.table named vector of p values from an analysis e.g. mutSig, the names of the genes are
#' @param min.perms minimum number of permutations to do in the adaptive permutation test
#' @param max.perms maximum number of permutations to do in the adaptive permutation test
#' @param length.range length 2 integer vector specifying min and max gene set size to score after intersection with genes in sig.table default: c(5,50)
#' @export
#' @author Marcin Imielinski
pmGSEA = function(gene.sets, sig.table, min.perms = 1e2, max.perms = 1e5,
  positive.selection = T, # if positive.selection = F will look at genes enriched in high p values (ie negative selection)
  length.filter = F,
  length.range = c(5, 50), # lengths of gene sets to consider if length filtering
  mc.cores = 1,
  plot.hist = F, verbose = F, bootstrap = T, rank.test = F)
{
  PERM.THRESH = 10;
  if (is.character(gene.sets))
      gene.sets = list(gene.sets)

  if (inherits(sig.table, 'data.frame'))
    {
      if (is.null(sig.table$p))
        tmp = sig.table[,1]
      else
        tmp = sig.table$p

      names(tmp) = rownames(sig.table)
      sig.table = tmp;
    }
  sig.table = sig.table[!is.na(sig.table)]

  out = data.frame(name = names(gene.sets), p = rep(NA, length(gene.sets)), fdr = NA, genes = '', stringsAsFactors = F)

  if (length.filter)
    {
      gene.sets = lapply(gene.sets, function(x) intersect(names(sig.table), x))
      out$length = sapply(gene.sets, length)
      ix = out$length>=length.range[1] & out$length<=length.range[2]
      out = out[ix, ]
      gene.sets = gene.sets[ix]
    }

  if (nrow(out)==0)
    return(out)

  out = rbindlist(mclapply(1:nrow(out), function(i)
#  for (i in 1:nrow(out))
  {
    if (verbose) print(sprintf('Starting %s of %s lists', i, nrow(out)))

      gene.set = gene.sets[[i]];
      perm = c();
      gene.set = intersect(gene.set, names(sig.table));

      if (length(gene.set)>0)
        {

                                        # used to calculate K-S statistic
          D_plus = seq(1/length(gene.set), 1, 1/(length(gene.set)));
          D_minus = D_plus- 1/length(gene.set)

          obs.dat = sort(sig.table[gene.set])

          score = c(obs.dat-D_plus, obs.dat-D_minus)
          if (positive.selection)
            {
              obs = max(-pmin(0, score))
              leading.edge = unique(names(score[1:which.max(-score)]))
            }
          else
            {
              obs = max(pmax(0, score))
              leading.edge = unique(names(score[1:which.max(score)]))
            }

          for (this.perm in 10^(seq(log10(min.perms), log10(max.perms))))
            {
              if (verbose) print(sprintf('Finished %s permutations .. ', length(perm)))
              n.perm = this.perm - length(perm)

              if (bootstrap)
                tmp.ind = sample(length(sig.table), length(gene.set)*n.perm, replace = TRUE)
              else
                tmp.ind = unlist(sapply(1:ceiling(length(gene.set) * n.perm / length(sig.table)), function(x) sample(length(sig.table), length(sig.table))))

              data.vec = sig.table[tmp.ind[1:(n.perm*length(gene.set))]]
              data.vec = data.vec[order(rep(1:n.perm, each = length(gene.set)), data.vec)]; # sorts each "permuted gene list")
              data.mat = matrix(data.vec, nrow = n.perm, byrow = T)

              if (positive.selection)
                perm = c(perm, apply(data.mat, 1, function(x) max(-pmin(c(x-D_plus, x-D_minus)))))
              else
                perm = c(perm, apply(data.mat, 1, function(x) max(pmax(c(x-D_plus, x-D_minus)))))

              if (sum(perm>obs) > PERM.THRESH) break()
            }

          out = data.table(
            name = names(gene.sets)[[i]],
            len = length(gene.sets[[i]]),
            leading.edge = paste(paste(names(obs.dat[leading.edge]), ' (', obs.dat[leading.edge], ')', sep = ''), collapse = ', '),
            genes = paste(paste(names(obs.dat), ' (', obs.dat, ')', sep = ''), collapse = ', '),
            p = signif((1+sum(perm>obs))/(1+length(perm)), 2),
            obs = obs)

          ## out[i, 'leading.edge'] = paste(paste(names(obs.dat[leading.edge]), ' (', obs.dat[leading.edge], ')', sep = ''), collapse = ', ')
          ## out[i, 'genes'] = paste(paste(names(obs.dat), ' (', obs.dat, ')', sep = ''), collapse = ', ')
          ## out[i, 'p'] = signif((1+sum(perm>obs))/(1+length(perm)), 2);
          ## out[i, 'obs'] = obs;
                                        #      if (verbose) print(out[i,  setdiff(names(out), 'genes')])
          if (verbose) print(out$name)
          return(out)
        }
    }, mc.cores = mc.cores))

  out$fdr = p.adjust(out$p, 'BH');

  out = out[order(p), ]

  return(out)
}


####################
#' @name gsea_leading_edge
#' @title gsea_leading_edge
#'
#' @description
#' Draws gsea plot for an input gene set and outputs leading edge
#'
#' @author Marcin Imielinski
#' @export
####################
gsea_leading_edge = function(gene.set, sig, draw.plot = T, cex = 1, asp = 2, eps = 1e-16, name = '')
  {
    p = sig[gene.set, pval]+eps;
    names(p) = gene.set
    D_plus = seq(1/length(gene.set), 1, 1/(length(gene.set)));
    D_minus = D_plus- 1/length(gene.set)
    ix = order(rep(1:length(p), 2));
    p = sort(p);
    score = (-c(p-D_minus, p-D_plus)[ix]);
    max.score = max(c(score, 0));
    leading.edge = unique(names(score[1:which.max(score)]))

    if (draw.plot)
      {
        xx = c(p, p)[ix];
        yy = pmax(score, 0);
        names(xx) = names(yy) = NULL
        par(xpd = NA)
        xlim = c(min(p)/2, 1)
        ylim = c(0, max(score)*1.6)
        plot(xx, yy, xlim = xlim, ylim = ylim, xlab = 'P value', ylab = 'K-S Running sum', main = paste(name, '\nScore', round(max.score,2)), log = 'x')
        polygon(c(0, xx, 1, 0), c(0, yy, 0, 0), col = 'gray')
        ix.le = match(leading.edge, names(p))
        score.yy = -pmin(p-D_minus, p-D_plus)[ix.le]+0.02*max(score)
#        text.yy = pmin(0.85*par('usr')[4], rep(par('usr')[4], length(ix.le))*rep(c(0.1, 0.2, 0.3), length(ix.le))[1:length(ix.le)] + score.yy)
        text.xx = 10^seq(log10(xlim[1]), log10(xlim[2]), diff(log10(xlim)/(length(ix.le)-1)))
        text.yy = par('usr')[4]*0.8
        text(text.xx, text.yy, names(p)[ix.le], cex = cex, srt = 90, adj = c(0, 0.5))
        MID.POINT = 0.8;
        segments(p[ix.le], MID.POINT*text.yy, text.xx, -0.02*max(score)+text.yy, lty = 3)
        segments(p[ix.le], score.yy, p[ix.le], MID.POINT*text.yy, lty = 3)
      }

    return(leading.edge)
  }


################
#' @name ccf
#' @title ccf
#'
#' @description
#' computes fuzzy histogram of ccf across a set of n mutation calls
#' using alt allele count, total counts, cn, and purity, and provided
#' grid.size
#'
#' as per landau, carter et al
#'
#' altc, totc, and cn are of length n, purity is length
#' @author Marcin Imielinski
#' @export
################
ccf = function(altc, totc, cn, purity, grid.size = 0.01, verbose = F)
  {
    if (length(altc) != length(totc) | length(totc) != length(cn))
      stop('altc, totc, cn must have same length')

    EPS = 1e-16;

    fc = purity/((1-purity)*2 + purity*cn)
    x.ccf = seq(0, 1, grid.size);
    x.ccf = x.ccf[2:length(x.ccf)]
    p.ccf = rep(0, length(x.ccf))
    ix = !is.na(cn)
    ix[cn==0] = FALSE;
    if (any(ix))
      for (i in which(ix))
        {
          if (verbose == T & (i %% 1000)==0)
            cat(i, '\n')
          tmp = dbinom(altc[i], totc[i], fc[i]*x.ccf)
          if (any(is.na(tmp)))
            browser()
          if (any(is.na(p.ccf)))
            browser()
          p.ccf.last = p.ccf;
          p.ccf = p.ccf + tmp/(sum(tmp)+EPS)
        }
    else
      return(p.ccf)
    names(p.ccf) = x.ccf;
    return(p.ccf/sum(ix))
  }


################
#' @name af
#' @title af
#'
#' @description
#' computes prob density of af over a set of n mutation calls
#' using alt allele count, total counts, cn, and provided
#' grid.size
#'
#' as per landau, carter et al
#'
#' altc, totc, and cn are of length n, purity is length 1
#' @author Marcin Imielinski
#' @export
################
af = function(altc, totc, grid.size = 0.01, verbose = F)
  {
    if (length(altc) != length(totc))
      stop('altc, totc must have same length')

    EPS = 1e-16;

    x.af = seq(0, 1, grid.size);
    x.af = x.af[2:length(x.af)]
    p.af = rep(0, length(x.af))
    if (length(altc)>0)
      for (i in 1:length(altc))
        {
          if (verbose == T & (i %% 1000)==0)
            cat(i, '\n')
          tmp = dbinom(altc[i], totc[i], x.af)
          if (any(is.na(tmp)))
            browser()
          if (any(is.na(p.af)))
            browser()
          p.af.last = p.af;
          p.af = p.af + tmp/(sum(tmp)+EPS)
        }
    else
      return(p.af)

    names(p.af) = x.af;
    return(p.af/length(altc))
  }

################
#' @name af2
#' @title af2
#'
#' @description
#' computes 2D probability density af over a set of n mutation counts
#' two samples, output columns correspond to x and rows correspond to y
#' and dimnames correspond to amounts
#'
#' as per landau, carter et al
#'
#' altc1, totc1, altc2, totc2 are of length n
#'
#' @author Marcin Imielinski
#' @export
################
af2 = function(altc1, totc1, altc2, totc2, grid.size = 0.01, verbose = F, animate = NA)
  {
    if (length(altc1) != length(totc1))
      stop('altc, totc must have same length')

    EPS = 1e-16;

    x.af = seq(0, 1, grid.size)[-1];
    y.af = seq(0, 1, grid.size)[-1];
    p.af = matrix(0, ncol = length(x.af), nrow = length(y.af))
    if (animate)
      image(p.af, add = F, col = terrain.colors(100))
    ix = which(!is.na(altc1) & !is.na(altc2) & !is.na(totc1) & !is.na(totc2))
    if (length(ix)>0)
      for (i in ix)
        {
          if (verbose == T & (i %% 50)==0)
            cat(i, '\n')
          tmp1 = dbinom(altc1[i], totc1[i], x.af)
          tmp2 = dbinom(altc2[i], totc2[i], y.af)
          tmp = cbind(tmp2) %*% rbind(tmp1)
          p.af.last = p.af;
#          p.af = p.af + tmp/(sum(tmp)+EPS)
          p.af = p.af + tmp
          if (!is.na(animate) & (i %% animate)==0)
            image(p.af, add = T, col = terrain.colors(100))
        }
    else
      return(p.af)

    dimnames(p.af) = list(y.af, x.af)
    return(p.af/length(ix))
  }



################
#' @name clone_cluster
#' @title clone_cluster
#'
#' @description
#' determines "clone membership" using CCF threshold for n variants, given alt read count, total read count, and purity
#'
#' first fits k component mixture model (k pre-specified) to CCF histogram .. mixture model can also be given as input
#' returns cluster centers and membership probabilities for each mutation
#'
#' $mu k vector of means
#' $sigma k vector of sigma
#' $p.cluster n x k matrix of cluster probabilities
#' $lambda k vector of cluster membership frequencies
#' return n x k matrix of probabilities that CCF>ccf.thresh for each variant
#'
#' altc, totc, and cn are of length n, purity is length 1
#'
#' @author Marcin Imielinski
#' @export
clone_cluster = function(altc, totc, cn, purity, thresh = 0.95, k = 2, mix.model = NULL, grid.size = 0.01, verbose = F, nsamp = 1e4)
  {
    if (length(altc) != length(totc) | length(totc) != length(cn))
      stop('altc, totc, cn must have same length')

    p.ccf = NULL;

    if (is.null(mix.model))
      {
        cat('Generating CCF histogram\n')
        p.ccf = ccf(altc, totc, cn, purity, grid.size, verbose);
        cat('Clustering histogram\n')
        s.ccf = sample(as.numeric(names(p.ccf)), 1e4, prob = p.ccf, replace = T)
        mix.model = mixtools::normalmixEM(s.ccf, k = k)
      }
    else
      k = length(mix.model$mu)

    cat('Computing cluster membership\n')
    res = list(mu = mix.model$mu, sigma = mix.model$sigma);
    fc = purity/((1-purity)*2 + purity *cn)
    x.ccf = seq(0, 1, grid.size)
    x.ccf = x.ccf[2:length(x.ccf)]
    ix.clonal = x.ccf>=thresh
    p.cluster = matrix(NA, nrow = length(altc), ncol = length(res$mu))
    ix = altc[!is.na(altc)] & !is.na(cn) & !is.na(totc)
    ix[cn==0] = FALSE
    if (any(ix))
      for (i in which(ix))
        {
          if (verbose == T & (i %% 1000)==0)
            cat(i, '\n')
          p.ccf.i = dbinom(altc[i], totc[i], fc[i]*x.ccf)
          p.ccf.i = p.ccf.i/sum(p.ccf.i)
          p.mix.i = t(matrix(dnorm(rep(x.ccf, k), mean = rep(res$mu, each = length(x.ccf)), sd = rep(res$sigma, each = length(x.ccf))), ncol = k))
          p.cluster[i,] = p.mix.i %*% p.ccf.i
        }

    res$p.cluster = p.cluster / rowSums(p.cluster)
    res$lambda = colSums(res$p.cluster, na.rm = T);
    res$lambda  = res$lambda/sum(res$lambda, na.rm = T)
    res$p.ccf = p.ccf;
    return(res)
  }


##  ______   ______   __     __          __                          __
## |      \ /      \ |  \   |  \        |  \                        |  \
##  \$$$$$$|  $$$$$$\| $$   | $$       _| $$_     ______    ______  | $$  _______
##   | $$  | $$ __\$$| $$   | $$      |   $$ \   /      \  /      \ | $$ /       \
##   | $$  | $$|    \ \$$\ /  $$       \$$$$$$  |  $$$$$$\|  $$$$$$\| $$|  $$$$$$$
##   | $$  | $$ \$$$$  \$$\  $$         | $$ __ | $$  | $$| $$  | $$| $$ \$$    \
##  _| $$_ | $$__| $$   \$$ $$          | $$|  \| $$__/ $$| $$__/ $$| $$ _\$$$$$$\
## |   $$ \ \$$    $$    \$$$            \$$  $$ \$$    $$ \$$    $$| $$|       $$
##  \$$$$$$  \$$$$$$      \$              \$$$$   \$$$$$$   \$$$$$$  \$$ \$$$$$$$


##################
#' @name igv
#' @title igv
#'
#' @description
#' Controls IGV on localhost (or specified host, separate from where R session is running).  Igv application must be running and listening to a specified port.  Then if you configure this port
#' via environment variables (IGV_HOST, IGV_PORT) in the current R session then you can apply the following usages
#'
#' igv(fn) ## sends any given file(s) into igv (eg .bam, .wig, .bed)
#' igv(loci = cool.loci) ## plots the windows specified as GRanges or IGV-parsable strings (eg gene name)
#' igv(gr = cool.gr) ## sends granges object to igv session, Note: currently requires the ability to write to a public_html that is web viewable by computer running IGV
#' igv(snapshot = fn) ## sends current screen to file
#' igv(new = TRUE) ## refreshes current session
#' igv(reset = TRUE) ## resets connections, sometimes useful if IGV not responding
#'
#' If alternate file paths are present on server (where R is runinng) and computer running IGV, then can specify gsub.paths variable which is a list of length 2 vectors
#' specifying how to convert file paths from arguments given to igv to ones that can be loaded locally.
#'
#' @param paths file paths to display on current igv session (=NULL)
#' @param loci GRanges or IGV parsable string specifying what window(s) on genome to view (=NULL)
#' @param gr  GRanges or GRAngesList of numeric genomic data or interval genomic annotations to send to IGV session, if gr has field $score then data will be dumped to .bw otherwise to .bed or .gff (=NULL)
#' @param snapshot file path to store snapshot in (has to be interpretable on file system where IGV is running)
#' @param track.view command for setting the track display mode ("expand","squish" or "collapse")
#' @param new logical flag whether to start new IGV session
#' @param reset logical flag whether to reset connection between R and IGV (useful if IGV non responsive)
#' @param host character specifying host where IGV is running
#' @param port integer specifying port where IGV is running
#' @param mac logical flag specifying whether host is a local "mac" (ie then apply gsub.paths) otherwise keep paths the same
#' @param gsub.paths list of length 2 vectors specifying gsub args to apply to filenames when mac = TRUE
#' @export
#' @author Marcin Imielinski
##################
igv = function(
    paths = NULL, ## server file system paths to load
    gr = NULL, ## granges object
    loci = NULL, #can either be a list or data.frame with fields $chr, $pos1, $pos2 or a gene name / gene list
    snapshot = NULL, # file path to save image to
    track.view = NULL,
    new = FALSE,
    reset = FALSE,
    wkspace = 'PanLungWGS',
    host = Sys.getenv('IGV_HOST'),
    mac = !grepl('(^cga)|(node\\d+)', host), ##
    rawpaths = FALSE,
    sort.locus = NULL,
    gsub.paths = list(
        ## c('~/', '~/home/'),
        ## c('/gpfs/internal/', '/internal/'),
        ## c('/data/analysis/', '/analysis/'),
        ## c('/data/research/mski_lab/data', '/data'),
        ## c('/data/analysis/', '/analysis/'),
        ## c('/seq/picard_aggregation/', '/Volumes/seq_picard_aggregation/'),
        ## c('/xchip/singtex/', '/Volumes/xchip_singtex/'),
        ## c('/cga/meyerson/', '/Volumes/cga_meyerson/'),
        ## c('/seq/tier3b/', '/Volumes/seq_tier3b/'),
        ## c('/xchip/cle/', '/Volumes/xchip_cle/'),
        ## c('/xchip/cga/', '/Volumes/xchip_cga/'),
        ## c('/xchip/beroukhimlab/', '/Volumes/xchip_beroukhimlab/'),
        ## c('/cgaext/tcga/', '/Volumes/cgaext_tcga/'),
        ## c(Sys.getenv('HOME'), 'HOME')
        ),
    port = Sys.getenv('IGV_PORT')
)
{
    if (reset)
        {
            cat('Resetting - closing all connections\n')
            closeAllConnections()
        }

    if (nchar(port)==0)
        port = NULL

    if (is.null(port))
        port = 60151

    if (is.null(host))
        host = ''

    if (nchar(host)==0)
        {
            warning('IGV_HOST field is empty and host is not provided, defaulting to mskilab as host')
            host = 'mskilab'
        }

    con = tryCatch(suppressWarnings(socketConnection(host = host, port = port, open="r+", blocking = TRUE)), error =
           function(x) (stop(sprintf('IGV does not appear to be running on host %s, port %s.  Please start IGV (v2 or later) in external shell for host %s, then retry igv command.', host, port, Sys.getenv('HOST')))))

    WEB.DIR = '~/public_html'
    IGV.TMP = sprintf('.igvtmp.%s.%s', host, port)
    if (new)
        {
            if (file.exists(paste(WEB.DIR, IGV.TMP, sep = '/')))
                system(sprintf('rm -rf %s', paste(WEB.DIR, IGV.TMP, sep = '/')))
            igv.cmd('new', con)
        }

    if (!is.null(paths))  ## if annotation is specified will treat id as a sample id and try to load the bam pointed to by the corresponding annotation, using wkspac
        {
            if (is.list(paths))
                paths = unlist(paths)

            if (is.vector(paths))
                if (any(ix <- !is.na(paths) & !nchar(paths)==0))
                    {
                        paths = paths[ix]

                        if (rawpaths)
                            new.paths = paths
                        else
                            new.paths = normalizePath(paths);

                        if (mac & !rawpaths)
                            {
                                cat('converting server paths for MacOS IGV..\n')
                                for (g in gsub.paths)
                                    {
                                        new.paths = gsub(gsub('([^\\w])', '\\\\\\1', g[1], perl = TRUE), g[2], new.paths)
                                    }
                            }

                        for (p in new.paths)
                            igv.cmd(paste(c('load', p), collapse = ' '), con)
                    }
                else
                    warning('paths input must be character vector or list')
        }

    if (!is.null(gr))
        {
            system(paste('mkdir -p', IGV.TMP))
            grname = gsub('[^\\w]+', '_',deparse(substitute(gr)), perl = TRUE)
            if (is(gr, 'GRanges'))
                if (!is.null(gr$score))
                    {
                        gr = gr[!is.na(gr$score), 'score']
                        gr = gr[!is.infinite(gr$score), ]
                        gr = coverage(gr, weight = gr$score)
                        format = 'bw'
                    }
                else
                    format = 'gff3'
            else if (is(gr, 'GRangesList'))
                format = 'gff3'
            else if (is(gr, 'GAlignments'))
                format = 'bam'
            else
                {
                    close(con)
                    stop('Unrecognized type given for gr input')
                }

            dirname = paste(IGV.TMP, paste(gsub('\\.', '', runif(1)), sep = '.'), sep = '/')
            localdir = paste(normalizePath(WEB.DIR), dirname, sep = '/')
            webdir = paste(paste(BASE.URL, '/~', Sys.getenv('USER'), sep = ''), dirname, sep = '/')
            system(paste('mkdir -p', localdir))
            localfn = paste(localdir, paste(grname, format, sep = '.'), sep = '/')
            url = paste(webdir, paste(grname, format, sep = '.'), sep = '/')

            if (format != 'bw')
                if (!is.null(names(gr)))
                    if (is.null(gr$ID))
                        gr$ID = names(gr)

            if (length(gr) > 0) {
                mc <- mcols(gr)
                mcols(gr) <- mc[,!sapply(seq(ncol(mc)), function(x) class(mc[,x])) %in% c('DNAStringSet')]
            }

            if (length(gr)>0)
                export(gr, con = localfn, format = format)

            Sys.sleep(2)
            igv.cmd(paste('load', url), con)
        }

    if (!is.null(loci))
        {
            if (!is.character(loci))
                {
                    names(loci) = NULL
                    if (is(loci, 'GRangesList'))
                        loci = unlist(loci)

                    names(loci) = NULL

                    loci = as.data.frame(loci[, c()])
                    sort.locus = loci[1,]
                    sort.locus$start = sort.locus$end = (sort.locus$start + sort.locus$end)/2
                }
        }

    if (!is.null(loci))
        if (is.list(loci))
            {
                loci = standardize_segs(loci)
                igv.cmd(paste('goto', paste(paste(loci$chr, ':', loci$pos1, '-', loci$pos2, sep = ""), collapse = ", ")), con)
            }
        else
            igv.cmd(paste('goto', paste(loci, collapse = " | ")), con);

    ## trackview commands
    cat(sprintf("%s: %s\n",track.view,class(track.view)))
    if (!is.null(track.view)) {
        track.view <- tolower(track.view)
        if (track.view %in% c('expand','squish','collapse')) {
            igv.cmd(track.view, con)
        }
        else {
            warning(sprintf("unrecognixed track.view '$s'\n"))
        }
    }

    if (!is.null(sort.locus))
        {
            sort.locus = standardize_segs(sort.locus)
            if (is.null(sort.locus$pos2))
                sort.locus$pos2 = sort.locus$pos1

            if (is.null(sort.locus$option))
                sort.locus$option = 'base'

            cmd = paste('sort', sort.locus$option, paste(sort.locus$chr, ":", sort.locus$pos1, sep = ""))
                                        #     print(cmd)
            igv.cmd(cmd, con)
        }

    if (!is.null(snapshot))
        {
            igv.cmd(paste('snapshotDirectory', gsub('^\\~', '$HOME', dirname(snapshot))), con)
            igv.cmd(paste('snapshot', basename(snapshot)), con)
        }

    close(con)
}

#####################
# igv.cmd
#
# Low level utility function to send and echo commands to igv connection "con" via writeLines
#
#####################
igv.cmd = function(cmd, con, verbose = TRUE)
{
    if (verbose)
        cat(sprintf('Sending to IGV via connection %s: %s\n', summary(con)$description, cmd))
    writeLines(cmd, con)
    response <- readLines(con,n=1)
    if (verbose)
        cat(sprintf("IGV replies: %s\n",response))
}



#' @name igv.loci
#' @title igv.loci
#'
#' Wrapper for igv function to dump table + screenshots for individual GRanges loci
#' that have sample column (default Tumor_Sample_Barcode) that is a key into
#' ind data.table where columns matching col.string are fetched and plotted
#'
#' IGV host and port are taken from environment variables
#' @description
#' Wrapper for igv function to visualize a given (single) mutation in a maf file
#'
#' @export
#' @author Marcin Imielinski
igv.loci = function(mut, ## GRanges of loci
    ind, ## keyed individual table containing all the paths to send to igv
    out.path, ## where to dump files and mutations
    sample.key = 'Tumor_Sample_Barcode', ## column in mut which to use to key ind
    sleep = 30, ## how much to sleep per screenshot
    window = 400, # bp window around mutation
    host = Sys.getenv('IGV_HOST'),
    port = Sys.getenv('IGV_PORT'),
    overwrite = FALSE,
    snapshots = TRUE, ## if false then won't dump snapshots (ie if only want to dump table)
    verbose = TRUE
)
  {
      if (is.null(key(ind)))
          stop('ind must be keyed')

      if (!(sample.key %in% names(values(mut))))
          stop(sprintf('Sample key %s does not exist in mutations'))

      if (!file.exists(out.path))
          system(paste('mkdir -p', paste(out.path, 'imgs', sep = '/')))

      nms = names(mut)
      if (is.null(nms))
          nms = 1:length(mut)

      out.path = normalizePath(out.path)
      img.path = paste(out.path, '/imgs/', nms, '.png', sep = '')
      img.path.web = paste('imgs/', nms, '.png', sep = '')
      mut$img.link = paste("<a href = \"", img.path.web,  "\"> img </a>", sep = '')
      out.path.htab = paste(out.path, 'mutations.html', sep = '/')
      out.path.tab = paste(out.path, 'mutations.txt', sep = '/')

      if (verbose)
          {
              cat('Taking', length(mut), 'snapshots across', ncol(ind),
                  'annotations (including', paste(setdiff(names(ind), key(ind))[1:(min(ncol(ind)-1))], collapse = ', '))
              cat(' using IGV instance running on host', host, ', port(s)', paste(port, collapse = ' and '), '\nDumping tables to', out.path.htab, '\nand', out.path.tab, '\n')
              cat('\tGiving you a chance to think things over...\n')
          }

      write.htab(mut, out.path.htab)
      write.tab(mut, out.path.tab)
      mut$img.path = img.path
      mut$img.path.web = img.path.web

      if (snapshots)
          {
              port.j = 1 ## if several ports

              if (overwrite)
                  ix = 1:length(mut)
              else
                  {
                      ix = which(!file.exists(mut$img.path))
                      cat(length(mut)-length(ix), 'of', length(mut), 'image files already exist\n')
                  }

              if (length(ix)>0)
                  for (i in 1:length(ix))
                      {
                          if (!is.null(sleep))
                              if (!is.na(sleep))
                                  if (sleep>0)
                                      Sys.sleep(sleep)
                      port.j = ifelse(port.j==length(port), 1, port.j+1)
                      this.port = port[port.j]
                      cat('Plotting igv for mutation', ix[i], 'which is', i, 'of', length(ix), '\n')
                      this.mut = mut[ix[i]]
                      igv(reset = TRUE, host = host, port = this.port)
                      igv(new = TRUE, host = host, port = this.port)
                      igv(ind[values(this.mut)[, sample.key]], host = host, port = this.port)
                      igv(loci = this.mut + window, host = host, port = this.port)
                      cat('\ttaking snapshot and dumping to', this.mut$img.path, '\n')
                      igv(snapshot = this.mut$img.path, host = host, port = this.port)
                  }
          }
  }


#' @name dcast2
#' @title dcast.data.table but allows vector arguments for value.var,
#' @description
#' if value.var is a vector then will combine the right hand side column names with each element of value.var
#' in a merged cast table
#'
#' @export
#' @author Marcin Imielinski
dcast2 = function(data, formula, ..., value.var = NULL,
    fun.aggregate = function(x) if (length(x)<=1) x[1] else paste(x, collapse = ','), sep = '_')
{
    tmp = strsplit(as.character(formula), ' \\+ ')
    left.side = tmp[[2]]
    right.side = tmp[[3]]

        terms = sapply(unlist(as.list(attr(terms(formula), "variables"))[-1]), as.character)
        if (is.null(value.var))
            value.var = setdiff(colnames(data), c(left.side, right.side))
        dt = lapply(value.var, function(x)
            {
                if (is.data.table(data))
                    d = dcast.data.table(data, formula,..., fun.aggregate = fun.aggregate, value.var = x)
                else
                    {
                        d = as.data.table(dcast(data, formula,  ..., fun.aggregate = fun.aggregate, value.var = x))
                    }
                new.cols = setdiff(colnames(d), c(key(d), left.side, right.side))
                setnames(d, new.cols, paste(new.cols, x, sep = sep))
                return(d)
            })
        out = dt[[1]]
        if (length(dt)>1)
            for (d in dt[-1])
                out = merge(out, d)
        return(out)
    }

#' @name .igv_host
#' @title  .igv_host
#' @rdname igv_host
#' @description
#'
#' several settings for igv host
#'
#' @export
#' @author Marcin Imielinski
.igv_host = function(h)
    {
        if (grepl('^cga', h))
            Sys.setenv(IGV_HOST = h, IGV_PORT = "60199")
        else if (grepl('(retina)|(laptop)', h))
            Sys.setenv(IGV_HOST = 'wm19d-4a%!%1', IGV_PORT = "60151")
        else if (grepl('air', h))
            Sys.setenv(IGV_HOST = 'vpn5-124', IGV_PORT = "60151")
    }

#' Round a set of GRanges to another set
                                        #
#' "rounds" a set of query ranges Q to a given interval set S using the following rule:
#' 1) If q in Q is partially / fully within S then return intersection of q with S.
#' 2) If q intersects multiple ranges in S and \code{up = F} then return the "first" range, otherwise the last range
#' 3) If q in Q is fully outside of S (ie fully inside not S) then return the \code{start-1} (if \code{up = T}) or \code{end+1} (if \code{up = F})
#' of the matching range in not S
#'
#' @param Q Query \code{GRanges} (strand is ignored)
#' @param S Subject \code{GRanges} (strand is ignored)
#' @param up [default TRUE] See description.
#' @param parallel [default FALSE] If \code{TRUE}, assumes Q and S are same length and this analysis is only performed between the corresponding Q and S pairs.
#' @return Rounded \code{GRanges}
#' @examples
#' \dontrun{query   <- GRanges(1, IRanges(c(100,110),width=201), seqinfo=Seqinfo("1", 500))
#' subject <- GRanges(1, IRanges(c(160,170),width=201), seqinfo=Seqinfo("1", 500))
#' gr.round(query, subject)}
#' @export
gr.round = function(Q, S, up = TRUE, parallel = FALSE)
{
    str = strand(Q)
    Q = gr.stripstrand(Q)
    S = gr.stripstrand(S)
    nS = gaps(S);
    QS = gr.findoverlaps(Q, S)
    tmp = gr.findoverlaps(Q, nS)
    QnotS = nS[tmp$subject.id]
    QnotS$query.id = tmp$query.id

    if (parallel)
    {
        QS = QS[QS$query.id==QS$subject.id]
        QnotS = QnotS[QnotS$query.id==QnotS$subject.id]
    }

    if (up)
        suppressWarnings(end(QnotS) <- start(QnotS) <- end(QnotS)+1)
    else
        suppressWarnings(start(QnotS) <- end(QnotS) <- start(QnotS)-1)

    suppressWarnings(out <- sort(grbind(QS, QnotS)))

    if (up)
    {
        out = rev(out)
        out = out[!duplicated(out$query.id)]
    }
    else
        out = out[!duplicated(out$query.id)]

    out = out[order(out$query.id)]
    values(out) = values(Q)
    names(out) = names(Q)
    strand(out) = str;
    return(out)
}


#' Convert from chrXX to numeric format
#'
#' Convert from chrXX to numeric format
#' @param x factor, Rle or character vector with chromosome names
#' @param xy Flag to convert M to 25, Y to 24 and X to 23. Default FALSE
#' @return character vector with xy=FALSE, or numeric vector with xy=TRUE
#' @export
##########################
chr2num = function(x, xy = FALSE)
  {
      if (inherits(x, 'factor') | inherits(x, 'Rle'))
            x = as.character(x)

     out = gsub('chr', '', x);

     if (!xy)
            out = as.numeric(gsub('M', '25', gsub('Y', '24', gsub('X', '23', out))))

     return(out)
       }



#' gr.refactor
#'
#' Takes a pile of ranges gr and new seqnames "sn" (either of length 1 or
#' of length(gr)) and returns a gr object with the new seqnames and same
#' widths and new start coordinates.  These coordinates are determined by placing
#' each gr on the corresponding chromosome at 1 + gap after previous gr (or at 1)
#' @param gr \code{GRanges} to refactor
#' @param sn character vector of new seqnames
#' @param gap Default 0
#' @param rev Default FALSE
#' @name gr.refactor
#' @export
gr.refactor = function(gr, sn, gap = 0, rev = FALSE)
{
    if (is.factor(sn))
        slev = levels(sn)
    else
        slev = unique(sn);

    sn = cbind(as.character(start(gr)), as.character(sn))[,2]
    w = width(gr)
    gap = pmax(cbind(gap, w)[,1], 0);
                                        #    gap = pmax(gap, 0);
                                        #    starts = levapply(width(gr), sn, function(x) cumsum(gap+c(1, x[1:(length(x)-1)]))[1:length(x)]-gap)

    starts = levapply(1:length(w), sn, function(x) cumsum(gap[x] + c(1, w[x[1:(length(x)-1)]])[1:length(x)])-gap[x])
    ir = IRanges(starts, width = width(gr))

                                        # figure out seqlevels so that the order matches seqlevels of gr with
                                        # derivative chromosomes "next" to their original
    sl = aggregate(end(ir)+gap, by = list(sn), FUN = max); sl = structure(sl[,2], names = sl[,1])

                                        # reorder and add any missing levels
    oth.names = setdiff(slev, names(sl))
    if (length(oth.names)>0)
        sl[oth.names] = NA
    sl = sl[slev]

    out = GRanges(sn, ir, strand = strand(gr), seqlengths = sl)
    values(out) = values(gr);

    return(out)
}



#' grl.span
#'
#' Returns GRanges object representing the left / right extent of each GRL item.  In case of "chimeric" GRL items (ie that map
#' to two chromosomes) there are two options:
#' (1) specify "chr" chromosome as argument to subset GRL's that are on that chromosome, and compute GRL extents from this, any GRL
#'     full outside of that chromosome will get a 0 width GRL
#' (2) (default) allow chimeric GRL items to get an extent that is with respect to the first chromosome in that GRL
#'
#' If a grl item contains ranges that lie on different chromosomes, then corresponding grange will have chromosome "NA" and IRange(0, 0)
#' @param grl \link{GRangesList} to query
#' @param chr [Default NULL]
#' @param ir [Default FALSE]
#' @param keep.strand [Default TRUE]
#' @name grl.span
#' @export
grl.span = function(grl, chr = NULL, ir = FALSE, keep.strand = TRUE)
{
    if (is.null(names(grl)))
        names(grl) = 1:length(grl);

    tmp = tryCatch(as.data.frame(grl), error = function(e) e)

    if (inherits(tmp, 'error')) ## gr names are screwy so do some gymnastics
    {
        if (is.null(names(grl)))
            names.grl = 1:length(grl)
        else
            names.grl = names(grl);

        element = as.character(Rle(names.grl, sapply(grl, length)))
        tmp.gr = unlist(grl)
        names(tmp.gr) = NULL;
        tmp = as.data.frame(tmp.gr);
        tmp$element = element;
    }

    if (is.null(chr))
    {
        chrmap = aggregate(formula = seqnames ~ element, data = tmp, FUN = function(x) x[1]);
        chrmap = structure(as.character(chrmap[,2]), names = chrmap[,1])

        if (keep.strand)
        {
            strmap = aggregate(formula = as.character(strand) ~ element, data = tmp, FUN =
                                                                                         function(x) {y = unique(x); if (length(y)>1) return('*') else y[1]})
            strmap = structure(as.character(strmap[,2]), names = strmap[,1])
            str = strmap[names(grl)];
        }
        else
            str = '*'

        tmp = tmp[tmp$seqnames == chrmap[tmp$element], ]; ## remove all gr from each GRL item that don't map to the chr of the first gr
        chr = chrmap[names(grl)];
        out.gr = GRanges(chr, IRanges(1,0), seqlengths = seqlengths(grl), strand = str)
    }
    else
    {
        if (length(chr)>1)
            warning('chr has length greater than 1, only the first element will be used')
        tmp = tmp[tmp$seqnames == chr[1], ]
        out.gr = rep(GRanges(chr, IRanges(1, 0)), length(grl)) # missing values
    }

    if (nrow(tmp)>0)
    {
        tmp = split(GRanges(tmp$seqnames, IRanges(tmp$start, tmp$end)), tmp$element)
        out.gr[match(names(tmp), names(grl))] = GRanges(chr[names(tmp)],
                                                        IRanges(sapply(start(tmp), min), sapply(end(tmp), max)), strand = strand(out.gr)[match(names(tmp), names(grl))]);
        names(out.gr) = names(grl)
    }
    return(out.gr)
}

#' Create GRanges from Seqinfo
#'
#' Creates a genomic ranges from seqinfo object
#' ie a pile of ranges spanning the genome
#' @param si Seqinfo object
#' @param strip.empty Don't know. Default FALSE
#' @examples
#' \dontrun{si <- Seqinfo(names(hg_seqlength(), hg_seqlengths()))
#' seqinfo2gr(si)}
#' @export
seqinfo2gr <- function(si, strip.empty = FALSE)
{
    if (is(si, 'vector')) ## treat si as seqlengths if vector
        si = Seqinfo(seqlengths = si, seqnames = names(si))
    else if (!is(si, 'Seqinfo'))
        si = seqinfo(si)

    sl = seqlengths(si)
    sn = seqnames(si);
    sl[is.na(sl)] = 0;

    if (strip.empty)
    {
        sn = sn[sl!=0];
        sl = sl[sl!=0];
    }

    sigr = GRanges(sn, IRanges(rep(1, length(sl)), width = sl), seqlengths = seqlengths(si), strand = rep('+', length(sl)))
    names(sigr) = sn;

    return(sigr)
}

#' gr.tostring
#'
#' dumps out a quick text representation of a gr object (ie a character vector)
#' @param gr \code{GRanges}
#' @param places Number of decimal places. Default 2
#' @param interval Default 1e6
#' @param unit Default "MB"
#' @param prefix Default "chr"
#' @return text representation of input
#' @name gr.tostring
#' @export
gr.tostring = function(gr, places = 2, interval = 1e6, unit = 'MB', prefix = 'chr')
{
    p1 = round(start(gr)/interval, places);
    p2 = round(end(gr)/interval, places);
    return(paste(prefix, as.character(seqnames(gr)), ':', p1, '-', p2, ' ', unit, sep = ''));
}

#' Convert data.table to GRanges
#'
#' Takes as input a data.table which must have the fields: start, end, strand, seqnames.
#' All of the remaining fields are added as meta data to the GRanges
#' @param dt data.table to convert to GRanges
#' @return GRanges object of length = nrow(dt)
#' @examples
#' \dontrun{r <- dtgr(data.table(start=1, seqnames="X", end=2, strand='+'))}
#' @export
dtgr <- function(dt) {

    rr <- IRanges(dt$start, dt$end)
    if (!'strand' %in% colnames(dt))
        dt$strand <- '*'
    sf <- factor(dt$strand, levels=c('+', '-', '*'))
    ff <- factor(dt$seqnames, levels=unique(dt$seqnames))
    out <- GRanges(seqnames=ff, ranges=rr, strand=sf)
    if (inherits(dt, 'data.table'))
        mc <- as.data.frame(dt[, setdiff(colnames(dt), c('start', 'end', 'seqnames', 'strand')), with=FALSE])
    else if (inherits(dt, 'data.frame'))
        mc <- as.data.frame(dt[, setdiff(colnames(dt), c('start', 'end', 'seqnames', 'strand'))])
    else
        warning("Needs to be data.table or data.frame")
    if (nrow(mc))
        mcols(out) <- mc
    return(out)
}

#' Filters GRangesList to only include ranges in the specified window
#'
#' (this is different from %in% which does not remove non matching ranges from the grls)
#'
#' does not return list in necessarily same order
                                        # @param grl \link{GRangesList} to filter
                                        # @param windows \link{GRanges} windows to keep
#' @name grl.filter
#' @export
grl.filter = function(grl, windows)
{
    tmp = as.data.frame(grl);
    tmp = tmp[seg.on.seg(tmp, windows), ]
    FORBIDDEN = c('seqnames', 'start', 'end', 'strand', 'ranges', 'seqlevels', 'seqlengths', 'isCircular', 'genome', 'width', 'element');
    gr.metadata = tmp[, setdiff(colnames(tmp), FORBIDDEN)];

    if (!is.null(dim(gr.metadata)))
        out.grl = split(GRanges(tmp$seqnames, IRanges(tmp$start, tmp$end), seqlengths = seqlengths(grl), gr.metadata,
                                strand = tmp$strand), tmp$element)
    else
        out.grl = split(GRanges(tmp$seqnames, IRanges(tmp$start, tmp$end), seqlengths = seqlengths(grl),
                                strand = tmp$strand), tmp$element);

    values(out.grl) = values(grl)[match(names(out.grl), names(grl)), ]
    return(out.grl);
}
#' grl.split
#'
#' splits GRL's with respect to their seqnames and strand (default), returning
#' new grl whose items only contain ranges with a single seqname / strand
#'
#' can also split by arbitrary (specified) genomic ranges value fields
#' @param grl \code{GRangesList} to split
#' @param seqname Default TRUE
#' @param strand Default TRUE
#' @param values columns of values field in grl
#' @name grl.split
#' @export
grl.split = function(grl, seqname = TRUE, strand = TRUE,
                     values = c() # columns of values field in grl
                     )
{
    ele = tryCatch(as.data.frame(grl)$element, error = function(e) e)
    if (inherits(ele, 'error'))
    {
        if (is.null(names(grl)))
            nm = 1:length(names(grl))
        else
            nm = names(grl)

        ele = unlist(lapply(1:length(grl), function(x) rep(nm[x], length(grl[[x]]))))
    }

    gr = unlist(grl)
    names(gr) = NULL;

    by = ele;
    if (seqname)
        by = paste(by, seqnames(gr))

    if (strand)
        by = paste(by, strand(gr))

    values = intersect(names(values(gr)), values);
    if (length(values)>0)
        for (val in values)
            by = paste(by, values(gr)[, val])

    out = split(gr, by);
    names(out) = ele[!duplicated(by)]

    values(out) = values(grl[ele[!duplicated(by)]])

    return(out)
}

#' Get GC content from reference genome
#'
#' Uses BSgenome package to compute gc content for a collection of segments in seg data frame ($chr, $start, $end or $chr, $pos1, $pos2 or $chr, $begin, $end)
#' Returns vector of gc content of length nrow(segs).
#' @param segs Segment data frame to pull gc from
#' @param bs_genome A \code{\link{BSgenome}} object. Perhaps \code{BSgenome.Hsapiens.UCSC.hg19::Hsapiens}
#' @export
#' @name gc_content
gc_content = function(segs, bs_genome) ##build = 'hg19')
{
    segs = standardize_segs(segs, chr = TRUE);

    ## NEW
    tmp = getSeq(bs_genome, segs$chr, segs$pos1, segs$pos2, as.character = TRUE)
    ##     if (build == 'hg19') {
    ##       if (requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly=TRUE)) {
    ##         tmp = getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, segs$chr, segs$pos1, segs$pos2, as.character = TRUE)
    ##       }
    ##     }
    ##     else if (build == 'hg18') {
    ##       if (requireNamespace("BSgenome.Hsapiens.UCSC.hg18", quietly=TRUE)) {
    ##         tmp = getSeq(BSgenome.Hsapiens.UCSC.hg18::Hsapiens, segs$chr, segs$pos1, segs$pos2, as.character = TRUE)
    ##       }
    ##     }
    ##     else
    ##       stop('gc_content: hg build not recognized');

    ## OLD
    ##    if (build == 'hg19')
    ##      library(BSgenome.Hsapiens.UCSC.hg19)
    ##    else if (build == 'hg18')
    ##     library(BSgenome.Hsapiens.UCSC.hg18)
    ##  else
    ##   stop('gc_content: hg build not recognized');

    ## tmp = getSeq(Hsapiens, segs$chr, segs$pos1, segs$pos2, as.character = T)

    return(as.numeric(sapply(gregexpr('[GC]', tmp), length)/sapply(tmp, nchar)))
}

#' import.ucsc
#'
#' wrapper around rtracklayer import that
#' (1) handles "as" formats
#' (2) has additional flag chrsub to sub in 'chr' in selection, and then sub it out of the output
#' @name import.ucsc
#' @export
import.ucsc = function(con, selection = NULL, text, chrsub = TRUE, verbose = FALSE, as = NULL, ...) {
    si = NULL;

    if (verbose)
        cat('importing', as.character(con), '\n')

    if (grepl('(\\.bw)|(\\.bigwig)', con, ignore.case = TRUE))
    {

        if (is.character(con))
            f = BigWigFile(normalizePath(con))
        else
            f = con

        si = tryCatch(seqinfo(f), error = function(con) NULL)
    }
    else if (grepl('\\.wig', con, ignore.case = TRUE))
    {

        if (is.character(con))
            f = WIGFile(normalizePath(con))
        else
            f = con

        si = tryCatch(seqinfo(f), error = function(con) NULL)
    }
    else if (grepl('\\.bed', con, ignore.case = TRUE))
    {

        if (is.character(con))
            f = BEDFile(normalizePath(con))
        else
            f = con
                                        #                            si = tryCatch(seqinfo(f), error = function(con) NULL)
        bed.style = T
    }
    else if (grepl('\\.gff', con, ignore.case = TRUE))
    {

        if (is.character(con))
            f = GFFFile(normalizePath(con))
        else
            f = con

        si = tryCatch(seqinfo(f), error = function(con) NULL)
    }
    else if (grepl('\\.2bit', con, ignore.case = T))
    {
        if (is.character(con))
            f = TwoBitFile(normalizePath(con))
        else
            f = con

        si = tryCatch(seqinfo(f), error = function(con) NULL)
    }
    else if (grepl('\\.bedgraph', con, ignore.case = T))
    {

        if (is.character(con))
            f = BEDGraphFile(normalizePath(con))
        else
            f = con
                                        #                           si = tryCatch(seqinfo(f), error = function(con) NULL)
        bed.style = T
    }
    else
        f = con

    if (chrsub & !is.null(si) & !is.null(selection))
        selection = gr.fix(gr.chr(selection), si, drop = T)

    if (class(f) %in% c('BEDFile'))
    {
        if (!is.null(selection))
            out = import(f, selection = selection, asRangedData = FALSE, ... )
        else
            out = import(f, asRangedData = FALSE, ...)
    }
    else
    {
        if (!is.null(selection))
            out = import(f, selection = selection )
        else
            out = import(f)
    }

    if (!is(out, 'GRanges'))
        out = as(out, 'GRanges')

                                        #    if (chrsub & !is.null(si))
    if (chrsub)
        out = gr.sub(out, 'chr', '')

    if (verbose)
        cat('Finished importing', as.character(con), '\n')

    return(out)
}


#####################################################################
#
#
# $$\      $$\ $$\                                      $$\     $$\ $$\
# $$$\    $$$ |\__|                                     $$ |    \__|$$ |
# $$$$\  $$$$ |$$\  $$$$$$$\  $$$$$$$\       $$\   $$\$$$$$$\   $$\ $$ |
# $$\$$\$$ $$ |$$ |$$  _____|$$  _____|      $$ |  $$ \_$$  _|  $$ |$$ |
# $$ \$$$  $$ |$$ |\$$$$$$\  $$ /            $$ |  $$ | $$ |    $$ |$$ |
# $$ |\$  /$$ |$$ | \____$$\ $$ |            $$ |  $$ | $$ |$$\ $$ |$$ |
# $$ | \_/ $$ |$$ |$$$$$$$  |\$$$$$$$\       \$$$$$$  | \$$$$  |$$ |$$ |
# \__|     \__|\__|\_______/  \_______|       \______/   \____/ \__|\__|
#
#
# Misc util
################

#' Improved rbidn for intersecting columns of data.frames or data.tables
#'
#' like rbind, but takes the intersecting columns of the dfs
#' rrbind = function(df1, df2, [df3 ... etc], )
#' @param ... list of data frames to concatenate
#' @param union if union flag is used then will take union of columns (and put NA's for columns of df1 not in df2 and vice versa). Default TRUE
#' @param as.data.table [Default FALSE] return as a \link{data.table}
#' @export
rrbind2 = function(..., union = T, as.data.table = FALSE)
{
    dfs = list(...);  # gets list of data frames
    dfs = dfs[!sapply(dfs, is.null)]
    dfs = dfs[sapply(dfs, ncol)>0]
    names.list = lapply(dfs, names);
    cols = unique(unlist(names.list));
    unshared = lapply(names.list, function(x) setdiff(cols, x));
    ix = which(sapply(dfs, nrow)>0)
    ## only call expanded dfs if needed
    if (any(sapply(unshared, length) != 0))
        expanded.dts <- lapply(ix, function(x) {
            tmp = dfs[[x]]
            if (is.data.table(dfs[[x]]))
                tmp = as.data.frame(tmp)
            tmp[, unshared[[x]]] = NA;
            return(as.data.table(as.data.frame(tmp[, cols])))
        })
    else
        expanded.dts <- lapply(dfs, function(x) as.data.table(as.data.frame(x)))

    ## convert data frames (or DataFrame) to data table.
    ## need to convert DataFrame to data.frmae for data.table(...) call.
    ## structure call is way faster than data.table(as.data.frame(...))
    ## and works on data.frame and DataFrame
                                        #    dts <- lapply(expanded.dfs, function(x) structure(as.list(x), class='data.table'))
                                        #   rout <- data.frame(rbindlist(dts))

    rout <- rbindlist(expanded.dts)
    if (!as.data.table)
        rout = as.data.frame(rout)

    if (!union)
    {
        shared = setdiff(cols, unique(unlist(unshared)))
        rout = rout[, shared];
    }

    return(rout)
}

#' Identify matches between query and dictionary
#'
#' Wrapper around matchPdict to identify matches between a query
#' string query and dictionary dict (both BString objects or subclasses)
#'
#' @param query Query
#' @param dict Dictionary
#' @param midpoint Flag for output the coordinates of the match as the location,
#'   where the midpoint of the dict string matches the given query. Default FALSE
#' @return a vector of indices of length width(query) that contains
#' indices of the (starting) dictionary in the query string
#' @export
match.bs = function(query, dict, midpoint = FALSE)
{
    names(dict) = as.character(1:length(dict))

    tmp = sort(unlist(matchPDict(dict, query)))
    out = rep(NA, length(query))

    if (!midpoint)
        out[start(tmp)] = as.numeric(names(tmp))
    else
        out[floor((start(tmp)+end(tmp))/2)] = as.numeric(names(tmp))

    return(out)
}

#' Wrapper around BSgenome call
#'
#' Retreives either the BSgenome hg18 or hg19 genome by default.  Requires packages
#' BSgenome.Hsapiens.UCSC.hg19 for hg19 and BSgenome.Hsapiens.UCSC.hg19 for hg18.
#'
#' If fft = TRUE, can also also return the hg19 ffTrack (requires that the file exists)
#' Requires the existence of environment variable HG.FFT pointing to ffTrack .rds file..
#'
#' @param hg19 Logical whether to return hg18 or hg19 BSgenome. Default TRUE
#' @param fft Logical whether to return an ffTrack. Default FALSE
#' @return BSgenome or ffTrack of the genome
#' @export
read_hg = function(hg19 = T, fft = F)
{
    if (fft)
    {
        if (file.exists(Sys.getenv('HG.FFT')))
            REFGENE.FILE.HG19.FFT = Sys.getenv('HG.FFT')
        else if (file.exists('~/DB/ffTracks/hg19.rds'))
            REFGENE.FILE.HG19.FFT = '~/DB/ffTracks/hg19.rds'
        else
            stop("Need to supply environment variable to FFtracked genome or load BSGenome. Env Var: HG.FFT")

        return(readRDS(REFGENE.FILE.HG19.FFT))
    }
    else
    {
        require(BSgenome)
        if (hg19)
            library(BSgenome.Hsapiens.UCSC.hg19)
        else
            library(BSgenome.Hsapiens.UCSC.hg18)
        return(Hsapiens)
    }
}

#' Filter reads by average PHRED score
#' Defines a cutoff score for the mean PHRED quality of a read
#' in a GRanges.
#' @param gr GRanges or data.table of reads that has a \code{qname} and \code{qual} field
#' @param cutoff cutoff score for mean PHRED quality. Default "+"
#' @return GRanges or data.table where reads have mean quality score >= cutoff
#' @export
gr.readfilter <- function(gr, cutoff = '+') {

    cutoff <- as.numeric(charToRaw(cutoff))
    qual   <- as.character(gr$qual)
    logvec <- sapply(qual, function(x) mean(as.numeric(charToRaw(x))) < cutoff)

    qn <- unique(gr$qname[logvec])

    gr <- gr[!(gr$qname %in% qn)]

    return(gr)
}

#' Check if reads are clipped
#'
#' Returns a logical vector of length of the input GRanges that
#' that classifies a read as clipped or not. The user can specify
#' a cutoff value for how many bases need to be clipped.
#' @param gr Granges OR data.table that has \code{cigar} field and \code{qname} field
#' @param clip.cutoff Minimum number of bases that are clipped to call the reads as clipped
#' @return logical of length of input, denoting whether that read is part of a clipped read pair.
#' @export
gr.isclip <- function(gr, clip.cutoff=10) {
    if (inherits(gr, 'GRanges') && length(gr)==0)
        return(logical(0))
    if (inherits(gr, 'data.table') && nrow(gr) == 0)
        return(logical(0))

    if (inherits(gr, 'GRanges'))
        nm <- names(mcols(gr))
    else
        nm <- colnames(gr)
    if (any(!('cigar' %in% nm)))
        stop('gr.isclip: reads need flag and cigar')
    cig <- countCigar(gr$cigar)
    return(cig[,"S"] >= clip.cutoff)
    ##logvec <- grepl('[0-9][0-9]S', gr$cigar)
    ##logvec[is.na(logvec)] <- FALSE
    ##return(logvec)
}


#' Checks if reads are discordant
#'
#' Returns a logical vector denoting if a read is discordant.
#' There is only a minimum absolute isize, and any read below this is
#' not considered discordant. This will return logicals based on read pairs
#' @param gr Granges OR data.table that has \code{isize} field and \code{qname} field
#' @param isize Minimum insert size required to call dis<cordant. Default 1000
#' @param unmap.only Find only pairs with an unmapped read
#' @return logical vector of length of input, denoting each read as discordant or not
#' @export
gr.isdisc <- function(gr, isize=1000, unmap.only=FALSE) {

    if (inherits(gr, 'GRanges') && length(gr)==0)
        return(logical(0))
    if (inherits(gr, 'data.table') && nrow(gr)==0)
        return(logical(0))

    if (inherits(gr, 'GRanges'))
        nm <- names(mcols(gr))
    else
        nm <- colnames(gr)

    if (any(!(c('isize') %in% nm )))
        stop('gr.isdigsc: reads need flag and cigar')

    if (inherits(gr, 'GRanges'))
        st <- start(gr) == 1
    else
        st <- gr$start == 1
    if (unmap.only)
        logvec <- bitAnd(gr$flag, 8) != 0 | st  # last two get reads with unmapped mates and reads that are unmapped, respectively
    else
        logvec <- abs(gr$isize) >= isize | gr$isize==0 | bitAnd(gr$flag, 8) != 0 | st  # last two get reads with unmapped mates and reads that are unmapped, respectively
    logvec[is.na(logvec)] <- FALSE
    qn <- gr$qname[logvec]
    isdisc <- gr$qname %in% qn
    return(isdisc)
}



#' Return windows with minimal coverage
#'
#' Takes a set of GRanges and removes any ranges that
#' don't have a minimal coverage value. If you give it
#' a GRangesList, you will get back an unlisted GRanges.
#'
#' @param gr \code{GRanges} to filter
#' @param min.cov Minimum number of overlaps to keep. Default 2
#' @param buffer Add a buffer to the ranges when computing overlaps. Default 0
#' @param ignore.strand Ignore the strand when comparing overlaps. Default TRUE
#' @param pintersect Force the pintersect option for \link{gr.findoverlaps}
#' @return GRanges
#' @export
gr.mincov <- function(gr, min.cov=2, buffer=0, ignore.strand=TRUE, pintersect=FALSE) {

    if (inherits(gr, 'GRangesList'))
        gr <- unlist(gr)
    if (!inherits(gr, 'GRanges'))
        stop('gr.mincov: Requires a GRanges input')

    gr2 <- gr.pad(gr, buffer)

    tab <- table(gr.findoverlaps(gr2, gr2, ignore.strand=ignore.strand, pintersect=pintersect)$subject.id)
    winkeep <- as.numeric(names(tab[tab >= min.cov]))

    return(gr[winkeep])

}


#' Nice padding
#'
#' @return GRanges
#' @keywords internal
gr.pad = function(gr, pad)
{
    start(gr) = pmax(1, start(gr)-pad)
    en = pmin(seqlengths(gr)[as.character(seqnames(gr))], end(gr)+pad)
    end(gr) = ifelse(is.na(en), end(gr)+pad, en)
    return(gr)
}

#' gr2grl
#' Quick way to make grl from list of indices into a GRanges gr
#'
#' @param gr \code{GRanges} to split
#' @param ix vector to split on
#' @export
gr2grl = function(gr, ix)
{
    out = split(gr[unlist(ix)], sapply(1:length(ix), function(x) rep(x, length(ix[[x]]))))
    if (!is.null(names(ix)))
        names(out) = names(ix)
    return(out)
}


#' Wrapper to base \code{system} function to call system (e.g. perl) from R.
#' The only benefit to this wrapper is a more controlled verbose argument.
#'
#' @author Jeremiah Wala \email{jwala@@broadinstitute.org}
#' @param syscall string containing the system call
#' @param verbose print the syscall to screen, and it's stdout
#' @export
#' @examples
#' # system.call('perl s/[0-9]+//g file1 > file2')
system.call <- function(syscall, verbose=T) {
    if (verbose)
        print(syscall)
    if (verbose)
        system(syscall)
    else {
        system(syscall, intern=TRUE) #, stderr=FALSE, stdin=FALSE)
    }
                                        #system(syscall, intern=T, ignore.stdout=TRUE, ignore.stderr=TRUE)
}

#'
#' identifies events that are in ra1 that do not exist in ra2.
#' Aside from ra1 and ra2, all arguments are sent to ra.overlaps
#'
#' @name ra.overlaps
#' @export
ra.setdiff <- function(ra1, ra2, ...) {

    ro <- ra.overlaps(ra1, ra2, ...)

    in.ra1.only <- setdiff(seq_along(ra1), ro[, 'ra1.ix'])
    return(ra1[unique(in.ra1.only)])

}

#' ra.union
#'
#' returns events in ra1 that are in ra2 also
#'
#' @name ra.union
#' @export
ra.union <- function(ra1, ra2, ...)
    rebturn(ra1[unique(ra.overlaps(ra1, ra2, ...)[, 'ra1.ix'])])

#' gr.all
#'
#' Return a GRanges that holds interavals for all of HG19
#'
#' @param unmap [default F] Optinally add a "unmapped" chr
#' @param M [default F] Include mitochondrial chr
#' @param Y [default T] Include Y chr
#' @return \code{GRanges} object with one element per chromosome
gr.all <- function(unmap=FALSE, M=FALSE, Y=TRUE) {
    gr <- si2gr(gr.tfix(GRanges(1, IRanges(1,1))))

    if (!M)
        gr <- gr[seqnames(gr) != 'M']
    if (!Y)
        gr <- gr[seqnames(gr) != 'Y']


    if (!unmap)
        return(gr[!seqnames(gr) %in% "Unmapped"])
    else
        return(gr)
}


setGeneric('%|%', function(gr, ...) standardGeneric('%|%'))
setMethod("%|%", signature(gr = "GRanges"), function(gr, df) {
    if (is.data.table(df))
        df = as.data.frame(df)
    else if (inherits(df, 'GRanges'))
        df = values(df)
    values(gr) = cbind(values(gr), df)
    return(gr)
})

#' @name toggle_grfo
#' @title toggle data.table vs IRanges find overlaps
#' @description
#'
#' toggles global setting of whether to use data.table vs IRanges find overlaps machinery
#'
#' @export
#' @author Marcin Imielinski
.toggle_grfo = function()
{
    old.val = as.logical(Sys.getenv('GRFO_FOVERLAPS'))
    if (is.na(old.val))
        old.val = FALSE
    Sys.setenv(GRFO_FOVERLAPS = !old.val)
    cat('GRFO_FOVERLAPS is', Sys.getenv('GRFO_FOVERLAPS'), '\n\t...Default gr.findoverlaps behavior will use',
        ifelse(Sys.getenv('GRFO_FOVERLAPS'), 'data.table foverlaps', 'IRanges findOverlaps'), '\n')
}


setGeneric('%|%', function(gr, ...) standardGeneric('%|%'))
setMethod("%|%", signature(gr = "GRanges"), function(gr, df) {
    if (is.data.table(df))
        df = as.data.frame(df)
    else if (inherits(df, 'GRanges'))
        df = values(df)
    values(gr) = cbind(values(gr), df)
    return(gr)
})


subset2 <- function(x, condition) {
    condition_call <- substitute(condition)
    r <- eval(condition_call, x)
    browser()
    x[r, ]
}


##################################
#' @name vaggregate
#' @title vaggregate
#'
#' @description
#' same as aggregate except returns named vector
#' with names as first column of output and values as second
#'
#' Note: there is no need to ever use aggregate or vaggregate, just switch to data.table
#'
#' @param ... arguments to aggregate
#' @return named vector indexed by levels of "by"
#' @author Marcin Imielinski
#' @export
##################################
vaggregate = function(...)
  {
    out = aggregate(...);
    return(structure(out[,ncol(out)], names = do.call(paste, lapply(names(out)[1:(ncol(out)-1)], function(x) out[,x]))))
  }


####################################
#' @name modix
#' @title modix
#'
#' @description
#' Takes integer input ix and projects on to 1-based modulus over base l
#'
#' ie modix(1, 5) -> 1, modix(5, 5) -> 5, modix(6, 5) -> 1
#'
#' @param ix input indices to apply module
#' @param l base of ix
#' @return ((ix-1) mod l) - 1
#' @author Marcin Imielinski
#' @export
####################################
modix = function(ix, l)
  {
    return(((ix-1) %% l)+1)
  }




#' @name affine.map
#' @title affine.map
#' @description
#'
#'
#' affinely maps 1D points in vector x from interval xlim to interval ylim,
#' ie takes points that lie in
#' interval xlim and mapping onto interval ylim using linear / affine map defined by:
#' (x0,y0) = c(xlim(1), ylim(1)),
#' (x1,y1) = c(xlim(2), ylim(2))
#' (using two point formula for line)
#' useful for plotting.
#'
#' if cap.max or cap.min == T then values outside of the range will be capped at min or max
#' @rdname affine-map-methods
#' @author Marcin Imielinski
#' @export
#' @keywords internal
affine.map = function(x, ylim = c(0,1), xlim = c(min(x), max(x)), cap = F, cap.min = cap, cap.max = cap, clip = T, clip.min = clip, clip.max = clip)
{
  #  xlim[2] = max(xlim);
  #  ylim[2] = max(ylim);

  if (xlim[2]==xlim[1])
    y = rep(mean(ylim), length(x))
  else
    y = (ylim[2]-ylim[1]) / (xlim[2]-xlim[1])*(x-xlim[1]) + ylim[1]

  if (cap.min)
    y[x<min(xlim)] = ylim[which.min(xlim)]
  else if (clip.min)
    y[x<min(xlim)] = NA;

  if (cap.max)
    y[x>max(xlim)] = ylim[which.max(xlim)]
  else if (clip.max)
    y[x>max(xlim)] = NA;

  return(y)
}


#' @name alpha
#' @title alpha
#' @description
#' Give transparency value to colors
#'
#' Takes provided colors and gives them the specified alpha (ie transparency) value
#'
#' @author Marcin Imielinski
#' @param col RGB color
#' @keywords internal
#' @export
alpha = function(col, alpha)
{
  col.rgb = col2rgb(col)
  out = rgb(red = col.rgb['red', ]/255, green = col.rgb['green', ]/255, blue = col.rgb['blue', ]/255, alpha = alpha)
  names(out) = names(col)
  return(out)
}

#' Blends colors
#'
#'
#' @param cols colors to blend
#' @keywords internal
#' @name blend
#' @export
blend = function(cols)
{
  col.rgb = rowMeans(col2rgb(cols))
  out = rgb(red = col.rgb['red']/255, green = col.rgb['green']/255, blue = col.rgb['blue']/255)
  return(out)
}


#' @name col.scale
#' @title col.scale
#' @description
#'
#' Assign rgb colors to numeric data
#'
#' Assigns rgb colors to numeric data values in vector "x".. maps scalar values
#' in val.range (default c(0,1)) to a linear color scale of between col.min (default white)
#' and col.max (default black), each which are length 3 vectors or characters.  RGB values are scaled between 0 and 1.
#' Values below and above val.min and val.max are mapped to col.max and col.max respectively
#'
#' @author Marcin Imielinski
#' @export
#' @keywords internal
col.scale = function(x, val.range = c(0, 1), col.min = 'white', col.max = 'black', na.col = 'white',
                     invert = F # if T flips rgb.min and rgb.max
)
{

  ## NOTE fix
  error = NULL

  if (!is.numeric(col.min))
    if (is.character(col.min))
      col.min = col2rgb(col.min)/255
    else
      error('Color should be either length 3 vector or character')

    if (!is.numeric(col.max))
      if (is.character(col.max))
        col.max = col2rgb(col.max)/255
      else
        error('Color should be either length 3 vector or character')

      col.min = as.numeric(col.min);
      col.max = as.numeric(col.max);

      x = (pmax(val.range[1], pmin(val.range[2], x))-val.range[1])/diff(val.range);
      col.min = pmax(0, pmin(1, col.min))
      col.max = pmax(0, pmin(1, col.max))

      if (invert)
      {
        tmp = col.max
        col.max = col.min
        col.min = tmp
      }

      nna = !is.na(x);

      out = rep(na.col, length(x))
      out[nna] = rgb((col.max[1]-col.min[1])*x[nna] + col.min[1],
                     (col.max[2]-col.min[2])*x[nna] + col.min[2],
                     (col.max[3]-col.min[3])*x[nna] + col.min[3])

      return(out)
}


#' @name lighten
#' @title lighten
#' @description
#' lighten
#'
#' lightens / darkens colors by brighness factor f (in -255 .. 255) that will make lighter if > 0 and darker < 0
#' @author Marcin Imielinski
#' @keywords internal
lighten = function(col, f)
{
  M = col2rgb(col)
  return(apply(matrix(pmax(0, pmin(255, M + f*matrix(rep(1, length(M)), nrow = nrow(M)))), ncol = length(col))/255, 2, function(x) rgb(x[1], x[2], x[3])))
}


#' @name plot.blank
#' @title plot.blank
#' @description
#' Make a blank plot
#'
#' Shortcut for making blank plot with no axes
#' @author Marcin Imielinski
#' @keywords internal
plot.blank = function(xlim = c(0, 1), ylim = c(0,1), xlab = "", ylab = "", axes = F, bg.col = "white", ...)
{
  par(bg = bg.col)
  plot(0, type = "n", axes = axes, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
  #    par(usr = c(xlim, ylim))
}


#' standardize_segs
#'
#' (data frame seg function)
#'
#' Takes and returns segs data frame standardized to a single format (ie $chr, $pos1, $pos2)
#'
#' if chr = TRUE will ensure "chr" prefix is added to chromossome(if does not exist)#'
#' @export
standardize_segs = function(seg, chr = FALSE)
{
  #if (inherits(seg, 'IRangesList'))
  #  seg = irl2gr(seg);

  if (is(seg, 'matrix'))
    seg = as.data.frame(seg, stringsAsFactors = FALSE)

  # if (inherits(seg, 'RangedData') | inherits(seg, 'GRanges') | inherits(seg, 'IRanges'))
  # {
  #   val = as.data.frame(values(seg));
  #   values(seg) = NULL;
  #   seg = as.data.frame(seg, row.names = NULL);  ## returns compressed iranges list
  #   seg$seqnames = as.character(seg$seqnames)
  # }
  # else
  val = NULL;

  field.aliases = list(
    ID = c('id', 'patient', 'Sample'),
    chr = c('seqnames', 'chrom', 'Chromosome', "contig", "seqnames", "seqname", "space", 'chr', 'Seqnames'),
    pos1 = c('start', 'loc.start', 'begin', 'Start', 'start', 'Start.bp', 'Start_position', 'pos', 'pos1', 'left', 's1'),
    pos2 =  c('end', 'loc.end', 'End', 'end', "stop", 'End.bp', 'End_position', 'pos2', 'right', 'e1'),
    strand = c('strand', 'str', 'strand', 'Strand', 'Str')
  );

  if (is.null(val))
    val = seg[, setdiff(names(seg), unlist(field.aliases))]

  seg = seg[, intersect(names(seg), unlist(field.aliases))]

  for (field in setdiff(names(field.aliases), names(seg)))
    if (!(field %in% names(seg)))
      names(seg)[names(seg) %in% field.aliases[[field]]] = field;

  if (chr)
    if (!is.null(seg$chr))
      if (!grepl('chr', seg$chr[1]))
        seg$chr = paste('chr', seg$chr, sep = "");

  if (is.null(seg$pos2))
    seg$pos2 = seg$pos1;

  missing.fields = setdiff(names(field.aliases), c(names(seg), c('chr', 'ID', 'strand')));

  if (length(missing.fields)>0)
    warning(sprintf('seg file format problem, missing an alias for the following fields:\n\t%s',
                    paste(sapply(missing.fields, function(x) paste(x, '(can also be', paste(field.aliases[[x]], collapse = ', '), ')')), collapse = "\n\t")));

  if (!is.null(val))
    seg = cbind(seg, val)

  return(seg)
}


#' @name qstat
#' @title qstat
#' @description
#'
#' Tabulates SGE cluster usage (qstat()) or if full = TRUE flag given will dump out
#' all running jobs in a data.table
#'
#' @author Marcin Imielinski
#' @export
qstat = function(full = FALSE, numslots = TRUE, resources = full)
    {
      nms = c('jobid','prior','ntckt','name','user','project','department','state','cpu','mem','io','tckts','ovrts','otckt','ftckt','stckt','share','queue','slots')
      cmd = 'qstat -u "*" -ext'
      if (resources)
        cmd = paste(cmd, "-r")
      p = pipe(cmd)
      tmp = readLines(p)
      if (resources) ## parse resources text
      {
        resource.line = grepl('^    \\s+', tmp) ## this is a bit fragile since we are depending on the "resource lines" having more than 4 spaces
        lines = tmp[!resource.line]
        line.num = cumsum(!resource.line)
        rmap = cbind(
          data.table(lnum = line.num[resource.line], has.colon = grepl('\\:', tmp[resource.line])),
          as.data.table(matrix(unlist(lapply(strsplit(tmp[resource.line],
                                                      '(\\:\\s+)'),
                                             function(x) c(x[1], x[length(x)], length(x)==2))),
                               ncol = 3, byrow = TRUE)))
        rmap = rmap[, .(lnum,
                        tag = fill.blanks(ifelse(as.logical(V3) | has.colon, gsub('\\W+', '_', str_trim(V1)), NA)),
                        val = ifelse(!as.logical(V3) & has.colon, '', str_trim(V2)))]
        ## reformat the "=" tags
        rmap[grepl('=', val), ":="(tag = paste(sapply(strsplit(val, '=|\\s+'), '[', 1), sep = '_'), val = paste(sapply(strsplit(val, '=|\\s+'), '[', 2), sep = '_'))]
        rtab = dcast.data.table(rmap, lnum ~ tag, value.var = 'val', fun.aggregate =
                                                                       function(x) paste(x, collapse= ';'))
      }
      else
        lines = tmp

      tab = strsplit(str_trim(lines), '\\s+')
      close(p)
      iix = sapply(tab, length)<=length(nms) & sapply(tab, length)>14
      if (sum(iix)==0)
        return(data.table())
      tab = lapply(tab, function(x) x[1:length(nms)])
      tmp = as.data.table(matrix(unlist(tab[iix]), ncol = length(nms), byrow = TRUE))
      setnames(tmp, nms)

      tmp[, host := queue]
      tmp[, queue := sapply(strsplit(queue, '@'), '[', 1)][is.na(queue), host := '']
      tmp[, slots := as.numeric(slots)]
      tmp[, mem := as.numeric(mem)]

      if (resources)
      {
        tmp.lnum = which(iix)
        tmp = cbind(tmp, rtab[.(tmp.lnum), ])
        tmp$lnum = NULL

        if (!is.null(tmp$Hard_Resources))
        {
          tmpl = lapply(strsplit(tmp$Hard_Resources, ';'), strsplit, '=')
        }
      }

        if (!full)
            {
                states = unique(c('r', 'qw', sort(unique(tmp$state))))
                tmp$state = factor(tmp$state, states)
                if (numslots)
                    melted = tmp[, sum(pmax(1, as.numeric(slots), na.rm = TRUE)), by = list(user, state)]
                else
                    melted = tmp[, length(name), by = list(user, state)]
                p = pipe('whoami')
                whoami = readLines(p)
                close(p)
                out = dcast2(melted, user ~ state, value.var = "V1", fun.aggregate = sum)
                setnames(out, gsub('_V1', '', names(out)))
                jcount = rowSums(as.matrix(out[, -1, with = FALSE]))
                out$user = factor(out$user, unique(c(whoami, out$user[order(-jcount)])))
                setkey(out, user)
                return(out)
            }
        else
            return(tmp)
    }



#' @name sstat
#' @title sstat
#' @description
#'
#' Tabulates Slurm cluster usage (sstat()) or if full = TRUE flag given will dump out
#' all running jobs in a data.table
#'
#' @author Zoran Gajic
#' @export
sstat = function(full = FALSE, numslots = TRUE, resources = T){
    asp = "username,groupname,state,name,jobid,associd"
    if(resources){
        asp = c(asp, "timelimit,timeused,submittime,starttime,endtime,eligibletime,minmemory,numcpus,numnodes,priority,nice,reason,reboot")
    }
    p = pipe(paste('squeue -O', paste(asp, collapse = ',')))
    res = readLines(p)
    close(p)

    ## need to do some fixed width parsing training our parser on the header
    header = res[1]
    res = res[-1]
    bks = c(0, str_locate_all(header, '\\s++')[[1]][, 'end'], nchar(header))
    bkdt = data.table(start = bks[-length(bks)]+1, end = bks[-1])
    nms = strsplit(header, '\\s+')[[1]] %>% tolower    
    tmp = lapply(1:length(nms), function(i) str_trim(substr(res, bkdt[i,start], bkdt[i,end])))
    names(tmp) = nms
    out = do.call(data.table, tmp)
   
    out$state = factor(out$state, unique(c(out$state, 'RUNNING'))) %>% relevel('RUNNING')

    if (!full)
    {
      if (numslots)
        out = dcast.data.table(out[, sum(as.numeric(cpus)), by = .(user, state)], user ~ state, fill = 0, value.var = 'V1')[rev(order(RUNNING)), ]
      else
        out = dcast.data.table(out[, .N, by = .(user, state)], user ~ state, fill = 0, value.var = 'N')[rev(order(RUNNING)), ]
    }

    return(out)
}


#' @name gigs
#' @title gigs
#' @description
#'
#' Takes string representing memory and returns numeric string representing number of GB
#' represented by this string or NA
#'
#' @author Marcin Imielinski
#' @export
gigs = function(x)
{
  x = str_trim(x)
  xn = as.numeric(gsub('[A-Za-z]', '', x))

  multiplier =
    ifelse(grepl('Y(i)?(B)?$', x, ignore.case = TRUE), 1e15,
    ifelse(grepl('Z(i)?(B)?$', x, ignore.case = TRUE), 1e12,
    ifelse(grepl('E(i)?(B)?$', x, ignore.case = TRUE), 1e9,
    ifelse(grepl('P(i)?(B)?$', x, ignore.case = TRUE), 1e6,
    ifelse(grepl('T(i)?(B)?$', x, ignore.case = TRUE), 1e3,
    ifelse(grepl('G(i)?(B)?$', x, ignore.case = TRUE), 1,
    ifelse(grepl('M(i)?(B)?$', x, ignore.case = TRUE), 1e-3,
    ifelse(grepl('K(i)?(B)?$', x, ignore.case = TRUE), 1e-6,
    ifelse(grepl('(B)?$', x, ignore.case = TRUE), 1e-9,
           NA)))))))))

  return(xn*multiplier)
}

#' @name seconds
#' @title seconds
#' @description
#'
#' Takes string representing time elapsed and returns numeric string representing number of seconds
#' represented by this string or NA
#'
#' @author Marcin Imielinski
#' @export
seconds = function(x)
{
  naix = !is.na(x)
  out = rep(NA, length(x))
  ll = strsplit(x[naix], "\\:")
  bad = sapply(ll, length)!=4
  naix[naix][bad]= FALSE
  out[naix] = matrix(as.numeric(unlist(ll[!bad])), ncol = 4, byrow = TRUE) %*% cbind(c(3600, 60, 1, 0.01))
  return(out)
}

#' @name qhost
#' @title qstat
#' @description
#'
#' Tabulates per host cluster load
#'
#' @author Marcin Imielinski
#' @export
qhost = function(full = FALSE, numslots = TRUE)
    {
        nms = c('HOSTNAME','ARCH', 'NCPU' ,'NSOC', 'NCOR' ,'NTHR',  'LOAD' , 'MEMTOT'  ,'MEMUSE',  'SWAPTO',  'SWAPUS')
        p = pipe('qhost')
        tab = strsplit(readLines(p), '\\s+')
        close(p)
        iix = sapply(tab, length)<=length(nms) & sapply(tab, length)>6
        if (sum(iix)==0)
            return(data.table())
        tab = lapply(tab, function(x) x[1:length(nms)])
        tmp = as.data.table(matrix(unlist(tab[iix][-c(1:2)]), ncol = length(nms), byrow = TRUE))
        numnms = c('NCPU' ,'NSOC', 'NCOR' ,'LOAD', 'NTHR', 'SWAPTO',  'SWAPUS')
        setnames(tmp, nms)
        for (x in numnms)
            tmp[[x]] = suppressWarnings(as.numeric(tmp[[x]]))
        tmp$MEMUSE = suppressWarnings(pmax(as.numeric(gsub('G', '', tmp$MEMUSE)), 0, na.rm = TRUE))
        tmp$MEMTOT = suppressWarnings(pmax(as.numeric(gsub('G', '', tmp$MEMTOT)), 0, na.rm = TRUE))
        return(tmp)
    }

#' @name qviz
#' @title qviz
#' @description
#'
#' Plots resources (default h_vmem requestes) across cluster either for
#' provided slice of qstat(full = TRUE) output or a fresh call to qstat(full = TRUE)
#'
#' @author Marcin Imielinski
#' @export
qviz = function(res = NULL, queue = NULL, field = "global_mem", frac = FALSE, all = FALSE, plot = TRUE)
{
  if (is.null(res))
    res = qstat(full = TRUE)

  if (!all)
    res = res[state == 'r', ]

  if (!is.null(queue))
  {
    tmp.q = queue
    res = res[res$queue %in% tmp.q, ]
  }

  res$val = res[[field]]

  if (is.character(res$val))
    res$val = gigs(res$val)

  summ = res[, .(val = sum(val, na.rm = TRUE)), keyby = .(host, user)]

  if (frac)
    summ[, val := round(val / sum(val), 2), by = host]

  tmp = dcast.data.table(summ, host ~ user, value.var = "val", fill = 0)

  mat = as.matrix(tmp[,-1])
  rownames(mat) = tmp$host

  if (plot)
    d3heatmap::d3heatmap(mat, scale = 'none', Rowv = TRUE, labRow = rownames(mat), cexCol = 1, cexRow = 0.8)
  else
    mat
}


#' @name ddd
#' @title ddd
#' @description
#'
#' shortcut to gr2dt
#'
#' @export
ddd = function(x)
    {
        if (is.data.frame(x))
            as.data.table(x)
        else
            gr2dt(x)
    }

#' @name relib
#' @title relib
#' @description
#'
#' Reload library
#'
#'
#' @export
relib = function(lib = 'Flow')
    {
        if (sprintf("package:%s", lib) %in% search())
        {
                txt = sprintf("unload('%s')", lib)
                eval(parse(text = txt))

                txt = sprintf("muffle(detach('package:%s', force = TRUE))", lib)
                eval(parse(text = txt))
            }
        txt = sprintf("library(%s)", lib)
        eval(parse(text = txt))
    }



#' @name chron
#' @title chron
#' @description
#'
#' Repeat a command periodically, e.g. every 10 seconds
#'
#' @export
chron = function(expr, period = 5)
    {
        while (TRUE)
            {
                print(eval(expr))
                Sys.sleep(period)
            }
    }

#' @name queues
#' @title queues
#' @description
#' Lists all available queues
queues = function()
    {
        p = pipe('qconf -sql ')
        out = readLines(p)
        close(p)
        return(out)
    }

#' @name loud
#' @title loud
#' @description
#'
#' Runs a system command but prints a message with the output
#'
#' @export
loud = function(x)
    {
        message(x)
        system(x)
        cat('')
    }

#' @name rpipe
#' @title read pipe
#'
#' readsLines from pipe and then closes the pipe
rpipe = function(cmd)
{
    p = pipe(cmd)
    out = readLines(p)
    close(p)
    return(out)
}


#' @name camerplot
#' @title cameraplot
#' @description
#' plots the results of CAMERA in limma package
#' @export
#' @param camera.res output of camera from limma, or data.table with fields $name, $P, $Direction, $FDR
#' @param gene.sets gene set input to camera (named list of indices into the voom.res gene expression matrix)
#' @param genes alternate to voom.res and design can just provide a named numeric vector of effect sizes, named by genes
#' @param voom.res output of voom from limma (default NULL)
#' @param design design matrix input to camera (default NULL)
#' @param max.genes max genes to draw in "leading edge" of gene set
#' @param min.corr minimal abs(correlation) value for leading edge definition
#' @param cex.space  label spacing expansion factor (use if labels get too crowded
#' @param cex.slab  set label cex
#' @param cex.glab  gene label cex
#' @param lwd.notch  notch thickness
#' @param text.shift  amount to shift text from notches (>0, <1)
#' @param text.shift  minimal distance between labels
#' @param height.wf height of the topmost correlation waterfall plot
#' @param col.axis axis color character
#' @param col.ramp ramp from lowest to highest expression to phenotype correlation (default blue, red)
#' @importFrom stringr str_trim
#' @author Marcin Imielinski
cameraplot = function(camera.res, gene.sets, voom.res = NULL, design = NULL, genes = NULL, contrast = ncol(design),
    title = 'Camera Gene set notch plot',
    cex.space = 1,
    col.axis = alpha('gray20', 0.8),
    col.ramp = c('blue', 'red'),
    cex.slab = 1,
    cex.glab = 1,
    lwd.notch = 1,
    tick.w = 0.1,
    max.genes = 10,
    text.shift = 0.5,
    height.wf = 0.1,
    min.corr = 0.1,
    min.dist = 10,
    middle = NULL,
    max.gene.sets= 20,
    gtext.shift = 0.2)
    {

        if (height.wf>1 | height.wf<0)
            {
                warning('Waterfall height should be between 0 and 1, defaulting to 0.2')
                height.wf = 0.2
            }

        if (nrow(camera.res)>max.gene.sets)
            {
                warning(sprintf('Entered %s gene sets .. only will plot the %s topmost (change max.genes.sets param value to increase)', nrow(camera.res), max.gene.sets))
                camera.res = camera.res[1:max.gene.sets, ]
            }

      if (is.null(camera.res$name))
        camera.res$name = rownames(camera.res)
      setnames(camera.res, gsub('\\W', '_', names(camera.res), perl = TRUE))

      if (is.null(camera.res$PValue))
        camera.res$PValue = camera.res$p

      if (is.null(camera.res$PValue))
        camera.res$PValue = camera.res$P

      if (is.null(camera.res$FDR))
        camera.res$FDR = camera.res$fdr

      camera.res = as.data.table(camera.res)[order(-PValue), ]

      if (is.null(camera.res$Direction))
        camera.res$Direction = 'Up'

      my.sets = gene.sets[camera.res$name]
      if (!is.null(voom.res))
        my.corr = apply(voom.res$E, 1, cor, y = design[,contrast])
      else if (!is.null(genes))
        my.corr = genes
      else
        stop('Must provide either voom.res / design pair or named numeric vector of genes->effect sizes')

        my.rank = rank(-my.corr, ties.method = 'first')
        my.set.rank = sapply(my.sets, function(x) my.rank[x][!is.na(my.rank[x])])
      cm = function(x, range = c(-0.5, 0.5), middle = 0)
      {
        width = diff(range)
        x = pmin(pmax(range[1], x), range[2])
        xt = ifelse(x<middle, affine.map(x, xlim = c(range[1], middle), ylim = c(0, 0.5)), affine.map(x, xlim = c(middle, range[2]), ylim = c(0.5, 1)))
        out = rgb(colorRamp(col.ramp)(xt), maxColorValue = 256)
      }
        text.space = 0.01*cex.space*length(my.rank)
        gene.xcoord = length(my.rank)*seq(0, 0.5, length.out = max.genes)
        notch.coord = cbind(unlist(my.set.rank), rep(1:length(my.set.rank), sapply(my.set.rank, length)))
        ## genes to label
        gene.list = lapply(1:length(my.set.rank), function(y)
            {
              x = my.set.rank[[y]]
              if (camera.res[y, Direction == 'Up'])
                x = sort(x[my.corr[names(x)]>min.corr])
              else
                x = rev(sort((x[my.corr[names(x)]<(-min.corr)])))
              
              if (length(x)>0)
                return(x[1:min(length(x), max.genes)])
              else
                return(c())
            })
        gene.nm = unlist(lapply(gene.list, names))
        gene.corr = my.corr[gene.nm]
        gene.coord.top = gene.coord.bot = cbind(unlist(gene.list), rep(1:length(gene.list), sapply(gene.list, length))+tick.w/2)
        gene.coord.top[,2] = gene.coord.top[,2] + gtext.shift
        ## spread out labels
        min.dist = min.dist + nchar(gene.nm)*text.space
        for (i in 2:nrow(gene.coord.top))
        {
            if (nrow(gene.coord.bot)>0)
                {
                    if (camera.res[gene.coord.bot[i, 2], Direction] == 'Up') ## left to right
                    {
                        if (gene.coord.top[i,2] == gene.coord.top[i-1,2] & gene.coord.top[i,1]-gene.coord.top[i-1,1]<min.dist[i])
                            gene.coord.top[i,1] = gene.coord.top[i-1,1] + mean(c(min.dist[i-1], min.dist[i]))
                    }
                    else ## right to left
                    {
                        if (nrow(gene.coord.top)>0)
                            if (gene.coord.top[i,2] == gene.coord.top[i-1,2] & gene.coord.top[i-1,1]-gene.coord.top[i,1]<min.dist[i])
                                gene.coord.top[i,1] = gene.coord.top[i-1,1] - mean(c(min.dist[i-1], min.dist[i]))
                    }
                }
        }
                                        #gene.coord.top[,1] = position.labels(gene.coord.bot[,1], groups = gene.coord.top[,2], min.dist.ll = min.dist, min.dist.pl = 0)
        par(mar = 0.5*c(0,20, 5, 5))
        par(xpd = NA)
        rownames(notch.coord) = unlist(lapply(my.set.rank, names))
      graphics::layout(c(1, 2), heights = c(height.wf, 1-height.wf))
      if (!is.null(voom.res))
        {
          my.corr.range = c(-0.5, 0.5)
          if (is.null(middle))
            middle = 0
        }
      else
        {
          my.corr.range = range(my.corr, na.rm = TRUE)
          if (is.null(middle))
            middle = median(my.corr, na.rm = TRUE)
        }
      plot(0, type ="n", xlim = c(0, length(my.rank)), ylim = my.corr.range, ylab = '', xlab = "", axes = F, main = title)
#        plot(0, type ="n", xlim = c(0, length(my.rank)), ylim = c(-0.5,0.5), ylab = '', xlab = "", axes = F, main = title)
      par(mar = 0.5*c(5,20, 0, 5))
      axis(2, at = signif(seq(my.corr.range[1], my.corr.range[2], length.out = 3),2), col.axis = col.axis)
#        axis(2, at = seq(-0.5, 0.5, length.out = 3), labels = as.character(signif(seq(my.corr.range[1], my.corr.range[2], length.out = 3),2)), col.axis = col.axis)
      mtext(side = 2, 'Gene Correlations', line = 3, col = col.axis)
                                        #      my.corr.transformed = (my.corr-my.corr.range[1])/diff(my.corr.range) - 0.5
      my.corr.transformed = my.corr
      lines(my.rank, my.corr.transformed, type = 'h', col = cm(my.corr, my.corr.range, middle))
        plot(0, type ="n", xlim = c(0, length(my.rank)), ylim = c(0,length(my.set.rank)+1), xlab = "", ylab = "", axes = F)
        ## draw set labels
        set.labs = camera.res[, sprintf('%s (Dir = %s, P = %s, FDR = %s)',
            gsub('REACTOME', '', gsub('_', ' ', name)),
            Direction,
            signif(PValue,2), signif(FDR, 2)),
            , by = 1:length(name)][, str_trim(V1)]
        text(length(my.rank)/2, 1:length(my.set.rank) + text.shift, set.labs, adj = c(0.5, 0), cex = cex.slab*0.5, srt = 0)
                                        # draw gene labels
        text(gene.coord.top[,1],
             gene.coord.top[,2]+0.1*gtext.shift,
             gene.nm,
             col = cm(gene.corr, my.corr.range, middle), cex = cex.glab*0.3, adj = c(0.5, 0))
                                        # draw lines linking gene labels to notches
        segments(gene.coord.top[,1],
                 gene.coord.top[,2],
                 gene.coord.bot[,1],
                 gene.coord.bot[,2]+gtext.shift*0,
                 col = alpha(cm(gene.corr, my.corr.range, middle), 0.2), lty = 1)
                                        # draw "axis" backbone of notches
        segments(rep(0, length(my.set.rank)),
                 1:length(my.set.rank),
                 rep(length(my.rank), length(my.set.rank)),
                 1:length(my.set.rank),
                 col = col.axis, lwd = 0.3, lty = 3)

        ## draw notches
        segments(notch.coord[,1],
                 notch.coord[,2]-tick.w/2,
                 notch.coord[,1],
                 notch.coord[,2]+tick.w/2, col = cm(my.corr[rownames(notch.coord)], my.corr.range, middle), lwd = 2*lwd.notch)

        axis(3, pos = nrow(camera.res)+1, at = c(seq(0, length(my.corr), 5000), length(my.corr)), col.axis = col.axis, col = col.axis, lwd = 0.8, cex.axis = 0.8)
        axis(1, pos = 0.5, at = c(seq(0, length(my.corr), 5000), length(my.corr)), col.axis = col.axis, col = col.axis, lwd = 0.8, cex.axis = 0.8)
        mtext("Gene ranks", side=1, cex.lab=1,col= col.axis)
        mtext("Gene ranks", side=3, cex.lab=1,col= col.axis)
    }



## #' @name position.labels
## #' @title  position.labels
## #' @description
## #'
## #' Given 1D / 2D coordinates of "points", returns coordinates of optimal labels to those
## #' points minimizing distance between points and labels while obeying point-label and
## #' point to point constraints.  Uses Rcplex to solve QP
## #'
## #' @importFrom nloptr cobyla
## #' @param x matrix or vector of numeric coordinates
## #' @param min.dist scalar numeric of minimal distance between labels and labels to points
## #' @param min.dist.lp  scalar or length(x) vector numeric of minimal point to point distance
## #' @param min.dist.ll  scalar numeric of minimal label to label distance
## #' @export
## position.labels = function(x, min.dist = 0, min.dist.pl = min.dist, min.dist.ll = min.dist, x0 = x + runif(length(x)), groups = rep(1, length(x)), ftol = 0.01, xtol = 1e-3, maxeval = 2000)
##     {
##         #library('nloptr')
##         x = cbind(as.numeric(x))
##         if (length(min.dist.pl)==1)
##             min.dist.pl = rep(min.dist.pl, nrow(x))

##         hin = function(y)
##             {
##                 ##for all k and N(k) x[k] - x[N(k)] >= min.dist.ll[k]
##                 dists1 = do.call(c, lapply(split(y, groups), function(yy)
##                     {
##                         ## distance matrix minus diagonal
##                         if (length(yy)==1)
##                             return(min.dist.pl+100)
##                         tmp = as.matrix(dist(yy)) + diag(rep(NA, length(yy)))
##                         tmp[!is.na(tmp)]
##                     }))-min.dist.ll
##                 dists2 = sqrt(rowSums((x-cbind(y))^2))- min.dist.pl
## #                message('max dist ', signif(min(dists1),0), ' ' , signif(min(dists2),0))
##                 return(c(dists1, dists2))
##             }

##         fin = function(y)
##             {
##                 d = sum(sqrt(rowSums((x-cbind(y))^2)))
## #                if (is.na(d))
## #                    browser()
## #               message(d)
##                 return(d)
##             }

##         res = suppressWarnings(slsqp(x0, fn = fin, hin = hin, control = nl.opts(list(xtol_rel = xtol, ftol_abs = ftol, maxeval = maxeval))))
##         return(res$par)
##     }



#' @name parsesnpeff
#' @title parsesnpeff
#'
#' @description
#' parses vcf file containing SnpEff annotations on Strelka calls
#'
#' @param vcf path to vcf
#' @param id id of case
#' @return GRanges object of all variants and annotations
#' @author Kevin Hadi
#' @export
########
parsesnpeff = function(vcf, id = NULL)
           {
            print(vcf)
            fn = c('allele', 'annotation', 'impact', 'gene', 'gene_id', 'feature_type', 'feature_id', 'transcript_type', 'rank', 'variant.c', 'variant.p', 'cdna_pos', 'cds_pos', 'protein_pos', 'distance')
            out = read_vcf(vcf)
            ##            out$ALT = sapply(out$ALT, as.character)
            out$ALT = as.character(unstrsplit(vcf$ALT))
            out$REF = sapply(out$REF, as.character)
            out$vartype = ifelseong(nchar(out$REF) == nchar(out$ALT), 'SNV',
                ifelse(nchar(out$REF) < nchar(out$ALT), 'INS', 'DEL'))
            tmp = lapply(out$ANN, function(y) do.call(rbind, strsplit(y, '\\|'))[, 1:15, drop = FALSE])
            tmpix = rep(1:length(out), sapply(tmp, nrow))
            meta = as.data.frame(do.call(rbind, tmp))
            colnames(meta) = fn
            meta$varid = tmpix
            meta$file = vcf
            meta$pair = id
            out2 = out[tmpix]
            rownames(meta) = NULL
            values(out2) = cbind(values(out2), meta)
            names(out2) = NULL
            out2$ANN = NULL
            vcf$modifier = !grepl('(HIGH)|(LOW)|(MODERATE)', vcf$eff)
            return(out2)
        }



summary.df = function(x, nm = '', last = FALSE)
{
    if (is.null(x))
        out = data.frame(name = nm, method = as.character(NA), p = as.numeric(NA), estimate = as.numeric(NA), ci.lower = as.numeric(NA),  ci.upper = as.numeric(NA), effect = as.character(NA))
    else if ('lm' %in% class(x))
    {
        coef = as.data.frame(summary(x)$coefficients)
        if (nchar(nm)==0)
            nm = rownames(coef)

        colnames(coef) = c('estimate', 'se', 'stat', 'p')
        if (last)
            coef = coef[nrow(coef), ]
        coef$ci.lower = coef$estimate - 1.96*coef$se
        coef$ci.upper = coef$estimate + 1.96*coef$se
        if (summary(x)$family$link %in% c('log', 'logit'))
        {
                    coef$estimate = exp(coef$estimate)
                    coef$ci.upper= exp(coef$ci.upper)
                    coef$ci.lower= exp(coef$ci.lower)
                }
            out = data.frame(name = nm, method = summary(x)$family$family, p = signif(coef$p, 3), estimate = coef$estimate, ci.lower = coef$ci.lower, ci.upper = coef$ci.upper, effect = paste(signif(coef$estimate, 3), ' [',  signif(coef$ci.lower,3),'-', signif(coef$ci.upper, 3), ']', sep = ''))
        }
    else
        out = data.frame(name = nm, method = x$method, p = signif(x$p.value, 3), estimate = x$estimate, ci.lower = x$conf.int[1], ci.upper = x$conf.int[2], effect = paste(signif(x$estimate, 3), ' [',  signif(x$conf.int[1],3),'-', signif(x$conf.int[2], 3), ']', sep = ''))

    out$effect = as.character(out$effect)
    out$name = as.character(out$name)
    out$method = as.character(out$method)
    rownames(out) = NULL
    return(out)
}


#' @name sc.context
#' @title sc.context
#' @description
#'
#' Computes strand collapsed k-nucleotide contexts and representative string for various mutations represented as granges object
#' corresponding to reference coordinate that is being mutated (e.g 0 width for insertion, >=1 width for a del, and 1 width for a SNV, >1 width for a MNV)
#' with alt fields $ALT representing alternate sequences [ACGT]* where any non ACGT sequence is treated
#' as a blank (i.e. for a del) using reference genome hg that is either a BSGenome, ffTrack, or 2bit (i.e. an input into read_seq
#'
#' A strand collapsed k-nucleotide context involves a base that is being altered and (k-1)/2 nucleotides around it, where
#' k is an odd positive integer.
#'
#' Examples:
#'
#' A[T>A]G + represents a T>A mutation happening on the positive strand with an A
#' in 5' position and G in the 3' position
#'
#' A[>ATTTT]G - represents a >ATTTT insertion on the negative strand with an A
#' in 5' position and G in the 3' position
#'
#' A[ATTTT>]G + represents a ATTTT> deletion on the positive strand with an A
#' in 5' position and G in the 3' position
#'
#' default is k = 3
#'
#' @author Marcin Imielinski
#' @export
sc.context = function(mut, hg, k = 3, alt.field = 'ALT', mc.cores = 1, verbose = FALSE)
{
    if (!is(mut, 'GRanges'))
        stop(sprintf('Input argument mut must be granges with field $%s', alt.field))

    if (!(alt.field %in% names(values(mut))))
        stop(sprintf('Input argument mut must be granges with field $%s', alt.field))

    if (k %% 2 != 1)
        stop('k must be an odd number eg 1,3,5 etc')

    ## check if alt field represents actual DNA
    ALT = as.character(DNAStringSet(values(mut)[, alt.field]))

        ## apply strand collapse to tri nuc mutation context
    strand(mut) = '+'
    mut$ref = as.character(get_seq(hg, mut, mc.cores = mc.cores, verbose = verbose)) ## flank gets upstream 5' sequence
    mut$context5 = as.character(get_seq(hg, flank(mut, (k-1)/2), mc.cores = mc.cores, verbose = verbose)) ## flank gets upstream 5' sequence
    mut$context3 = as.character(get_seq(hg, gr.flipstrand(flank(gr.flipstrand(mut), (k-1)/2)), mc.cores = mc.cores, verbose = verbose)) ## flip of flank of flip gets downstream 3' sequence
    mut$alt = as.character(ALT)

    nuc5 = nuc3 = as.character(DNAStringSet(mkAllStrings(DNA_BASES,(k-1)/2)))
    ref = unique(c(as.character(DNAStringSet(mkAllStrings(DNA_BASES, 1))), mut$ref)) ## supplement ref with any "observed refs"
    alt = unique(c(as.character(DNAStringSet(mkAllStrings(DNA_BASES, 1))), as.character(ALT))) ## supplement ref with any "observed alt"
    dict = as.data.table(expand.grid(context5 = nuc5, ref = ref, alt = alt, context3 = nuc3))
    dict = dict[order(factor(ref, levels = c('A', 'C', 'G', 'T'))), ]
    dict[, num := 1:length(context5)]
    setkeyv(dict, c('context5', 'ref', 'alt', 'context3'))
    dictr = dict[, .(context5 = as.character(reverseComplement(DNAStringSet(context3))),
                      ref = as.character(reverseComplement(DNAStringSet(ref))),
                      alt = as.character(reverseComplement(DNAStringSet(alt))),
                     context3 = as.character(reverseComplement(DNAStringSet(context5))))]
    dict$rnum = dict[dictr, num]

    dict[, context := paste(context5, ref, context3, sep = '')]
    dict[, sign := 0]
    dict[is.na(rnum), sign := ifelse(as.character(context)<as.character(reverseComplement(DNAStringSet(context))), 1, -1)] ## if not in dictionary (ie indel) then use the byte rep to arbitrarily but reproducibly call one event mutant
    dict[!is.na(rnum), sign := ifelse(num<rnum, 1, -1)]

    tmp = dict[list(mut$context5, mut$ref, mut$alt, mut$context3), list(context5, ref, alt, context3, sign)]
    tmp[is.na(sign), sign := 0]

    strand(mut) = ifelse(is.na(tmp$sign), '*', ifelse(tmp$sign>0, '+', '-'))
    mut$context5 = as.character(tmp$context5)
    mut$context3 = as.character(tmp$context3)
    mut$alt = as.character(tmp$alt)
    mut$ref = as.character(tmp$ref)

    ## strand flip for negative strand mutations
    if (any(fix <- tmp$sign<0))
    {
        mut$context5[fix] = as.character(reverseComplement(DNAStringSet(as.character(tmp$context3[fix])))) ## keep track of reverse complement contexts
        mut$context3[fix] = as.character(reverseComplement(DNAStringSet(as.character(tmp$context5[fix]))))
        mut$ref[fix] = as.character(reverseComplement(DNAStringSet(as.character(tmp$ref[fix]))))
        mut$alt[fix] = as.character(reverseComplement(DNAStringSet(as.character(tmp$alt[fix]))) )
    }

    ## fill out reciprocal contexts
    mut$context5r = as.character(reverseComplement(DNAStringSet(as.character(tmp$context3)))) ## keep track of reverse complement contexts
    mut$context3r = as.character(reverseComplement(DNAStringSet(as.character(tmp$context5))))
    mut$altr = as.character(reverseComplement(DNAStringSet(as.character(mut$alt))))
    mut$refr = as.character(reverseComplement(DNAStringSet(as.character(mut$ref))) )

    ## put together into signature
    mut$ref.context = paste(mut$context5, mut$ref, mut$context3, sep = '')
    mut$ref.contextr = paste(mut$context5r, mut$refr, mut$context3r, sep = '')
    mut$mut.context = paste(mut$context5, '(', mut$ref, '>', mut$alt, ')', mut$context3, sep = '')
    mut$mut.contextr = paste(mut$context5r, '(', mut$refr, '>', mut$altr, ')', mut$context3r, sep = '')
    mut
}



#' @name staveRDS
#' @title staveRDS
#' @description
#'
#' Stamps and saves RDS file .. i.e. saving datestamped filename and
#' and soft link to the datestamped file
#'
#' @export
staveRDS = function(object, file, note = NULL, ..., verbose = FALSE)
{
  stamped.file = gsub('.rds$', paste('.', timestamp(), '.rds', sep = ''), file, ignore.case = TRUE)
  saveRDS(object, stamped.file, ...)

  if (file.exists(file))
  {
    if (verbose)
      message('Removing existing ', file)
    system(paste('rm', file))
  }

  if (verbose)
    message('Symlinking ', file, ' to ', stamped.file)

  system(paste('ln -sfn', normalizePath(stamped.file), file))

  if (!is.null(note))
  {
    writeLines(note, paste0(stamped.file, '.readme'))
  }
}

#' @name stavePDF
#' @title stavePDF
#' @description
#'
#' Stamps and saves PDF file .. i.e. saving datestamped filename and
#' and soft link to the datestamped file
#'
#' @export
stavePDF = function(object, file, note = NULL, ..., verbose = FALSE)
{
  stamped.file = gsub('.pdf$', paste('.', timestamp(), '.pdf', sep = ''), file, ignore.case = TRUE)
  ppdf(object, stamped.file, ...)

  if (file.exists(file))
  {
    if (verbose)
      message('Removing existing ', file)
    system(paste('rm', file))
  }

  if (verbose)
    message('Symlinking ', file, ' to ', stamped.file)

  system(paste('ln -sfn', normalizePath(stamped.file), file))

  if (!is.null(note))
  {
    writeLines(note, paste0(stamped.file, '.readme'))
  }
}

#' @name label.runs
#' @title label.runs
#' @description
#'
#' For logical input labels all instances of "TRUE" with a unique label and everything else as false
#'
#' For non-logical (e.g. character) input labels, labels each contiguous runs of the same value with a unique label
#' (note: even subsequent runs of an earlier used value in the vector will be given a new unique label)
#' 
#' 
#' @author Marcin Imielinski
#' @export
label.runs = function(x)
{
  if (!is.logical(x))
    {
      cumsum(abs(diff(as.numeric(c(0, as.integer(factor(x))))))>0)
    }
  else ## note will label all runs of FALSE with NA
  {
    as.integer(ifelse(x, cumsum(diff(as.numeric(c(FALSE, x)))>0), NA))
  }
}

#' @name fill.blanks
#' @title fill.blanks
#' @description
#'
#' Takes vector with NAs and "fill in blank" positions i with the value of the last non NA position
#' @export
#' @author Marcin IMielinski
fill.blanks = function(x)
{
  na.ix = is.na(x)

  ## none NA return
  if (!any(na.ix))
    return(x)

  ## all NA return NA
  if (all(na.ix))
    return(x)

  nna.lab = cumsum(!na.ix)
  unna.lab = unique(nna.lab)
  map = structure(match(unna.lab, nna.lab), names = unna.lab)
  x[na.ix] = x[map[nna.lab[na.ix]]]
  return(x)
}



#' @name dodo.call
#' @title dodo.call
#' @description
#' do.call implemented using eval parse for those pesky (e.g. S4) case when do.call does not work
#' @export
dodo.call = function(FUN, args)
{
    if (!is.character(FUN))
        FUN = substitute(FUN)
    cmd = paste(FUN, '(', paste('args[[', 1:length(args), ']]', collapse = ','), ')', sep = '')
    return(eval(parse(text = cmd)))
}


#' @name d3igraph
#' @title d3igraph
#' @description
#' Wrapper around network d3 package to quickly convert igraph to D3 and visualize
#' Note: note can send output to plot.html via vij
#'
#' @export
d3igraph = function(g)
{
    wc <- cluster_walktrap(g)
    members <- membership(wc)
    karate_d3 <- networkD3::igraph_to_networkD3(g, group = members)
    forceNetwork(Links = karate_d3$links, Nodes = karate_d3$nodes,
                 Source = 'source', Target = 'target',
                              NodeID = 'name', Group = 'group')
}





#' @name kill.zombies
#' @title kill.zombies
#' @description
#' Kill R zombies.  Needs to be run from R session that spawned the zombies.
#'
#' @export
kill.zombies = function(x)
{
  includes <- '#include <sys/wait.h>'
  code <- 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
  wait <- inline::cfunction(body=code, includes=includes, convention='.C')
  wait()
  message("Zombie kill complete")
}



#' @name brew
#' @title brew
#' @description
#'
#' Takes factor or character categorical vector and returns same length vector of colors one representing each category
#'
#' @export
brew = function(x, palette = "Accent")
{
  if (!is.factor(x))
    x = factor(x)

  ucols = structure(brewer.master(length(levels(x))), names = levels(x))
  return(ucols[x])
}

#' @name dflm
#' @title dflm
#' @description
#'
#' Formats lm, glm, or fisher.test outputs into readable data.table
#'
#' @export
dflm = function(x, last = FALSE, nm = '')
{
  if (is.null(x))
    out = data.frame(name = nm, method = as.character(NA), p = as.numeric(NA), estimate = as.numeric(NA), ci.lower = as.numeric(NA),  ci.upper = as.numeric(NA), effect = as.character(NA))
  else if (any(c('lm', 'betareg') %in% class(x)))
  {

    coef = as.data.frame(summary(x)$coefficients)
    colnames(coef) = c('estimate', 'se', 'stat', 'p')
    if (last)
      coef = coef[nrow(coef), ]
    coef$ci.lower = coef$estimate - 1.96*coef$se
    coef$ci.upper = coef$estimate + 1.96*coef$se
    if (!is.null(summary(x)$family))
    {
      fam = summary(x)$family$family
        if (summary(x)$family$link %in% c('log', 'logit'))
        {
          coef$estimate = exp(coef$estimate)
          coef$ci.upper= exp(coef$ci.upper)
          coef$ci.lower= exp(coef$ci.lower)
        }
    }
    else
      fam = 'Unknown'

    if (!last)
      nm = paste(nm, rownames(coef))
    out = data.frame(name = nm, method = fam, p = signif(coef$p, 3), estimate = coef$estimate, ci.lower = coef$ci.lower, ci.upper = coef$ci.upper, effect = paste(signif(coef$estimate, 3), ' [',  signif(coef$ci.lower,3),'-', signif(coef$ci.upper, 3), ']', sep = ''))
  }
  else if (class(x) == 'htest')
  {
    if (is.null(x$estimate))
      x$estimate = x$statistic
    if (is.null(x$conf.int))
      x$conf.int = c(NA, NA)
    out = data.table(name = nm, method = x$method, estimate = x$estimate, ci.lower = x$conf.int[1], ci.upper = x$conf.int[2], effect = paste(signif(x$estimate, 3), ' [',  signif(x$conf.int[1],3),'-', signif(x$conf.int[2], 3), ']', sep = ''), p = x$p.value)
  }
  else if (class(x) == 'polr')
  {
    coef = coef(summary(x)) %>% as.data.frame
    nm = paste(nm, rownames(coef))
    coef = as.data.table(coef)
    setnames(coef, c('estimate', 'se', 't'))
    out = data.table(name = nm) %>% cbind(coef)
    out$p =  pnorm(abs(out$t), lower.tail = FALSE) * 2
    out[, ci.lower := estimate-1.96*se]
    out[, ci.upper := estimate+1.96*se]
    out[, effect := paste(signif(estimate, 3), ' [',  signif(ci.lower,3),'-', signif(ci.upper, 3), ']', sep = '')]
  }
  else
  {
    out = data.frame(name = nm, method = x$method, p = signif(x$p.value, 3), estimate = x$estimate, ci.lower = x$conf.int[1], ci.upper = x$conf.int[2], effect = paste(signif(x$estimate, 3), ' [',  signif(x$conf.int[1],3),'-', signif(x$conf.int[2], 3), ']', sep = ''))
  }

  out$effect = as.character(out$effect)
  out$name = as.character(out$name)
  out$method = as.character(out$method)
  rownames(out) = NULL
  return(as.data.table(out))
}



#' @name match.seq
#' @title match.seq
#' @description
#'
#' (Exact) matches a set of character query sequences against a set of (optionally named) subject sequences
#' Returning a GRanges seqnames and coordinates with query.id (and optionally query name) as meta data.
#'
#' @param query character or DNAStringSet
#' @param subject character or DNAStringSet
#' @param mc.cores multithreading for parsing
#' @param ... additional params to PDict
#' @export
match.seq = function(query, subject, mc.cores = 1, verbose = FALSE, ...)
{
  if (is.null(names(subject)))
    names(subject) = 1:length(subject)

  if (any(duplicated(names(subject))))
  {
    warning('Names of subject sequences have duplicates, deduping')
    names(subject) = dedup(names(subject))
  }

  if (!is(query, 'PDict'))
    pdict = Biostrings::PDict(query, ...)

  if (!is(subject, 'DNAStringSet'))
    subject = Biostrings::DNAStringSet(subject)

  res = mclapply(1:length(subj), function(i)
  {
    if (verbose)
      message('Processed subject ', i)
    ll = lapply(Biostrings::matchPDict(pdict, subject[[i]]), as.data.table)
    out = rbindlist(lapply(1:length(ll), function(j) ll[[j]][, query.id := j][, seqnames := i]))
  },  mc.cores = mc.cores)

  if (verbose)
    message('Consolidating results into data.table and GRanges')

  res = rbindlist(res)
  res = res[order(seqnames, query.id, start), ]
  res[, seqnames := names(subject)[seqnames]]

  return(dt2gr(res, seqlengths = seqlengths(subject)))
}


#' @name grok_vcf
#' @title grok_vcf
#' @description
#'
#' Does additional processing of annotated vcf output and produces
#' a more readable granges output.
#'
#' @param x GRanges input
#' @export
grok_vcf = function(x, label = NA, keep.modifier = TRUE, long = FALSE, oneliner = FALSE, verbose = FALSE)
{
  fn = c('allele', 'annotation', 'impact', 'gene', 'gene_id', 'feature_type', 'feature_id', 'transcript_type', 'rank', 'variant.c', 'variant.p', 'cdna_pos', 'cds_pos', 'protein_pos', 'distance')

  if (is.character(x))
    {
      out = suppressWarnings(skidb::read_vcf(x))
      if (is.na(label))
        label = x    }
  else
    out = x

  if (is.na(label))
    label = ''

  if (verbose)
    message('Grokking vcf ', label)

  if (!long)
  {
        vcf = out
        if (length(vcf)>0)
        {
          if (!is.null(vcf$ANN))
          {
            vcf$eff = unstrsplit(vcf$ANN)
            vcf$modifier = !grepl('(HIGH)|(LOW)|(MODERATE)', vcf$eff)
            if (!keep.modifier)
              vcf = vcf[!vcf$modifier]
          }
          vcf$ref = as.character(vcf$REF)
          vcf$alt = as.character(unstrsplit(vcf$ALT))
          vcf = vcf[, sapply(values(vcf), class) %in% c('factor', 'numeric', 'integer', 'logical', 'character')]
          vcf$var.id = 1:length(vcf)
          vcf$type = ifelse(nchar(vcf$ref)==nchar(vcf$alt), 'SNV',
                     ifelse(nchar(vcf$ref)<nchar(vcf$alt),
                            'INS', 'DEL'))
          vcf$label = label
        }
        return(vcf)
  }
  else if (length(out)>0)
  {
        out$REF = as.character(out$REF)
        out$ALT = as.character(unstrsplit(out$ALT))
        out$vartype = ifelse(nchar(out$REF) == nchar(out$ALT), 'SNV',
                      ifelse(nchar(out$REF) < nchar(out$ALT), 'INS', 'DEL'))
        if (is.null(out$ANN))
          stop('no $ANN column, check to see if annotated VCF is formatted in the SnpEff style')
        else
          out$eff = unstrsplit(out$ANN)

        out$modifier = !grepl('(HIGH)|(LOW)|(MODERATE)', out$eff)
        if (!keep.modifier)
          out = out[!out$modifier]
        if (inherits(out$ANN, 'character'))
          annlist = strsplit(out$ANN, ',')
        else
          annlist = out$ANN %>% as.list
        tmp = lapply(annlist, function(y) do.call(rbind, lapply(strsplit(y, '\\|'), '[', 1:15)))
        tmpix = rep(1:length(out), elementNROWS(tmp))
        meta = as.data.frame(do.call(rbind, tmp))
        colnames(meta) = fn
        meta$varid = tmpix
        meta$file = label
        out2 = out[tmpix]
        rownames(meta) = NULL
        values(out2) = cbind(values(out2), meta)
        names(out2) = NULL
        out2$ANN = NULL
        precedence = c('trunc', 'cnadel', 'cnadup', 'complexsv', 'splice', 'inframe_indel', 'fusion', 'missense', 'promoter', 'regulatory', 'noncoding', 'inv', 'synonymous', '')
        eff = readRDS(system.file('extdata', 'snpeff_ontology.rds', package = 'skitools'))[, short := factor(short, precedence)][!is.na(short), ]

        .short = function(vcf)
        {
          tmp = strsplit(as.character(vcf$annotation), '\\&')
          dtl = data.table(eff = unlist(tmp), id = rep(1:length(tmp), lengths(tmp)))  %>% merge(eff, by = 'eff', allow.cartesian = TRUE) %>% unique(by = 'id')
          setkey(dtl, id)
          vcf$short = dtl[.(1:length(vcf)), short]
          return(vcf)
        }
        
        out2 = .short(out2)

        if (oneliner)
          out2$oneliner = paste(
            ifelse(!is.na(out2$gene),
                   as.character(out2$gene),
                   as.character(out2$annotation)),
            ifelse(nchar(as.character(out2$variant.p))>0,
                   as.character(out2$variant.p),
                   as.character(out2$variant.c)))
    }
    return(out2)
}


#' @name grok_bcf
#' @rdname
#' @title Reads and parses bcf via bcftools call
#' @param bcf path to bcf file
#' @param gr optional granges to query
#' @param bpath path to bcftools binary executable
#' @export
grok_bcf = function(bcf, gr = NULL, bpath = "/nfs/sw/bcftools/bcftools-1.1/bcftools", label = NA, filter = 'PASS', snv = FALSE, indel = FALSE, het = FALSE, hom = FALSE, keep.modifier = TRUE, long = FALSE, oneliner = FALSE, verbose = FALSE)
{
  cmd = sprintf('%s view %s', bpath, bcf)

  if (is.na(label))
    label = bcf

  if (!is.null(gr))
  {
    wins = paste(gr.string(gr.stripstrand(gr)), collapse = ',')
    cmd = paste(cmd, '-r', wins)
  }

  if (!is.null(filter) && !is.na(filter))
  {
    cmd = paste(cmd, sprintf('-i \'FILTER="%s"\'', filter))
  }

  if (het)
  {
    cmd = paste(cmd, '-g het')
  }

  if (indel)
  {
    cmd = paste(cmd, '-v indels')
  }

  if (snv)
  {
    cmd = paste(cmd, '-v snps')
  }

  if (het)
  {
    cmd = paste(cmd, '-g het')
  }

  if (hom)
  {
    cmd = paste(cmd, '-g hom')
  }

  ## quick dunlist
  .dunlist = function(x)
  {
    if (is.null(names(x)))
      names(x) = 1:length(x)
    out = data.table(listid = rep(names(x), elementNROWS(x)), V1 = unlist(x))
    setkey(out, listid)
    return(out)
  }
  
  if (verbose)
    message(cmd)
    
  p = pipe(cmd)
  lines = readLines(p)
  close(p)

  is.header = grepl('^\\#', lines)
  header = lines[is.header]
  contigs = strsplit(gsub('^\\#\\#', '', grep('contig', header, value = TRUE)), ',')
  sl = structure(names = gsub('.*ID=\\<', '', sapply(contigs, '[', 1)),
                 as.numeric(gsub('>$', '', gsub('.*length=\\<', '', sapply(contigs, '[', 2)))))
  
  other = lines[!is.header]
  if (length(other))
    {
      out = fread(paste(c(other, ''), collapse = '\n'), sep = "\t", header = FALSE)
      sn = unlist(strsplit(gsub('^\\#', '', header[length(header)]), '\\t'))
      sfields = sn[-c(1:9)]
      setnames(out, sn)
      out[, seqnames := as.character(CHROM)]
      out[, start := POS]
      out[, end := POS]
      ## unpack bcf "format" + sample fields
      fdat = .dunlist(strsplit(out$FORMAT, ':'))
      setnames(fdat,2,'field')
      out$FORMAT = NULL
      for (sfield in sfields) ## can be more than one sample field
      {
        fdat$value = .dunlist(strsplit(out[[sfield]], ':'))$V1
        fdatc = dcast.data.table(copy(fdat)[, field := paste(sfield, field, sep = '_')], listid ~ field, value.var = 'value')
        out[, listid := as.character(1:.N)]
        out = merge(out, fdatc, by = 'listid', all.x = TRUE, allow.cartesian = TRUE)
      }
      
      ## unpack "info" field
      idat = .dunlist(strsplit(out$INFO, ';'))
      idat = cbind(idat, colsplit(idat$V1, pattern = "=", names = c("field","value")))
      idatc = dcast.data.table(idat, listid ~ field, value.var = 'value')
      out$INFO = NULL
      out[, listid := as.character(1:.N)]
      mcols = setdiff(names(idatc), c('REF', 'ALT'))
      out = merge(out, idatc[, mcols, with = FALSE], by = 'listid', all.x = TRUE, allow.cartesian = TRUE)
      out = dt2gr(out, seqlengths = sl)
      out = grok_vcf(out, keep.modifier = keep.modifier, long = long, oneliner = oneliner, verbose = verbose, label = label)
    }
  else
    out = GRanges(seqlengths = sl)
  return(out)
}



#' multicoco
#'
#' multi-scale coverage correction
#'
#' Given gc and mappability coverage correction at k "nested" scales finds the coverage
#' assignment at the finest scale that yields the best correction at every scale
#'
#' FUN = is a function that takes in a data frame / granges with
#' $reads and other covariate functions $gc, $map and outputs a vector of corrected read counts
#'
#' cov is a constant with GRanges of coverage samples with (by default) fields $reads, $map, $gc
#'
#' base = is the multiple with which to perform "numlevs" additional scales of correction
#'
#####################################
multicoco = function(cov,
    numlevs = 1, ## numbers of scales at which to correct
    base = 100, ## scale multipler
    fields = c("gc", "map"), # fields of gc which to use as covariates
    iterative = TRUE,
    presegment = TRUE, ## whether to presegment
    min.segwidth = 5e6, ## when presegmenting, min seg width
    mono = FALSE, ## just do single iteration at finest scale
    verbose = T,
    imageroot = NULL, ## optional file root to which to dump images of correction
    FUN = NULL, ## function with which to correct coverage (by default loess correction modified from HMMcopy that takes in granges with fields $reads and other fields specified in "fields"
    ..., ## additional args to FUN
    mc.cores = 1)
    {
        require(Rcplex)
        if (verbose)
           cat('Converting to data.table\n')

        WID = max(width(cov))

        cov.dt = as.data.table(as.data.frame(cov))

        sl = structure(as.numeric(1:length(seqlevels(cov))), names = seqlevels(cov))

        if (verbose)
            cat('Grouping intervals\n')

        ## compute level means
        ## lev 0 = raw data
        ## lev 1 = base-fold collapsed
        ## lev 2 = base^2-fold collapsed
        ## etc
        parentmap= list() ## data.tables that map lev labels at level k  to parent lev labels
        cov.dt[, lev0:=as.character(1:length(seqnames))]
        for (k in 1:numlevs)
            {
                if (verbose)
                    cat('Defining', base^k, 'fold collapsed ranges\n')
                cov.dt[, eval(paste("lev", k, sep = '')) := as.character(sl[seqnames] + as.numeric(Rle(as.character(1:length(start)), rep(base^k, length(start)))[1:length(start)])/length(start)), by = seqnames]
                parentmap[[k]] = data.table(parent = cov.dt[, get(paste("lev", k, sep = ''))], child = cov.dt[, get(paste("lev", k-1, sep = ''))], key = 'child')[!duplicated(child), ]
            }

        if (presegment) ## perform rough segmentation at highest level
            {
                sl = seqlengths(cov)
                if (verbose)
                    cat('Presegmenting at ', as.integer(WID*base^(numlevs)), ' bp scale \n')
                require(DNAcopy)
                tmp.cov = seg2gr(cov.dt[,list(chr = seqnames[1], start = min(start), end = max(end), strand = strand[1], reads = mean(reads, na.rm = T)), by = get(paste("lev", numlevs, sep = ''))], seqlengths = sl)
                ix = which(!is.na(values(tmp.cov)[, 'reads']))
                cna = CNA(log(values(tmp.cov)[, 'reads'])[ix], as.character(seqnames(tmp.cov))[ix], start(tmp.cov)[ix], data.type = 'logratio')
                tmp = print(segment(smooth.CNA(cna), alpha = 1e-5, verbose = T))
                tmp = tmp[!is.na(tmp$loc.start) & !is.na(tmp$chrom) & !is.na(tmp$loc.end), ]
                seg = sort(seg2gr(tmp, seqlengths = sl))
                seg = seg[width(seg)>min.segwidth] ## remove small segments
                seg.dt = as.data.table(as.data.frame(seg))
                seg = seg2gr(seg.dt[, list(seqnames = seqnames,
                    start = ifelse(c(FALSE, seqnames[-length(seqnames)]==seqnames[-1]), c(1, start[-1]), 1),
                    end = ifelse(c(seqnames[-length(seqnames)]==seqnames[-1], FALSE), c(start[-1]-1, Inf), seqlengths(seg)[seqnames]))], seqlengths = sl)
                seg = gr.val(seg, tmp.cov, 'reads') ## populate with mean coverage
                seg$reads = seg$reads/sum(as.numeric(seg$reads*width(seg))/sum(as.numeric(width(seg)))) ## normalize segs by weigthed mean (so these become a correction factor)
            }
        else
            seg = NULL

        if (verbose)
            cat('Aggregating coverage within levels \n')

        ## list of data frames showing scales of increasing collapse

        cov.dt[, ix := 1:nrow(cov.dt)]

        cmd1 = paste('list(ix.start = ix[1], ix.end = ix[length(ix)], reads = mean(reads, na.rm = T),', paste(sapply(fields, function(f) sprintf("%s = mean(%s, na.rm = T)", f, f)), collapse = ','), ')', sep = '')

        cmd2 = paste('list(lab = lev0, reads,', paste(fields, collapse = ','), ', seqnames, start, end)', sep = '')

        if (mono)
            {
                if (verbose)
                    cat('Mono scale correction \n')
                 grs = list(cov.dt[, eval(parse(text=cmd2))])
                 numlevs = 1
             }
        else
            {
                grs = c( list(cov.dt[, eval(parse(text=cmd2))]),
                    lapply(1:numlevs, function(x)
                        {
                            if (verbose)
                                cat('Aggregating coverage in level', x,  '\n')
                            out = cov.dt[, eval(parse(text=cmd1)), keyby = list(lab = get(paste('lev', x, sep = '')))]
                            out[, ":="(seqnames = cov.dt$seqnames[ix.start], end = cov.dt$end[ix.start], start = cov.dt$start[ix.start])]
                            out[, ":="(ix.start= NULL, ix.end = NULL)]
                            return(out)
                        }))
            }

        setkey(grs[[1]], 'lab')

        ## modified from HMMCopy to
        ## (1) take arbitrary set of covariates, specified by fields vector
        ## (2) employ as input an optional preliminary (coarse) segmentation with which to adjust signal immediately prior to loess
        ## NOTE: (this only impacts the loess fitting, does not impose any segmentation on thed ata)
        ##
        if (is.null(FUN))
            FUN = function(x, fields = fields, samplesize = 5e4, seg = NULL, ## seg is a Granges with meta data field $reads
                verbose = T, doutlier = 0.001, routlier = 0.01)
                {
                    if (!all(fields %in% names(x)))
                        stop(paste('Missing columns:', paste(fields[!(fields %in% names(x))], collapse = ',')))

                    x$valid <- TRUE
                    for (f in fields)
                        {
                            x$valid[is.na(x[, f])] = FALSE
                            x$valid[which(is.infinite(x[, f]))] = FALSE
                        }

                    if (verbose)
                        cat('Quantile filtering response and covariates\n')

                    range <- quantile(x$reads[x$valid], prob = c(routlier, 1 - routlier), na.rm = TRUE)

                    if (verbose)
                        cat(sprintf("Response min quantile: %s max quantile: %s\n", round(range[1],2), round(range[2],2)))

                    domains = lapply(fields, function(f) quantile(x[x$valid, f], prob = c(doutlier, 1 - doutlier), na.rm = TRUE))
                    names(domains) = fields

                    x$ideal <- x$valid
                    x$ideal[x$reads<=range[1] | x$reads>range[2]] = FALSE

                    for (f in fields)
                        x$ideal[x[, f] < domains[[f]][1] | x[, f] > domains[[f]][2]] = FALSE

                    if (verbose)
                        cat(sprintf('Nominated %s of %s data points for loess fitting\n', sum(x$ideal), nrow(x)))

                    set <- which(x$ideal)

                    if (length(set)<=10)
                        {
                            warning("Not enough samples for loess fitting - check to see if missing or truncated data?")
                            return(x$reads)
                        }

                    for (f in fields)
                        {
                            if (verbose)
                                message(sprintf("Correcting for %s bias...", f))

                            select <- sample(set, min(length(set), samplesize))

                            x2 = x[, c('reads', f)]
                            x2$covariate = x2[, f]

                            x2s = x2[select, ]

                            if (!is.null(seg)) ## here we apply a prelmiinary segmentation to correct for large scale copy variation allow more power to reveal the covariate signal
                                {
                                    if (verbose)
                                        message('Applying preliminary segmentation prior to loess fitting')

                                    x.grs = gr.val(seg2gr(x[select, ], seqlengths = NULL), seg, 'reads')
                                    x2s$reads = x2s$reads/x.grs$reads
                                }

                            fit = loess(reads ~ covariate, data = x2s, span = 0.3)
                            if (is.na(fit$s))
                                {
                                    warning("Using all points since initial loess failed")
                                    fit = loess(reads ~ covariate, data = x2[select, ], span = 1)
                                }

                            tryCatch(
                                {
                                    if (!is.na(fit$s))
                                        {
                                            domain = domains[[f]]

                                            yrange <- quantile(x2s$reads, prob = c(routlier, 1 - routlier), na.rm = TRUE)
                                            df = data.frame(covariate = seq(domain[1], domain[2], 0.001))

                                            if (!is.null(imageroot))
                                                {
                                                    out.pdf = paste(imageroot, ifelse(grepl("/$", imageroot), '', '.'), f,'_correction.pdf', sep = '')
                                                    if (verbose) {
                                                        cat("Dumping figure to", out.pdf, "\n")
                                                    }

                                                    pdf(out.pdf, height = 6, width = 6)
                                                    plot(x2s$covariate, x2s$reads, col = alpha('black', 0.1), pch = 19, cex = 0.4, xlim = domain, ylim = yrange, ylab = sprintf('signal before %s correction', f), xlab = f);
                                                    lines(df$covariate, predict(fit, df), col = 'red', lwd = 2)
                                                    dev.off()
                                                }
                                            x$reads = x2$reads/predict(fit, x2) ## apply correction
                                        }
                                    else
                                        print("Loess failed, yielding NA loess object, continuing without transforming data")
                                }, error = function(e) print("Unspecified loess or figure output error"))
                        }
                    return(x$reads)
                }

        if (verbose)
            cat('Correcting coverage at individual scales\n')

        ## level 1,2, ..., k corrections
        ## these are the computed corrected values that will be input into the objective function

        if (iterative) ## iteratively correct
            {
                correction = NULL
                for (i in rev(1:length(grs)))
                    {
                        cat('Correcting coverage at ', WID*base^(i-1), 'bp scale, with', nrow(grs[[i]]), 'intervals\n')
                        if (i != length(grs))
                            grs[[i]]$reads = grs[[i]]$reads/correction[parentmap[[i]][grs[[i]]$lab, parent], cor]

                        if (WID*base^(i-1) > 1e5) ## for very large intervals do not quantile trim, only remove NA
                            grs[[i]]$reads.corrected = FUN(as.data.frame(grs[[i]]), fields, doutlier = 0, seg = seg)
                        else
                            grs[[i]]$reads.corrected = FUN(as.data.frame(grs[[i]]), fields, seg = seg);


                        if (is.null(correction))
                            correction = data.table(lab = grs[[i]]$lab, cor = grs[[i]]$reads / grs[[i]]$reads.corrected, key = 'lab')
                        else
                            {
                                ## multiply new correction and old correction
                                old.cor = correction[parentmap[[i]][grs[[i]]$lab, parent], cor]
                                new.cor = grs[[i]]$reads / grs[[i]]$reads.corrected
                                correction = data.table(lab = grs[[i]]$lab,  cor = old.cor * new.cor, key = 'lab') ## relabel with new nodes
                            }
                    }

                cov.dt$reads.corrected = grs[[1]][cov.dt$lev0, ]$reads.corrected
                rm(grs)
                gc()
#                cov.dt[, reads.corrected := grs[[1]][lev0, reads.corrected]]
            }
        else ## parallel mode
            {

                ## compute marginal values / ratios for reads and covariates
                if (verbose)
                    cat('Computing marginal values of read coverage and covariates\n')

                ## now compute marginal ratio relative next-level mean for all levels except for top level
                for (k in 1:numlevs)
                    {
                        if (verbose)
                            cat('Defining marginal coverage for', base^k, 'fold collapsed ranges\n')
                        ix = parentmap[[k]][grs[[k]]$lab, parent]
                        grs[[k]]$reads = grs[[k]]$reads / grs[[k+1]][ix, reads]
                        print('fucken bitch yeah yeah')

                        for (f in fields)
                            grs[[k]][, eval(parse(text = sprintf("%s := grs[[k]]$%s / grs[[k+1]][ix, %s]",f, f, f)))]
                    }

                grs = mclapply(grs, function(x) {
                    if (verbose)
                        cat('Correcting coverage at ', base^(k-1), 'fold collapse, with', nrow(grs[[k]]), 'intervals\n')
                    x$cor = FUN(as.data.frame(x), fields);
                    return(x)
                }, mc.cores = mc.cores)

                gc()

                for (k in 1:length(grs))
                    {
                        cov.dt[, eval(paste('cor', k-1, sep = '')) := as.numeric(NA)]
                        cov.dt[, eval(paste('cor', k-1, sep = ''))] = grs[[k]][cov.dt[, get(paste('lev', k-1, sep = ''))], ]$cor
                    }

                ulev = unique(cov.dt[, get(paste('lev', numlevs, sep = ''))])

                cov.dt$lid = factor(cov.dt[, get(paste('lev', numlevs, sep = ''))])
                cov.dt[, gid := 1:length(start)]
                cov.dt[, eval(parse(text = paste("reads.corrected := ", paste('cor', 0:numlevs, sep = '', collapse = "*"))))]

                ##        out = rbindlist(mclapply(ulev, function(x) .optimize_subtree(cov.dt[get(paste('lev', numlevs,sep = '')) == x], numlevs), mc.cores = mc.cores))
                ## .optimize_subtree = function(this.chunk, numlevs)   {
                ##     browser()
                ##     this.chunk.og = this.chunk
                ##     this.chunk.og[, reads.corrected := NA]
                ##     this.chunk.og$id = as.character(1:nrow(this.chunk.og))
                ##     setkey(this.chunk.og, 'id')
                ##     this.chunk = this.chunk.og[!is.na(cor0), ]

                ##     if (verbose)
                ##         cat('chunk', as.integer(this.chunk.og$lid)[1], 'of', length(levels(this.chunk.og$lid)), '\n')


                ##     if (nrow(this.chunk)>0)
                ##         {
                ##             y_lev = rbindlist(
                ##                 c(lapply(0:(numlevs-1), function(x) this.chunk[, list(lev = x, cor = get(paste('cor', x, sep = ''))[1], parent_lab = get(paste('lev', x+1, sep = ''))[1]), by = list(lab = get(paste('lev', x, sep = '')))]), list(this.chunk[, list(lev = numlevs, cor = get(paste('cor', numlevs, sep = ''))[1], parent_lab = NA), by = list(lab = get(paste('lev', numlevs, sep = '')))])))

                ##                         #                    y_lev = y_lev[!is.na(cor), ]

                ##             ## build tree structure of rows of y_lev i.e. the variables in our optimization
                ##             ## map nodes to their parents
                ##             y_lev[, parent.ix := match(paste(lev+1, parent_lab), paste(lev, lab))]
                ##             y_lev[, id := 1:length(parent.ix)]

                ##             if (any(!is.na(y_lev$parent.ix)))
                ##                 {
                ##                     ## make constraints matrix - one constraint per unique parent
                ##                     ##
                ##                     ## this defines parents in terms of their children (mean function)
                ##                     A = y_lev[!is.na(parent.ix), sparseMatrix(as.integer(factor(c(parent.ix, parent.ix))), c(parent.ix, id), x = c(rep(-1, length(parent.ix)), rep(1, length(id))))]
                ##                     A = cBind(A, A*0) ## add room for residual variables
                ##                     Arhs = rep(0, nrow(A)) ## right hand side of A is 0

                ##                     ## residual constraints relate nodes, their parents and the fits contained in the "cor" columns
                ##                     ## if q_ki is the fit for the ith node in the kth level
                ##                     ## then x_ki - q_ki*x_p(ki) + r_ki = 0

                ##                     ## this encodes basic residual matrix whose rows have the form:  x_ki - q_ki*x_p(ki) = residual
                ##                     ## except for the top level, which is a single parentless node, and has the form:  x_ki - q_ki = residual
                ##                     ##
                ##                     ## the variables are indexed 1:nrow(y_lev), and the residuals are the next nrow(y_lev) indices
                ##                     R = sparseMatrix(rep(1:nrow(y_lev),2), c(1:nrow(y_lev), nrow(y_lev) + 1:nrow(y_lev)), x = rep(c(1,-1), each = nrow(y_lev)))
                ##                     rownames(R) = as.character(y_lev$id)
                ##                     R[y_lev[!is.na(parent.ix), cbind(id, parent.ix)]] = -y_lev[!is.na(parent.ix), cor]
                ##                     Rrhs = y_lev[, ifelse(is.na(parent.ix), cor, 0)]
                ##                     R = R[!is.na(y_lev$cor), ]
                ##                     Rrhs = Rrhs[!is.na(y_lev$cor)]

                ##                     ## make objective function
                ##                     ##
                ##                     ## problem:
                ##                     ## find x, r minimizing ||r||
                ##                     ##
                ##                     ## Ax = 0 (encodes node-parent relationships where parent = mean(children)
                ##                     ## x_ki - q_ki*x_p(ki) - r_ki = 0 (for k<numlevs)
                ##                     ## x_ki - r_ki = q_ki (for k = numlevs, a single node)
                ##                     ##
                ##                     y_lev[, obj.weight := 1/length(id), by = lev]

                ##                     Q = diag(c(y_lev$obj.weight*0, 1 + 0*y_lev$obj.weight)) ## only put 1 weights on the matrix entries corresponding to the residual

                ##                     vtype = rep('C', ncol(A)) ## all variables are continuous
                ##                     sense = rep('E', nrow(A) + nrow(R)) ## all constraints specify equality
                ##                     tilim = 10;
                ##                     cvec = rep(0, ncol(A)) ## all variables are continuous

                ##                     ## now we need to encode the sum relationships between x0 etc.
                ##                     sol = Rcplex(cvec = cvec, Amat = rBind(A, R), bvec = c(Arhs, Rrhs), sense = sense, Qmat = Q, lb = rep(c(0, -Inf), each = nrow(y_lev)), ub = Inf, n = 1, objsense = "min", vtype = vtype, control = list(tilim = tilim))
                ##                     y_lev$opt = sol$xopt[1:nrow(y_lev)]
                ##                     setkey(this.chunk, 'lev0')
                ##                     this.chunk[y_lev[lev==0, lab], reads.corrected := y_lev[lev == 0, ]$opt]
                ##                     this.chunk[, reads.corrected := y_lev[lev == 0, ]$opt]
                ##                     this.chunk.og$reads.corrected = this.chunk[this.chunk.og$id, ]$reads.corrected
                ##                 }
                ##         }
                ##     return(this.chunk.og)
                ## }
            }

        if (verbose)
            cat('Converting to GRanges\n')

        gc()

        out = seg2gr(as.data.frame(cov.dt), seqlengths = seqlengths(cov))

        if (verbose)
            cat('Made GRanges\n')

        gc()
        return(out)
    }



#' @name print_eq
#' @title print_eq
#' @description
#' Prints and formats equations Ax = b in matrix Ab whose last column is the vector B
#'
#' @export
print_eq = function(Ab, xlabels = NULL)
  {
    Ab = as.matrix(Ab)

    A = Ab[, -ncol(Ab), drop = F]
    b = Ab[, ncol(Ab)]

    if (!is.null(xlabels))
      colnames(A) = xlabels
    else if (is.null(colnames(A)))
      colnames(A) = paste('x', 1:ncol(A), sep = '')
    else if (length(ix <- which(nchar(colnames(A)) == 0)) != 0)
      colnames(A)[ix] = paste('x', ix, sep = '')

    cat('', paste(sapply(1:nrow(A), function(x,y)
                        {
                          ix = which(A[x, ]!=0)
                          if (length(ix)>0)
                            paste(signif(as.numeric(A[x, ix]),2), y[ix], sep = ' ', collapse = ' +\t')
                          else
                            '0'
                        },
                        colnames(A)), '\t=', b, '\n'))
  }



#' @name gr2json
#' @title gr2json
#'
#' @description
#'
#' Dumps GRanges into JSON with metadata features as data points in  "intervals"
#'
#'
#'
#' @param GRange input granges object
#' @param file output json file
#' @author Marcin Imielinski
gr2json = function(intervals, file, y = rep("null", length(intervals)), labels = '', maxcn = 100, maxweight = 100)
{

    #' ++ = RL
    #' +- = RR
    #' -+ = LL
    qw = function(x) paste0('"', x, '"')

    ymin = 0;
    ymax = maxcn;

    nodes = intervals
    id = rep(1:length(nodes), 2)

    node.dt = data.table(
        iid = 1:length(nodes),
        chromosome = qw(as.character(seqnames(nodes))),
        startPoint = as.character(start(nodes)),
        strand = as.character(strand(nodes)),
        endPoint = as.character(end(nodes)),
        y = y,
        title = labels)

    oth.cols = setdiff(names(values(nodes)), colnames(node.dt))
    node.dt = as.data.table(cbind(node.dt, values(nodes)[, oth.cols]))

    oth.cols = union('type', oth.cols)
    if (is.null(node.dt$type))
        node.dt$type = 'interval'

    intervals.json = node.dt[, paste0(
        c("intervals: [", paste(
                              "\t{",
                              "iid: ", iid,
                              ", chromosome: ", chromosome,
                              ", startPoint: ", startPoint,
                              ", endPoint: ", endPoint,
                              ", y: ", y,
                              ", title: ", qw(title),
                               ", strand: ", qw(strand),
                              eval(parse(text = ## yes R code making R code making JSON .. sorry .. adding additional columns
                                             paste0("paste0(",
                                                    paste0('", ', oth.cols, ':", qw(', oth.cols, ')', collapse = ','),
                                                    ")", collapse = ''))),
                              "}",
                              sep = "",
                              collapse = ',\n'),
          "]"),
        collapse = '\n')
        ]

    meta.json =
        paste('meta: {\n\t',
              paste(
                  c(paste('"ymin:"', ymin),
                  paste('"ymax:"', ymax)),
                  collapse = ',\n\t'),
              '\n}')

    out = paste(c("var json = {",
                  paste(
                      c(meta.json,
                        intervals.json
                        ),
                      collapse = ',\n'
                  ),"}"),
                sep = "")

    writeLines(out, file)
}



##################################
#' segment
#'
#' Wrapper around cumSeg to segment numeric data in an input GRanges with signal meta data field (e.g. $signal)
#' Returns a GRAnges of piecewise constant regions with their associated value
#'
##################################
cumseg = function(gr, field = 'signal', log = T, type = 'bic', alg = 'stepwise', S = 1, verbose = F, mc.cores = 1, ...)
  {
    require(cumSeg)

    if (!(field %in% names(values(gr))))
      stop(sprintf('Field "%s" not a meta data columnin the input GRanges', field))

    if (log)
      values(gr)[, field] = log(values(gr)[, field])

    good.ix = !(is.infinite(values(gr)[, field]) | values(gr)[, field]==0)
    good.ix[is.na(values(gr)[, field])] = F
    gr = gr[good.ix]

    if (length(gr) == 0)
      return(gr[c(), field])

    args = list(...)

    if (!("sel.control" %in% names(args)))
      args$sel.control = sel.control(type = 'bic', alg = alg, S = S)

    grl = split(gr, as.character(seqnames(gr)))

    out = do.call('c', mclapply(names(grl), function(chr)
      {
        if (verbose)
          cat('Starting ', chr, '\n')
        this.gr = grl[[chr]]

        if (length(this.gr)<=1)
          return(si2gr(this.gr)[chr])

        tmp = do.call(jumpoints, c(list(values(this.gr)[, field], start(this.gr)), args))
        out = GRanges(chr, IRanges(c(1, tmp$psi), c(tmp$psi, seqlengths(this.gr)[[chr]])), seqlengths = seqlengths(this.gr))
        values(out)[, field] = NA
        values(out)[, field] = tmp$est.means

        if (verbose)
          cat('\t... generated ', length(out), ' segs\n')

        return(out)
      }, mc.cores = mc.cores))

    if (log)
      values(out)[, field] = exp(values(out)[, field])

    return(out)
  }



############################
#' pp2gb
#'
#' converts purity / ploidy to gamma / beta (or reverse)
#'
#' takes in gr with signal field "field"
#'
#' @param purity value between 0 and 1
#' @param ploidy value nonnegative
#' @param mu vector of n segment averages
#' @param w vector of n segment widths
#' @param gamma non-negative value
#' @param beta non-negative value
#' @return
#' list with purity / ploidy / gamma / beta entries
#' @export
############################
pp2gb = function(purity = NA, ploidy = NA, mu = NA, w = NA, gamma = NA, beta = NA)
{
    if (all(is.na(mu)) & all(is.na(w)))
        stop('mu and w should be non-empty')

    if (length(mu) != length(w))
        stop('mu and w should match up in length')

    w = as.numeric(w)
    mu[is.infinite(mu)] = NA
    w[is.na(mu)] = NA
    sw = sum(w, na.rm = T)
    mutl = sum(mu * w, na.rm = T)
    ncn = rep(2, length(mu))
    ploidy_normal = sum(w * ncn, na.rm = T) / sw  ## this will be = 2 if ncn is trivially 2

    if (is.na(gamma) & is.na(beta) & !is.na(purity) & !is.na(ploidy))
    {
        gamma = 2*(1-purity)/purity
        beta = ((1-purity)*ploidy_normal + purity*ploidy) * sw / (purity * mutl)
    }
    else if (!is.na(gamma) & !is.na(beta) & is.na(purity) & is.na(ploidy))
    {
        purity = 2/(gamma+2)
        v = beta * mu - ncn * gamma / 2
        ploidy = sum(v*w, na.rm = TRUE)/sum(w, na.rm = TRUE)
    }
    else
        stop('Either gamma and beta are empty OR purity and ploidy are empty')
    return(list(purity = purity, ploidy = ploidy, gamma = gamma, beta = beta))
}





#' @name karyoSim
#' @title karyoSim
#' @description
#'
#' Simulate (random) evolution of rearrangements according to input junctions, which are provided as a GRangesList, and
#' grouped into "events" by events list (list of numeric vectors or of lists of numeric vectors indexing "junctions")
#'
#' Goal of the simulation is to instantiate a collection of junctions (+/- approach some copy number profile)
#' through a sequence of rearrangements and whole-chromosome copy changes
#'
#' Junctions are sequences of signed reference intervals that are contiguous on the
#' on the tumor genome (usually pairs)
#'
#' Each event consists of either
#' (1) a "quasi reciprocal sequence" (QRS) of junctions, implemented during a single "cell cycle", and are specified by vectors of junction indices,
##    where no indices are repeated (save the last and first, which specifies a cycle)
#' (2) a set of sets of junctions, specified as a list of list of junction indices, again without repetition, corresponding to complex events
#'     spanning multiple QRS's or "cell cycles" e.g. a BFB, which involve a replication step in between each QRS.  The subsequent QRS's
#'     (attempted to be) instantiated in cis to the last item in the previous QRS
#' (3) a GRanges object specifying a reference locus and meta data field $type = "loss" or "gain" specifying one or more pieces of reference genome that should
#'     be lost or gained at a given step.
#'
#' - Events are interpreted as strings of one or more "quasi reciprocal sequence" (QRS) of junctions
#'   which may be closed / cyclic (if they begin and end with the same junction index) or open. in which case they will result
#'   in at least some interval loss.  We restrict QRS's to contain at most one repeated junction, and this has to be the first
#'   and the last item in the sequence.  "Quasi" reciprocal means that we allow some sequence to be lost or gained in between breaks.
#' - Every QRS is instantiated in the current genome, by mapping junctions, which are specified in haploid reference
#'   coordinates to intervals on the current genome.  By default, instantiation is chosen so that the source interval of
#'   every junction in the QRS is on the same chromosome as the target interval on the previous junction in the QRS.  If this is the
#'   case, then we say that the current junction  "follows" the previous one in this WQRS instantiation.  Quasi-reciprocality is then applied
#'   by possibly adding intervals at the site of a break (i.e. if the target interval of the previous junction is upstream of the
#'   source interval of the next junction).   In situations where an instance of a subsequent junction cannot be found to followthe current
#'   junction, then the chain is either (1) interpreted as "broken", i.e. equiv to an unbalanced rearrangement (if strict.chain = F or (2)
#'   the event is discarded (if strict.chain = T)
#' - Junctions / events can have many possible instantiations at a given round of evolution.
#'   This is because a given haploid interval on the reference can be associated with many loci on the tumor chromosome
#'   (in the simplest case, two homologues of the same chromosome)
#'   By default, the following preferences are exercised for choosing junction instantiations:
#'   (1) if a chromosome strand can be found that contains all the intervals in the junction (2) a chromosome whose both strands
#'   contain all the intervals in the junction (3) a set of chromosome that instantiates the event as a chain of junctions
#'   these prefereences can be over-ridden by specifying instant.local and instant.chain flags
#' - After every cycle we do a "clean up" which involves (1) rejoining any pairs of broken ends that were partners at the previous
#'   iteration (2) removing any fragments that lack a telomere (if req.tel = T) or lack other req.gr (3) replacing reverse complements
#'   of chromosomes from previous iteration that were rearranged in the previous iteration with the reverse complements of their alteration
#'   products in the current iteration.
#' - Every junction is implemented <exactly once> during the evolutionary history, i.e. lightning does not strike twice, infinite
#'   sites model
#'
#' p.chrom = prob of chrom event at each simulation step
#' p.chromloss = probability of chromosomal loss | chrom event (default 0.5)
#' p.chromgain = probability of chromosomal gain | chrom event (default 0.5)
#' lambda.chrom = poisson lambda of number of different chromosomes gained or lost at a chromosomal event
#' lambda.chromgain = poisson lambda of number of chromosomes gained at each "gain" event (default lambda.chrom)
#' labmda.chromloss = poisson lambda of number of chromosomes lostd at each "loss" event (default lambda.chrom)
#'
#' p.wgd = prob of whole genome doubling at each simulation step
#'
#' Optionally can provide a copy profile cn.profile (GRanges tiling genome with $cn meta data field) and heuristic will be applied
#' to attempt to "evolve" the simulation towards the observed copy profile (to be implemented)
#'
#' Output is provided as
#' - (if full = F) list with fields
#'                 $chroms = Named GRangesList of final tumor chromosomes in terms of reference intervals
#'                 $gChain = gChain mapping reference to tumor genome
#'                 $cn = gRanges in reference genome of copy counts of reference intervals
#'                 $events = data frame of event indices with field $id (for event id), $desc (see below for description)
#' - (if full = T) List of lists, each item k corresponding to each stage k of evolution and
#'                 containing the following fields:
#'                 $chroms = Named GRangesList of tumor chromosomes at that stage in terms of reference intervals
#'                 $gChain = gChain mapping reference to current genome k
#'                 $gChain.last = gChain mapping last evolution step to current (from reference in first item of history)
#'                 $cn = gRanges in reference genome of copy counts of reference intervals
#'                 $event = list with $id, $desc that gives the id and description of event
#'                          for chromosomal loss / gain $id = 'chromgain', or 'chromloss', $desc = indices chromosomes
#'                          for ra event, $id event id, desc = junctions involved
#'
########################################
karyoSim = function(junctions, # GRangesList specifying junctions, NOTE: currently only allows "simple junctions", i.e. locus pairs, eventually will
                               # allow multi (i.e. two or more) range pairs
  events = NULL, # list of integer vectors or lists of integer vetcors corresponding to "events", list item can be GRanges with meta data field $type
                 # with values "loss" or "gain"
  p.chrom = 0, ## probability of chromosomal event at each time step
  p.wgd = 0,  ## conditional prob of wgd given chromosomal event
  p.chromgain = 0.5, ## conditional probability of chromosome gain given not WGD, chromosomal event
  cn = NULL, ## GRanges with 'cn' property
  req.gr = NULL, ## GRanges that every chromosome needs to overlap in order to make it to the next evolution time step
                  ## e.g. centromeres
  req.tel = TRUE, ## logical flag whether to require every chromosome to have telomeres at both ends at every evolution time step
  neo.tel = NULL, ## GRanges specifying intervals that qualify as neo-telomeres, these will only be applied if req.tel = TRUE
  haploid = T, ## tells us whether input genome is haploid, in which case we will begin the simulation with a "genome doubling"
  local.junction = T, ## if T, we prefer to instantiate intervals of a junction "locally" on the same chromosome,
                     # if F we allow all instantiations to be equally likely
  local.qrs = T, ## if T, we prefer QRS instantiations that operate on a single chromosome, if F then all are equally good
                 ## by "prefer", we mean that we score each instantiation, and then choose the best scoring (or "a best", if there
                 ## are ties)
  force.event = T, ## if T, will attempt to implement QRS / event even if QRS can only be instantiated partially or in fragments
  lambda.chrom = 0, lambda.chromgain = lambda.chrom, lambda.chromloss = lambda.chrom,
  full = F, ## full output?
  random.event = T, ## if not random.event, then provided event order will be followed, and blank events will trigger a (random) chromosomal event
  precedence = NULL, ## length(events) x length(events) binary matrix of DAG entries ij representing whether event i occurs before event j
  dist = 1000, ## distance at which to allow deletion breaks
  verbose = T,
  ... # other optional input to chromoplexy()
  )
  {
    kag = JaBbA:::karyograph(junctions)

    ## check events to make sure kosher

    if (!all(sapply(1:length(events), function(x) {
      if (is(events[[x]], 'GRanges'))
        {
          ev = events[[x]]
          if (is.null(values(ev)$type))
            stop('GRanges specifying copy number events must have $type field set to "gain" or "loss"')
          else if (any(!(values(ev)$type %in% c('gain', 'loss'))))
            stop('GRanges specifying copy number events must have $type field set to "gain" or "loss"')
          else
            T
        }
      else if (is(events[[x]], 'vector'))
        {
          ev = as.numeric(events[[x]])
          if (any(!(abs(ev) %in% 1:length(junctions))))
            stop(sprintf('Event %s has junction index out of bounds', x))
          else
            T
        }
      else if (is.list(events[[x]]))
        all(sapply(1:length(events[[x]], function(y)
               {
                 ev = events[[x]][[y]]
                 if (any(!abs(ev)) %in% 1:length(junctions))
                   stop(sprintf('Event %s, subevent %s has junction index out of bounds', x, y))
                 else
                   T
               })))
    })))
      stop('Some events are of the wrong type')

    ## junctions in terms of graph nodes
    junctions.kg = kag$ab.edges[, c(1:2), ]

    if (is.null(events))
      {
        pc = chromoplexy(kag, dist = dist, ...)
        events = c(pc$paths, lapply(pc$cycles, function(x) c(x, x[1])))
        singletons = setdiff(1:nrow(kag$ab.edges), unlist(events)) ## all events that are not part of a path or cycle
        if (length(singletons)>0)
          events = c(events, split(singletons, 1:length(singletons)))
      }

    ## this helps us keep track of how many junctions we have accounted for
    ## while choosing events
    cn.event = sapply(events, is, 'GRanges')
    event2junction = sparseMatrix(i = 1, j = 1, x = 0, dims = c(length(events), nrow(junctions.kg)))
    if (any(!cn.event))
        {
          ix = cbind(unlist(mapply(function(x,y) rep(x, y), which(!cn.event), sapply(events[!cn.event], length))), unlist(events[!cn.event]))
          ix = ix[!duplicated(ix), ]
          event2junction[ix] = 1
        }

    ## events are "done" if we have already used a junction / ra that belongs to that event
    done.events = rep(F, nrow(event2junction))
    done.junctions = rep(F, ncol(event2junction))
    k = 0 ## evolution time step

    #### some local utility functions
    ####
    ####

    ## makes list mapping reference signed intervals to chromosomal interval coordinates
    .rid2cid = function(intervals)
      {
        out = split(c(1:nrow(intervals), -(1:nrow(intervals))), c(intervals[, '+'], intervals[, '-']))
        out = out[1:max(c(intervals[, '+'], intervals[, '-']))]
        return(out)
      }

    ## updates state with chrom gain or loss
    .chrom_change = function(state, gain = NULL, loss = NULL)
      {
        if (is.null(gain))
          gain = c()

        if (is.null(loss))
          loss = c()

        keep = !(state$intervals[, 'i'] %in% loss)
        gain = state$intervals[, 'i'] %in% gain;

        gain.intervals = state$intervals[gain, , drop = F]
        gain.intervals[, 'i'] = gain.intervals[, 'i'] + 0.01 ## give these new intervals a unique new chrom name

        ## conatenate and rename
        tmp.intervals = rbind(state$intervals[keep, ], gain.intervals)
        tmp.intervals = .munlist(lapply(split(1:nrow(tmp.intervals), tmp.intervals[,'i']), function(x) tmp.intervals[x, 3:ncol(tmp.intervals), drop = F]))

        out = list(
          intervals = tmp.intervals,
          rid2cid = .rid2cid(tmp.intervals),
          cid2prev = c(which(keep), which(gain))
          )

        return(out)
      }

    ## unlists and cbinds matrices (if dim = 2) or rbinds matrices (if dim = 1)
    ## whose first column specifies the list item index of the entry
    ## and second column specifies the sublist item index of the entry
    ## and the remaining columns specifies the value(s) of the vector
    ## or matrices.
    .munlist = function(x, dim = 1)
      {
        if (dim == 2)
          return(t(rbind(i = unlist(lapply(1:length(x), function(y) rep(y, ncol(x[[y]])))),
                         j = unlist(lapply(1:length(x), function(y) 1:ncol(x[[y]]))),
                         do.call('cbind', x))))
        else
          return(cbind(i = unlist(lapply(1:length(x), function(y) rep(y, nrow(x[[y]])))),
                       j = unlist(lapply(1:length(x), function(y) 1:nrow(x[[y]]))),
                       do.call('rbind', x)))
      }


    ## takes k vectors of length n_1 , ... , n_k and outputs a matrix
    ## of dimension (n_1 x ... x n_k) x k representing their cartesian product
    .cartprod = function(...)
      {
        vecs = list(...)
        if (length(vecs)==0)
          return(NULL)
        out = matrix(vecs[[1]], ncol = 1)
        if (length(vecs)==1)
          return(out)
        if (length(vecs)>1)
          for (i in 2:length(vecs))
            {
              y = vecs[[i]]
              ix = cbind(rep(1:nrow(out), length(y)), rep(1:length(y), each = nrow(out)))
              out = cbind(out[ix[,1], ], y[ix[,2]])
            }
        return(out)
      }

    ## main data structure to keep track of current state of chromosomal evolution
    current.state = list(
      intervals = list(),     ## matrix of n signed reference intervals on k chromosomes of current reference genome
                              ## "i" maps current chromosome, and "j" maps position in that chromosome
                              ## '+' col has rids on positive strand and '-' rids on negative strand of current genome
      rid2cid = list(), ## list mapping reference interval ids to signed current ids
      cid2prev = c()   ## vector mapping current ids (cids) to signed cids of previous genome (prev)
      )

    ix = order(as.character(seqnames(kag$tile)), as.character(strand(kag$tile)))
    ix.pos = ix[which(as.logical( strand(kag$tile)[ix]=='+'))]
    ix.neg = ix[which(as.logical( strand(kag$tile)[ix]=='-'))]
    current.state$intervals = .munlist(mapply(function(x, y) cbind("+"=x, "-"=y),
      split(ix.pos, as.character(seqnames(kag$tile)[ix.pos])),
      split(ix.neg, as.character(seqnames(kag$tile)[ix.neg]))))
    current.state$rid2cid = .rid2cid(current.state$intervals)
    current.state$cid2prev = 1:nrow(current.state$intervals)

    ## map reference intervals to their rev comp
    int2rc = suppressWarnings(match(kag$tile, gr.flipstrand(kag$tile)))

    ## keep track of telomeric reference intervals (todo: specify centromeric intervals as input or other customizable characteristics
    ## that will specify chromosomes that are kept from timepoint to timepoint in the simulation)
    is.tel = kag$tile$is.tel

    if (!is.null(neo.tel))
      is.tel = is.tel || gr.in(kag$tile, gr.stripstrand(neo.tel))

    if (!is.null(req.gr))
      is.req = gr.in(kag$tile, gr.stripstrand(req.gr))

    ## if we are in haploid land, first step is a "whole genome doubling" that will give us homologues
    if (haploid)
      current.state = .chrom_change(current.state, gain = unique(current.state$intervals[, 'i']))


    ## history is a list of current.states
    history = list()

    step = 0;  ## step in evolution

#'    done.events = rowSums(event2junction)==0
    done.events = rep(F, length(events))
    done.this.round = rep(F, length(done.events)) ## also keep track of events that have been done this round

    while (!all(done.events | done.this.round)) ## only finish when we have done all events or at least tried them once in this round
      {
        if (verbose)
          cat(sprintf('Evolution step %s\n', step))

        history = c(history, list(current.state))
        step = step +1

        ## rearrangement event triggered either randomly based on p.chrom
        ## or non-randomly if our event[[step]] is non empty

        if (random.event)
          k = sample(which(!done.events), 1)
        else
          k = step

        force.cn = is(events[[k]], 'GRanges')

        ## if random.event = F and event is empty, then will trigger (random) chromosomal loss / gain (see else statement below)
        if (!force.cn & runif(1)>=p.chrom)
          {
            if (verbose)
              cat(sprintf('trying event %s: %s\n', k, sapply(done.events, function(x) paste('(', x, ')', collapse = ', '), collapse = ', ')))

            done.this.round[k] = T
            done.events = done.events | rowSums(event2junction[, events[[k]], drop = FALSE])!=0
            last.qrs = NULL ## will store last.qrs in current coordinates if we have a multi-qrs event

            this.event = events[[k]]

            if (!is.list(this.event))
              this.event = list(this.event)

            qrs.i = 0
            abort = F
            while (qrs.i < length(this.event) & !abort) ## iterate through qrs's in this event
              {
                if (verbose)
                  cat(sprintf('QRS %s of event %s\n', k, qrs.i))

                qrs.i = qrs.i+1
                this.qrs = this.event[[i]] ## vector of junctions
                is.cycle = this.qrs[length(this.qrs)] == this.qrs[1]

                ## instantiate this QRS
                ## i.e. assign the reference-centric junctions.kg with pairs of intervals in current genome

                ## first enumerate all paths involving instantiations of junctions in qrs
                ## qrs.paths k x m matrix of k QRS instantiations, each consisting of a sequence of
                ## (signed) cids m across n junction id (n < m)
                ##
                ## qrs.juncid maps the columns of qrs paths to junction id in the sequence
                ##
                qrs.paths = array()
                qrs.juncid = c()
                qrs.fid = c() ## keeps track of "fragments" (in case force.event = F, and we are not

                j = 0
                while (j < length(this.qrs) & !abort)
                  {
                    j = j+1

                    if (j == 1)
                      {
                        if (verbose)
                          cat(sprintf('junction  %s of qrs %s event %s in step %s\n', j, qrs.i, k, step))

                        ## take cartesian product of all instantiations of node 1 and node 2 in junction
                        qrs.paths = do.call(.matcart, current.state$rid2cid[junctions.kg[abs(this.qrs[j]), ,ifelse(this.qrs[j]>0, 1, 2)]])
                        qrs.juncid = rep(j, dim(qrs.paths)[3]) ## this keeps
                        qrs.fid = rep(1, dim(qrs.paths)[3])

                        ## if we have multi qrs event and last.qrs is defined, then we constrain
                        ## the first event to be in cis (i.e. on the same chromosome) as the previous
                        if (!is.null(last.qrs))
                          qrs.paths = qrs.paths[qrs.paths[, 1] %in% current.state$intervals[abs(last.qrs), 'i'], ]

                        if (verbose)
                          cat(sprintf('\t dim of qrs.paths: (%s, %s)\n', nrow(qrs.paths), ncol(qrs.paths)))

                        if (nrow(qrs.paths)==0)
                          {
                            abort = T
                            break
                          }
                      }
                    else
                      {
                        if (verbose)
                          cat(sprintf('junction %s of qrs %s event %s in step %s\n', j, qrs.i, k, step))

                        qrs.paths.old = qrs.paths
                        tmp = do.call(.cartprod, c(list(1:nrow(qrs.paths)), current.state$rid2cid[junctions.kg[this.qrs[j],]]))
                        qrs.paths = cbind(qrs.paths[tmp[,1], ], tmp[, 2:ncol(tmp)])

                        if (verbose)
                          cat(sprintf('\t dim of qrs.paths: (%s, %s)\n', nrow(qrs.paths), ncol(qrs.paths)))

                        ## ensure that instantiations of source of current .cid in qrs is in cis with previous (including strand)
                        ## we check if the instantation of the first interval in this cid
                        ## is on the same chromosome as the instantiation of the last interval in the
                        ## last cid.  If not, then we throw them out
                        ## TODO: will check whether instantiations are within some threshold distance of each other
                        ##
                        ## if we run out of instantiations (i.e. the qrs cannot be fully instantiated) then we can
                        ## either implement a fragmented qrs (if force.event = T) or abort and try a different event
                        ##

                        keep = current.state$intervals[abs(qrs.paths[,length(qrs.juncid)+1]), 'i'] ==
                          current.state$intervals[abs(qrs.paths[,length(qrs.juncid)]), 'i'] ## check which have same chromosome
                        keep = keep & sign(qrs.paths[,length(qrs.juncid)+1]) == sign(qrs.paths[,length(qrs.juncid)]) ## check which have same strand

                        qrs.juncid = c(qrs.juncid, rep(j, length(junctions.kg[this.qrs[j],, '+'])))
                        if (!any(keep))
                          {
                            if (!force.event)
                              {
                                if (verbose)
                                  cat(sprintf('Aborting event %s at qrs %s\n', k, qrs.i))
                                abort = T
                                break
                              }
                            else ## if can't keep any, then we keep all and just start a new "fragment"
                              qrs.fid = c(qrs.fid, rep(qrs.fid[length(qrs.fid)]+1, length(junctions.kg[this.qrs[j],, '+'])))
                          }
                        else
                          {
                            qrs.paths = qrs.paths[keep, ]
                            qrs.fid = c(qrs.fid, rep(qrs.fid[length(qrs.fid)], length(junctions.kg[this.qrs[j],, '+'])))
                          }
                      }
                  }

                if (abort)
                  break

                qrs.score = junction.score = rep(0, nrow(qrs.paths))

                if (local.junction | local.qrs)
                  {
                    ## score how many pairs are on same chromosome
                    tmp = do.call('cbind',
                      lapply(split(1:length(qrs.juncid), qrs.juncid),
                             function(x) current.state$intervals[abs(qrs.paths[,x[1]]), 'i'] == current.state$intervals[abs(qrs.paths[,x[length(x)]]), 'i']))

                    if (local.junction)
                      junction.score = junction.score + rowSums(tmp)

                    ## check if entire event is on a single chromosome
                    if (local.qrs)
                      qrs.score = qrs.score + apply(tmp, 1, prod)
                  }

                if (verbose)
                  cat(sprintf('\t final dim of qrs.paths: (%s, %s)\n', nrow(qrs.paths), ncol(qrs.paths)))

                ## sort with respect to score, keep best, sample uniformly from best
                ord.ix = order(qrs.score, junction.score, decreasing = T)
                keep = sample(ord.ix[qrs.score[ord.ix] == qrs.score[ord.ix[1]] & junction.score[ord.ix] == junction.score[ord.ix[1]]], 1)

                ## this.qrs.path is vector of signed cid's
                this.qrs.path = qrs.paths[keep,]

                ## apply junctions

                ## first "check out" the relevant chromosomes, we will replace these in the output genome
                ## these are only chromosomes on which junction instantiations (instantations of intervals
                ## at ends of junctions) lie, and not any internal junction interval
                tmp.ix = which(diff(c(0, qrs.juncid, 0))!=0)
                internal = !c(1:length(qrs.juncid) %in% c(tmp.ix, tmp.ix-1))
                qrs.adj = !(1:length(current.state$intervals) %in% c(tmp.ix, tmp.ix + 1))

                ## make new "fragments" data structure initially representing
                ## all strands of "checked out" chromosomes, which we will break and join

                ## chroms and strands to "check out" of the current genome
                check.out = cbind(chr = current.state$intervals[abs(this.qrs.path[!internal]), 'i'], str = sign(this.qrs.path[!internal]))
                check.out = check.out[!duplicated(check.out), ]

                if (verbose)
                  cat(sprintf('\t checked out chroms: (%s)\n', paste(sign(this.qrs.path)[!internal]*current.state$intervals[abs(this.qrs.path[!internal]), 'i'], collapse = ',')))

                ## cook up list of fragment ids
                tmp.fid = split(cumsum(unlist(lapply(current.state$intervals[check.out[,1]], function(x) rep(1, length(x))))),
                  unlist(lapply(1:nrow(check.out), function(x) rep(x, ncol(current.state$intervals[[check.out[x,1]]])))))

                ## fragments is n x 4 matrix representing n "checked out" single strand DNA intervals across k fragments
                ## and their fragment number (i), fragment pos (j), signed current genome interval id (cid),
                ## flags specifying whether it is an amp bridge (is.bridge) and to be broken (to.break)
                ##
                fragments = .munlist(mapply(function(chr, str) {
                  cid = which(current.state$intervals[, 'i'] == chr)
                  if (str==1)
                    cbind(cid = cid, to.break = F, is.bridge = F)
                  else ## rev prev strands
                    cbind(cid = -cid, to.break = F, is.bridge = F)
                  }, check.out[, 'chr'], to.break = F, check.out[, 'str']), dim = 1)

                if (verbose)
                  cat(sprintf('\t checked out %s fragments comprising %s intervals on %s chromosomes\n', length(unique(fragments[, 'i'])), nrow(fragments)), length(unique(check.out[, 'chr'])))

                ## this k x 2 matrix will keep track of left and right (current genome) neighbors of fragments
                ## "left" and "right" store either NA or the fragment number of the neighbor
                fragment.partners = array(NA, dim = c(unique(fragments[, 'i']), 2), dimnames = list(NULL, c('left', 'right')))

                ## current2frag = n x 2 matrix mapping signed cid --> fid
                ## we also keep mapping of current genome signed intervals to (non bridge) fragment interval,
                ## within a single QRS, this mapping will be one to one, since the only signed intervals
                ## that get duplicated are within amp bridges.
                ##
                ## note: BFB's are implemented by connecting signed intervals to their reciprocal, which will result
                ## in duplication only after we "strand complete" the fragments produced by the QRS,
                ## however, during the implementation of the QRS, there will be a one to one mapping
                ## between current genome signed intervals and single stranded fragments
                ##
                current2frag = matrix(NA, nrow = length(current.state$cid2prev), ncol = 2,
                  dimnames = list(NULL, c('+', '-')))
                current2frag[cbind(abs(fragments[, 'cid']), ifelse(fragments[, 'cid']>0, 1, 2))] = 1:nrow(fragments)
                current2frag.nna = !is.na(current2frag) ## will be useful for updating this

                is.start = c(T, diff(qrs.fid)!=0) ## vector of qrs fragment starts (will only be one T entry if force.event = F)
                qrs.iter = split(1:length(qrs.juncid), qrs.juncid)

                k = 0
                ## iterate through all the adjacencies / junctions in this qrs
                for (k in 1:(length(qrs.iter) - !is.cycle)) ## stop at next to last if is.cycle
                  {
                    ## if first in qrs (or first in qrs fragment), then apply both breaks, and no amp bridge
                    ## if last in qrs cycle, then apply no breaks, and possibly two amp bridges (from previous, to first)
                    ## otherwise apply one break, and possibly one amp bridge (from previous)

                    if (verbose)
                      cat(sprintf('current fragments: \n%s\n',
                                  paste(paste('[', sapply(split(fragments[, 'cid'], fragments[, 'i']), paste, collapse = ' '), ']', sep = ''), collapse = '\n')))

                    fragments[, 'to.break'] = F

                    ## junction.cid = signed cid of instantiations of current junction intervals on genome
                    this.junction = this.qrs.path[qrs.iter[[k]]]

                    ## we locate their fragment locations (flocs) fragment_id fragment_pos
                    this.junction.fids = current2frag[cbind(abs(this.junction), ifelse(this.junction>0, 1, 2))]

                    if (is.start[k])
                      {
                        fragments[this.junction.fids, 'to.break'] = T

                        if (verbose)
                          {
                            frag.ix = fragments[, 'i'] == fragments[this.junction.fids, 'i']
                            frag1 = fragments[frag.ix[fragments[frag.ix, 'j'] <= fragments[this.junction.fids, 'j']], 'cid']
                            frag2 = fragments[frag.ix[fragments[frag.ix, 'j'] > fragments[this.junction.fids, 'j']], 'cid']
                            cat(sprintf('[%s] --> [%s], [%s]\n',
                                        paste(fragments[frag.ix, 'cid'], collapse = ' '),
                                        paste(frag1, collapse = ' '),
                                        paste(frag2, collapse = ' ')))
                          }
                      }
                    else
                      {
                        ## if not start, check for backward amp bridge between target
                        ## adj of previous junction and source adj of current
                        ##
                        ## amp bridge occurs only if last adj targ is upstrand of this edge source on current genome
                        last.adj.targ.qix = qrs.iter[[k-1]][length(qrs.iter[[k-1]])]
                        last.adj.targ = this.qrs.path[, last.adj.targ.qix]
                        this.adj.source = this.qrs.path[, qrs.iter[[k]][1]]

                        ## sanity check: last adj targ and this edge source should be on the same
                        ## chromosome and strand in the current genome
                        if (!(last.adj.targ[1] == this.adj.source[1] & last.adj.targ[3] == this.adj.source[3]))
                          stop('something wrong: last.adj.targ and current.adj are not on same chr and strand')

                        ## amp.bridge only if targ at or before source, which is left on pos strand and right on neg
                        is.amp.bridge = ((last.adj.targ[3] == 1 & last.adj.targ[2] <= this.adj.source[2]) |
                                         (last.adj.targ[3] == 2 & last.adj.targ[2] >= this.adj.source[3]))

                        ## backward amp bridge will be added to fragment opposite last target
                        ## amp bridge consists of intervals between last target and current source (inclusive)
                        if (is.amp.bridge)
                          {
                            ## find interval where to add amp.bridge
                            ## this interval is opposite last target in current genome
                            ## which is left of last.adj.targ for neg strand and right of last.adj.targ for positive
                            this.chr = which(current.state$intervals[, 'i'] == last.adj.targ[1])
                            to.add.fid = current2frag[this.chr[last.adj.targ[2]], last.adj.targ[3]]

                            ## make amp.bridge frag (i j cid is.bridge to.break)
                            amp.bridge.frag = cbind(cid = c(1, -1)[last.adj.targ[3]] * this.chr[last.adj.targ[2]:this.adj.source[2]])

                            amp.bridge.frag = cbind(
                              i = fragments[to.add.fid, 'i'],
                              j = fragments[to.add.fid, 'j'] + 1:nrow(amp.bridge.frag), ## we are adding to the right so j will be the new index
                              amp.bridge.frag, is.bridge = T, to.break = F)

                            if (verbose)
                              cat(sprintf('Making amp bridge %s\n', amp.bridge.frag))

                            ## sanity check: is there something wrong, i.e. the interval that we are adding to
                            ## is not at the right end of its fragment
                            if (to.add.fid != nrow(fragments))
                              if (fragments[to.add.fid, 'i'] == fragments[to.add.fid+1, 'i'])
                                stop('something is wrong: right amp bridge to be added to internal fragment')

                            ## add amp.bridge to right of to.add.fid in fragments
                            aft.ix = c()
                            bef.ix = 1:to.add.fid
                            if (to.add.fid != nrow(fragments))
                              aft.ix = (to.add.fid+1):nrow(fragments)
                            fragments = cbind(fragments[bef.ix, ], amp.bridge, fragments[aft.ix, ])

                            ## update current2frag
                            pix = c(bef.ix, rep(NA, nrow(amp.bridge)), aft.ix)
                            current2frag[current2frag.nna] = match(current2frag[current2frag.nna], pix)

                            ## remove any reference neighbors from the right side of this frag
                            ## (i.e. it can only be attached through a subsequent adjacency)
                            fragment.partners[flocs[to.add.fid,'i'], 'right'] = NA
                          }
                      }

                    if (!(is.cycle & k == (length(qrs.iter)-1))) ## unless next to last in cycle, break target adjacency in junction
                      {
                        tmp.fid = this.junction.fids[length(this.junction.fids)]+1
                        if (fragments[tmp.fid, 'i'] != fragments[tmp.fid-1, 'i'])
                          stop('Something is wrong, we are breaking at the beginning of a fragment')
                        fragments[tmp.fid, ] = T
                      }
                    else
                      {
                        ## if we are next to last iteration of cycle QRS,
                        ## check for forward amp bridge between target
                        ## adj of this junction and source adj of next (i.e. first)
                        ##
                        ## amp bridge occurs only if this adj target is upstrand of next edge source
                        ## on current genome
                        this.adj.targ = this.qrs.path[, qrs.iter[[k]][1]]
                        next.adj.source = this.qrs.path[, qrs.iter[[k+1]][1]]

                        ## sanity check: this adj targ and next edge source should be on the same
                        ## chromosome and strand in the current genome
                        if (!(next.adj.source[1] == this.adj.targ[1] & next.adj.source[3] == this.adj.targ[3]))
                          stop('something wrong: next.adj.source and this.adj.targ are not on same chr and strand')

                        ## amp.bridge only if source at or before targ, which is left on pos strand and right on neg
                        is.amp.bridge = ((this.adj.targ[3] == 1 & this.adj.targ[2] <= next.adj.source[2]) |
                                         (this.adj.targ[3] == 2 & this.adj.targ[2] >= next.adj.source[3]))

                        ## forwards amp bridge extends fragment opposite <next> source with intervals from
                        ## <this target> until <next source>
                        if (is.amp.bridge)
                          {
                            ## find interval where to add amp.bridge
                            ## this interval is opposite last target in current genome
                            ## which is right of next.adj.source for neg strand and left of next.adj.source for positive
                            this.chr = which(current.state$intervals[, 'i'] == next.adj.source[1])
                            tmp.cid = ifelse(next.adj.source[3] == 1, this.chr[next.adj.source[2]-1], this.chr[next.adj.source[2]+1])
                            to.add.fid = current2frag[tmp.cid, next.adj.source[3]]

                            ## make amp.bridge frag (rid fid is.bridge to.break)
                            amp.bridge.frag = cbind(cid = c(1, -1)[this.adj.targ[3] * this.adj.targ[3]] * this.chr[this.adj.targ[2]:next.adj.source[2]])

                            amp.bridge.frag = cbind(
                              i = fragments[to.add.fid, 'i'],
                              j = 1:nrow(amp.bridge.frag), ## we are adding to the left
                              amp.bridge.frag, is.bridge = T, to.break = F)

                            ## check to see if there is something wrong
                            if (fragments[to.add.fid, 'j'] != 1)
                              warning('something is wrong: left amp bridge to be added to internal fragment')

                            ## shift j on the current fragments
                            tmp.ix = fragments[, 'i'] == fragments[to.add.fid, 'i']
                            fragments[tmp.ix, 'j'] = fragments[tmp.ix, 'j'] + nrow(amp.bridge.frag)

                            ## add amp.bridge to <left> of to.add.fid in fragments
                            bef.ix = c()
                            if (to.add.fid != 1 )
                              bef.ix = 1:(to.add.fid-1)
                            aft.ix = to.add.fid:nrow(fragments)

                            fragments = cbind(fragments[bef.ix, ], amp.bridge, fragments[aft.ix, ])

                            pix = c(bef.ix, rep(NA, nrow(amp.bridge)), aft.ix)
                            current2frag[current2frag.nna] = match(current2frag[current2frag.nna], pix)

                            ## remove reference partners from left side of this frag
                            fragment.partners[flocs[to.add.fid, 'i'], 'left'] = NA
                          }
                      }

                    ### so far, we have applied no breaks in <this> QRS iteration, we have only selected fragments to break
                    ### (and appended amp bridges to fragments broken in previous iterations)

                    ## to apply a break, we split and re-.munlist the old fragments list using to.break field and update fragment.partners

                    ## obtaining new fragments only requires relabeling the 'i' fields
                    old.fragments = fragments;
                    fragments = .munlist(lapply(split(1:ncol(fragments), fragments[, 'i'] + fragments[, 'to.break']/10),
                      function(x) fragments[x, ]), dim = 1)

                    ## updating fragment.partners requires letting separated partners
                    ## remember who they were just joined to, so they can be later rejoined
                    ## if they don't find a new partner
                    ## (unless they are amp bridges or their partner was re-fused)
                    new2old.frag = rep(NA, nrow(fragment.partners))
                    new2old.frag[fragments[,'i']] = old.fragments[, 'i']
                    fragment.partners = fragment.partners[new2old.frag, ]
                    fragment.partners[fragments[fragments[, 'to.break'], 'i'], 'right'] = fragments[fragments[, 'to.break'], 'i']+1
                    tmp.ix = which(fragments[, 'to.break'])+1
                    fragment.partners[fragments[tmp.ix, 'i'], 'left'] = fragments[tmp.ix, 'i'] - 1

                    ## to apply fusion, we only have to move rows around in the fragments matrix
                    ## moving the left fragment in the fusion immediately before the right fragment in the fusion

                    ## update this junction.fids to current fragments matrix
                    this.junction.fids = current2frag[this.junction[, c('cid','str')]]

                    ## junction connects fragment ending in interval junction.fids[1] to the fragment
                    ## beginning with junction.fids[length(junction.fids)] (with any other intervals connected in between)
                    ## (NOTE: current implementation treats junction only as interval pair)
                    left.ix = range(which(fragments[, 'i'] == fragments[junction.fids[1], 'i']))
                    right.ix = range(which(fragments[, 'i'] == fragments[junction.fids[length(junction.fids)], 'i']))

                    old.fragments = fragments;
                    fragments[c(left.ix, right.ix), ] =  fragments[left.ix[1], 'i'] ## rename fused fragment according to the left partner
                    pix = c(c(left.ix, right.ix), (1:nrow(fragments))[-c(left.ix, right.ix)]) ## save permutation

                    ## apply permutation to fragments and current2frag
                    fragments = fragments[pix, ]
                    current2frag[current2frag.nna] = match(current2frag[current2frag.nna], pix)

                    ## rename fragments so they are numbered in order (to make compatible with fragment partners)
                    fragments[, 'i'] = 1 + cumsum(c(0, diff(fragments[, 'i'])))

                    ## let the new fragment inheriting the name of the left partner inherit the partner of the right partner
                    fragment.partners[old.fragments[left.ix[1], 'i'], 'right'] = fragment.partners[old.fragments[right.ix[1], 'i'], 'right']

                    ## rewire fragment partners so that new fragment is first
                    fragment.partners = rbind(fragment.partners[old.fragments[left.ix[1], 'i'], ],
                      fragment.partners[-(old.fragments[c(left.ix[1], right.ix[1]), 'i']), ])
                  }

                ## now we have finished processing QRS

                ############
                ## Rejoin
                ############

                ## match up and rejoin "left" and "right" partners
                ## to.rejoin is adj.matrix of directed graph that should only have paths
                to.rejoin = matrix(0, nrow = nrow(fragment.partners), ncol = nrow(fragment.partners))
                to.rejoin[cbind(fragment.partners[, 'left'],
                                fragment.partners[, 'right'][match(fragment.partners[, 'left'], fragment.partners[, 'right'])])] = 1

                sources = which(colSums(to.rejoin,1)==0 & rowSums(to.rejoin,1)>0)
                rejoined.paths = list();
                for (s in 1:length(sources))
                  {
                    rejoined.paths[[s]] = sources[s]
                    next.v = which(to.rejoin[rejoined.paths[[s]][length(rejoined.paths[[s]])], ]!=0)
                    while (length(next.v)>0)
                      {
                        rejoined.paths[[s]] = c(rejoined.paths[[s]], next.v)
                        next.v = which(to.rejoin[rejoined.paths[[s]][length(rejoined.paths[[s]])], ]!=0)
                      }
                  }

                ## all other fragments will be in single fragment paths
                tmp.ix = setdiff(1:nrow(fragments.partners), unlist(rejoined.paths))
                non.rejoined.paths = split(tmp.ix, 1:length(tmp.ix))

                ## now permute and relabel fragment matrix according to paths
                all.paths = c(rejoined.paths, non.rejoined.paths)

                new.frag.ix = mapply(function(x, y) unlist(x[y]), split(1:nrow(fragments), fragments[, 'i']), all.paths, SIMPLIFY = F)
                fragments = .munlist(lapply(new.frag.ix, function(x) fragments[x, ]), dim = 1)

                ## update current2frag
                current2frag[current2frag.nna] = match(current2frag[current2frag.nna], unlist(new.frag.ix))

                ############
                ## Clean up
                ############

                ## now remove fragments that (1) do not have a "legal" telomeric interval at both ends and (2) (if is.req is not NULL)
                ## do not contain and do not contain an is.req interval
                frag.rid = split(fragments[, 'rid'], fragments[, 'i'])

                keep = sapply(frag.rid, function(x) is.tel[x[1]] & is.tel(x[length(x)]))
                if (!is.null(is.req))
                  keep = keep & sapply(frag.rid, function(x) any(is.req[x]))

                new.ix = fragments[,'i'] %in% which(keep)

                fragments = fragments[new.ix, ]
                current2frag[current2frag.nna] = match(current2frag[current2frag.nna], new.ix)

                if (verbose)
                  cat(sprintf('current fragments: \n%s\n',
                              paste(paste('[', sapply(split(fragments[, 'cid'], fragments[, 'i']), paste, collapse = ' '), ']', sep = ''), collapse = '\n')))

                ############
                ## update current.state
                ############

                ## we essentially are doing "DNA replication", since we are pulling both strands of fragments
                ## that we have joined here

                ## double stranded intervals are format i j + -
                remove = current.state$intervals[, 'i'] %in% check.out[, 'chr'] ## mark current state chromosomes to remove

                ## construct intervals to check in
                intervals.check.in = cbind(
                  i = fragments[, 'i']+0.1, ## provide fake chromosome label
                  j = fragments[, 'j'],
                  '+'  = sapply(fragments[, 'cid'], function(x) if (x>0) current.state$intervals[abs(x), '+'] else current.state$intervals[abs(x), '-']),
                  '-'  = sapply(fragments[, 'cid'], function(x) if (x>0) current.state$intervals[abs(x), '-'] else current.state$intervals[abs(x), '+'])
                  )

                ## renumber chromosomes
                tmp.intervals = rbind(intervals.check.in, current.state$intervals[!remove,])
                new.intervals = .munlist(split(1:nrow(tmp.intervals), tmp.intervals[, 'i']), function(x) tmp.intervals[x, c('+', '-')])

                if (verbose)
                  cat(sprintf('----------------------------\nChecking in chromosomes: \n%s\n-----------------------------\n',
                              paste(sapply(split(1:nrow(intervals.check.in), fragments[, 'i']),
                                           function(x) paste(new.intervals[x, '+'], new.intervals[x, '-'], sep = '\t', collapse = '\n'),
                                           collapse = '\n==========================\n'))))

                ## construct new state, including mapping to previous id's
                current.state = list(
                  intervals = new.intervals,
                  rid2cid = .rid2cid(new.intervals)
                  )

                ## if this is the second or later qrs in a multi-qrs event, then maintain "cid2prev" to map to the
                ## last current state in the history
                if (qrs.i>1)
                  current.state$cid2prev = current.state$cid2prev[c(fragments[, 'cid'], which(!remove))]
                else
                  current.state$cid2prev = c(fragments[, 'cid'], which(!remove))

                ## map last element in qrs.paths to new current.state cid
                ## this is where the next event in a multi qrs event will attempt to be instantiated
                last.qrs = current2frag[qrs.paths[length(qrs.paths)], ifelse(qrs.paths[length(qrs.paths)]>0, 1, 2)]
              }

            ## if we have aborted, then roll history back to last current.state
            ## since done.this.round[k] = T, we won't redo this event in this round
            if (abort)
              {
                current.state = history[[length(history)]]
                history = history[[-length(history)]]
                step = step - 1
              }
            else ## otherwise, update done.events with all events that intersect with an junction covered by this event
              {
                done.events = done.events | rowSums(event2junction[, which(event2junction[k, ])]) != 0
                done.this.round = rep(F, length(done.events)) ## reset set of events that have been done this round
              }
          }
        else ## do a chromosomal event
          {
            if (force.cn) ## then need to instantiate GRanges event
              {
                this.event = events[[k]]

                if (verbose)
                  cat(sprintf('Implementing event %s (CN event) on chroms %s\n', k, unique(as.character(seqnames(this.event)))))

                if (is.null(this.event$type))
                  stop('All GRanges specifying copy number events must have $type field set to "gain" or "loss"')

                ev2tile = gr.findoverlaps(this.event, kag$tile, pintersect = T)

                ## match each row of this event to possible instantiations
                ev2cid = lapply(split(current.state$rid2cid[ev2tile$subject.id], ev2tile$query.id), function(x) do.call('c', x))
                ev2cid.type = this.event$type[unique(ev2tile$query.id)] ## gain or loss

                ## now pick a random instantiation per event
                ev.cid = lapply(ev2cid, function(x) if (length(x)>0) x[sample(length(x))] else c())

                gain.cid = loss.cid = c()
                ev.gain = which(ev2cid.type == 'gain')
                ev.loss= which(ev2cid.type == 'loss')

                ## choose randomly among instantiations to choose chromosomes to lose or gain
                ## (choosing with replacement,so if a cid exists in multiple rows in this.event then they
                ## may be both chosen)
                if (length(ev.gain)>0)
                  gain.cid = unlist(lapply(ev.cid[ev.gain], function(x) if (length(x)>0) abs(x[sample(length(x),1)]) else c()))

                if (length(ev.loss)>0)
                  loss.cid = unlist(lapply(ev.cid[ev.loss], function(x) if (length(x)>0) abs(x[sample(length(x),1)]) else c()))

                gain  = current.state$intervals[gain.cid, 'i']
                loss  = current.state$intervals[loss.cid, 'i']

                current.state = .chrom_change(current.state, gain = gain, loss = loss)

                done.events[k] = T
              }
            else ## otherwise, random event
              {
                if (runif(1)<p.wgd & !force.cn)
                  current.state = .chrom_change(current.state, gain = 1:length(current.state$intervals))
                else
                  {
                    if (runif(1)<p.chromgain)
                      {
                        num.chrom = pmin(length(current.state$intervals), rpois(1, lambda = lambda.chromgain))
                        prob = sapply(current.state$intervals, function(x, w) sum(w[x[1, ]]), w = as.numeric(width(kag$tile)))
                        gain = sample(1:length(current.state$intervals), replace = F, prob = prob, size = num.chrom)
                        current.state = .chrom_change(current.state, gain = gain)
                      }
                    else
                      {
                        num.chrom = pmin(length(current.state$intervals), rpois(1, lambda = lambda.chromloss))
                        prob = sapply(current.state$intervals, function(x, w) sum(w[x[1, ]]), w = as.numeric(width(kag$tile)))
                        loss = sample(1:length(current.state$intervals), replace = F, prob = prob, size = num.chrom)
                        current.state = .chrom_change(current.state, loss = loss)
                      }
                  }
              }
          }
      }
    if (full)
      return(history)
    else
      return(current.state)
  }

#' @name all.paths
#' @title all.paths
#' @description
#'
#'
#' Low level function to enumerate all elementary paths and cycles through graph
#'
#' takes directed graph represented by n x n binary adjacency matrix  A and outputs all cycles and paths between source.vertices, sink.vertices
#'
#'
#' @param A nxn adjacency matrix
#' @param all logical flag, if all = T, will include all sources (parentless vertices) and sinks (childless vertices) in path computati
#' @param ALL logical flag, if ALL = T, will also include vertices without outgoing and incoming edges in paths
#' @param sources graph indices to treat as sources (by default is empty)
#' @param sinks graph indices to treat as sinks (by default is empty)
#' @param verbose logical flag
#' @return list of integer vectors corresponding to indices in A (i.e. vertices)
#' $paths = paths indices
#' $cycles = cycle indices
#' @export
all.paths = function(A, all = F, ALL = F, sources = c(), sinks = c(), source.vertices = sources, sink.vertices = sinks,
  exclude = NULL, ## specifies illegal subpaths, all such paths / cycles and
                  ## their supersets will be excluded, specified as k x nrow(A) matrix of vertex sets
  verbose = FALSE,...)
  {
    require(igraph)

    blank.vertices = which(Matrix::rowSums(A)==0 & Matrix::colSums(A)==0)

    if (ALL)
      all = T

    if (all)
      {
        source.vertices = which(Matrix::rowSums(A)>0 & Matrix::colSums(A)==0)
        sink.vertices = which(Matrix::colSums(A)>0 & Matrix::rowSums(A)==0)
      }

    out = list(cycles = NULL, paths = NULL)

    node.ix = which(Matrix::rowSums(A!=0)>0 | Matrix::colSums(A!=0)>0)
    if (length(node.ix)==0)
      return(out)

    A = A[node.ix, node.ix]

    if (!is.null(exclude))
      exclude = sign(abs(exclude[, node.ix]))

    ij = which(A!=0, arr.ind = T)
    B = sparseMatrix(c(ij[,1], ij[,2]), rep(1:nrow(ij), 2), x = rep(c(-1, 1), each = nrow(ij)), dims = c(nrow(A), nrow(ij)))
    I = diag(rep(1, nrow(A)))

    source.vertices = setdiff(match(source.vertices, node.ix), NA)
    sink.vertices = setdiff(match(sink.vertices, node.ix), NA)

    B2 = cBind(B, I[, source.vertices, drop = FALSE], -I[, sink.vertices, drop = FALSE])

    if (verbose)
      cat(sprintf('Computing paths for %s vertices and %s edges\n', nrow(B2), ncol(B2)))

    K = convex.basis(B2, verbose = verbose, exclude.range = exclude, ...)

    if (all(is.na(K)))
      return(out)

    K = K[, Matrix::colSums(K[1:ncol(B), ,drop = FALSE])!=0, drop = FALSE] ## remove any pure source to sink paths

    is.cyc = Matrix::colSums(B %*% K[1:ncol(B), ,drop = FALSE]!=0)==0

    out$cycles = lapply(which(is.cyc),
      function(i)
      {
        k = which(K[1:ncol(B), i]!=0)
        v.all = unique(as.vector(ij[k, , drop = FALSE]))
        sG = graph.edgelist(ij[k, , drop = FALSE])
        tmp.v = v.all[c(1,length(v.all))]
        p.fwd = get.shortest.paths(sG, tmp.v[1], tmp.v[2])
        p.bwd = get.shortest.paths(sG, tmp.v[2], tmp.v[1])
        return(node.ix[unique(unlist(c(p.fwd, p.bwd)))])
      })

    out$paths = lapply(which(!is.cyc),
      function(i)
      {
        k = K[1:ncol(B), i]
        eix = which(k!=0)
        v.all = unique(as.vector(ij[eix, , drop = FALSE]))
        sG = graph.edgelist(ij[eix, , drop = FALSE])
        io = B %*% k
        v.in = which(io<0)[1]
        v.out = which(io>0)[1]
        return(node.ix[unlist(get.shortest.paths(sG, v.in, v.out))])
      })

    if (length(out$cycles)>0)
      {
        tmp.cix = cbind(unlist(lapply(1:length(out$cycles), function(x) rep(x, length(out$cycles[[x]])))), unlist(out$cycles))
        out$cycles = out$cycles[!duplicated(as.matrix(sparseMatrix(tmp.cix[,1], tmp.cix[,2], x = 1)))]
      }

    if (length(out$paths)>0)
      {
        tmp.pix = cbind(unlist(lapply(1:length(out$paths), function(x) rep(x, length(out$paths[[x]])))), unlist(out$paths))
        out$paths = out$paths[!duplicated(as.matrix(sparseMatrix(tmp.pix[,1], tmp.pix[,2], x = 1)))]
      }

    if (ALL & length(blank.vertices)>0)
      out$paths = c(out$paths, lapply(blank.vertices, identity))

    return(out)
  }



#' @name collapse.paths
#' @title collapse.paths
#' @description
#'
#' collapse simple paths in a graph G (adjacency matrix or igraph object)
#' returns m x m new adjacency matrix and map of old vertex id's to new ones
#' $adj = m x m matrix
#' #map = length n with indices 1 .. m
#'
###############################################
collapse.paths = function(G, verbose = T)
{
  if (inherits(G, 'igraph'))
      G = G[,]

  out = G!=0

  if (verbose)
      cat('graph size:', nrow(out), 'nodes\n')

  ## first identify all nodes with exactly one parent and child to do initial collapsing of graph
  singletons = which(Matrix::rowSums(out)==1 & Matrix::colSums(out)==1)

  if (verbose)
      cat('Collapsing simple paths..\n')

  sets = split(1:nrow(G), 1:nrow(G))
  if (length(singletons)>0)
      {
          tmp = out[singletons, singletons]
          cl = clusters(graph(as.numeric(t(which(tmp, arr.ind = TRUE))), n = nrow(tmp)), 'weak')$membership
          dix = unique(cl)
          if (length(dix)>0)
              {
                  for (j in dix)
                      {
                          if (verbose)
                              cat('.')

                          ## grab nodes in this cluster
                          setj = singletons[which(cl == j)]

                          ## move all members into a single set
                          sets[setj[1]] = list(setj)
                          sets[setj[-1]] = list(NULL)

                          ## connect this node to the parent and child of the set
                          parent = setdiff(which(Matrix::rowSums(out[, setj, drop = FALSE])>0), setj)
                          child = setdiff(which(Matrix::colSums(out[setj, , drop = FALSE])>0), setj)
                          out[setj, c(setj, child)] = FALSE
                          out[c(setj, parent), setj] = FALSE
                          out[parent, setj[1]] = TRUE
                          out[setj[1], child] = TRUE
                      }
              }
      }

  if (verbose)
      cat('done\nnow fixing branches\n')

  todo = rep(FALSE, nrow(G))
  todo[Matrix::rowSums(out)==1 | Matrix::colSums(out)==1] = TRUE

  while (sum(todo)>0)
      {
        sets.last = sets
        out.last = out

        if (verbose)
            if ((sum(todo) %% 200)==0)
                cat('todo:', sum(todo), 'num sets:', sum(!sapply(sets, is.null)), '\n')

        i = which(todo)[1]

        todo[i] = F

        child = which(out[i, ])
        parent = which(out[,i])

        if (length(child)==1 & length(parent)==1) ## if there is exactly one child and one parent then we want to merge with one or both
            {
                ## if i-child has no other parents and i-parent has no other child
                ## then merge i, i-parent and i-child
                if (sum(out[,  child])==1 & sum(out[parent, ])==1)
                    {
                        grandch = which(out[child, ])
                        if (length(grandch)>0)
                            {
                                out[parent, grandch] = TRUE  ## parent inherits grandchildren of i
                                out[child, grandch] = FALSE
                            }
                        out[parent, i] = FALSE ## remove node i's edges
                        out[i, child] = FALSE
                        sets[[parent]] = c(sets[[parent]], sets[[child]], sets[[i]])
                        sets[c(i, child)] = list(NULL)
                        todo[child] = F ## no longer have to do i-child, since they have already been merged with parent
                    }
                ## otherwise if either i-child has no other parent or i-parent has no other children (but not both)
                ## then connect i-parent to i-child, but do not merge them (but merge ONE of them with i)
                else if (sum(out[,  child])==1 | sum(out[parent, ])==1)
                    {
                        ## if parent has no other children then merge with him
                        if (sum(out[parent, ])==1)
                            sets[[parent]] = c(sets[[parent]], sets[[i]])
                        else
                            sets[[child]] = c(sets[[child]], sets[[i]])

                        out[parent, child] = TRUE
                        out[parent, i] = FALSE ## remove node i's edges
                        out[i, child] = FALSE
                        sets[i] = list(NULL)
                    }
            }
        else if (length(child)==1 & length(parent)>1) ## if i has more than one parent but one child, we merge with child if child has no other parents
            {
                if (sum(out[, child])==1)
                    {
                        sets[[child]] = c(sets[[child]], sets[[i]])
                        out[parent, child] = TRUE
                        out[parent, i] = FALSE ## remove node i's edges
                        out[i, child] = FALSE ## remove node i's edges
                        sets[i] = list(NULL)
                    }


            }
        else if (length(child)>1 & length(parent)==1) ## if i has more than one child but one parent, then merge with parent if parent has no other children
            {
                if (sum(out[parent, ])==1)
                    {
                        sets[[parent]] = c(sets[[parent]], sets[[i]])
                        out[parent, child] = TRUE
                        out[parent, i] = FALSE ## remove node i's edges
                        out[i, child] = FALSE ## remove node i's edges
                        sets[i] = list(NULL)
                    }
            }
    }

  slen = sapply(sets, length)
  ix = which(slen>0)
  map = rep(NA, nrow(G))
  map[unlist(sets)] = match(rep(1:length(sets), slen), ix)
  out = out[ix, ix]
  colnames(out) = rownames(out) = NULL

  return(list(adj = out, map = map, sets = split(1:length(map), map)))
}


##############################################
#' @name sparse_subset
#' @title sparse_subset
#' @description
#'
#' given k1 x n matrix A and k2 x n matrix B
#' returns k1 x k2 matrix C whose entries ij = 1 if the set of nonzero components of row i of A is
#' a (+/- strict) subset of the nonzero components of row j of B
#'
sparse_subset = function(A, B, strict = FALSE, chunksize = 100, quiet = FALSE)
  {
    nz = Matrix::colSums(as.matrix(A)!=0, 1)>0

    if (is.null(dim(A)) | is.null(dim(B)))
      return(NULL)

    C = sparseMatrix(i = c(), j = c(), dims = c(nrow(A), nrow(B)))

    for (i in seq(1, nrow(A), chunksize))
      {
        ixA = i:min(nrow(A), i+chunksize-1)
        for (j in seq(1, nrow(B), chunksize))
          {
            ixB = j:min(nrow(B), j+chunksize-1)

            if (length(ixA)>0 & length(ixB)>0 & !quiet)
              cat(sprintf('\t interval A %s to %s (%d) \t interval B %d to %d (%d)\n', ixA[1], ixA[length(ixA)], nrow(A), ixB[1], ixB[length(ixB)], nrow(B)))
            if (strict)
              C[ixA, ixB] = (sign((A[ixA, , drop = FALSE]!=0)) %*% sign(t(B[ixB, , drop = FALSE]!=0))) * (sign((A[ixA, , drop = FALSE]==0)) %*% sign(t(B[ixB, , drop = FALSE]!=0))>0)
            else
              C[ixA, ixB] = (sign(A[ixA, nz, drop = FALSE]!=0) %*% sign(t(B[ixB, nz, drop = FALSE]==0)))==0
          }
      }

    return(C)
  }




##############################################
#' @name reads.to.walks
#' @title reads.to.walks
#' @description
#'
#' Utility function to realign reads to walks.
#'
#' Takes bam and collection of walks (GRanges list of signed intervals on hg19 or other BSgenome)
#' pulls reads in regions of walks, then uses bwa mem to realign reads to walks, returns paths to new bams
#' These bams are in "walk coordinates"
#'
#' Assumes bwa mem installed an runnable from command line.
#'
#' @param bam path to bam file
#' @param walks GRangesList of walks
#' @param outdir outdir to compute into
#' @param hg human genome sequence BSGenome or ffTrack
#' @param mc.cores number of cores
#' @param insert.size >= max insert size of library so that all relevant read pairs are extracted from the original bam
#' @return paths to new bams
#' These bams are in "walk coordinates"
#' @export
##############################################
reads.to.walks = function(bam, walks, outdir = './test', hg = skidb::read_hg(fft = T), mc.cores = 1, insert.size = 1e3, verbose = T)
  {
    system(paste('mkdir -p', outdir))

    .pairs2fq = function(x)
      {
        x.gr = grl.unlist(x)
        x.gr = x.gr[!is.na(x.gr$seq)]
        x.gr$first = bamflag(x.gr$flag)[,'isFirstMateRead']!=0
        x.gr$unmapped = bamflag(x.gr$flag)[,'isUnmappedQuery']!=0
        x.gr1 = x.gr[x.gr$first]
        x.gr2 = x.gr[!x.gr$first]
        x.gr1 = x.gr1[x.gr1$qname %in% x.gr2$qname]
        x.gr2 = x.gr2[match(x.gr1$qname, x.gr2$qname), ]
        if (any(ix <- as.logical(strand(x.gr1)=='-') & !x.gr1$unmapped)) ## rev comp sequence and reverse quality scores
          {
            x.gr1$seq[ix] = as.character(reverseComplement(DNAStringSet(x.gr1$seq[ix])))
            x.gr1$qual[ix] = sapply(x.gr1$qual[ix], function(x) rawToChar(rev(charToRaw(x))))
          }

        if (any(ix <- as.logical(strand(x.gr2)=='-') & !x.gr1$unmapped)) ## rev comp sequence and reverse quality scores
          {
            x.gr2$seq[ix] = as.character(reverseComplement(DNAStringSet(x.gr2$seq[ix])))
            x.gr2$qual[ix] = sapply(x.gr2$qual[ix], function(x) rawToChar(rev(charToRaw(x))))
          }
        fq1 = as.vector(t(cbind(paste('@', x.gr1$qname, sep = ''), x.gr1$seq, '+', x.gr1$qual)))
        fq2 = as.vector(t(cbind(paste('@', x.gr2$qname, sep = ''), x.gr2$seq, '+', x.gr2$qual)))
        return(list(fq1 = fq1, fq2 = fq2))
      }

    if (is.null(names(walks)))
      names(walks) = paste('walk', 1:length(walks))

    outdir = normalizePath(outdir)

    reads.fq = paste(outdir, '/reads.', 1:2, '.fq', sep = '')
    walk.fa = paste(outdir, '/', names(walks), '.fa', sep = '')
    walks.gff = paste(outdir, '/walks.gff', sep = '')

    if (!all(file.exists(walk.fa)))
      {
        if (verbose)
          cat('Fetching walk sequences\n')

        walk.seq = ffTrack::get_seq(hg, walks, mc.cores = mc.cores)
        names(walk.seq) = names(walks)

        sapply(1:length(walks), function(x) writeXStringSet(walk.seq[x], walk.fa[x]))

        if (is.list(walk.seq))
          {
            c(walk.seq[[1]], walk.seq[[2]]) ## weird R ghost
            walk.seq = do.call('c', walk.seq)
          }

        ## write compiled fasta
        writeXStringSet(walk.seq, paste(outdir, '/walks.fa', sep = ''))
      }

    if (!all(file.exists(walks.gff)))
      {
        tmp = gChain::spChain(walks)$y;
        export.gff(split(tmp, seqnames(tmp)), walks.gff)
      }

    ## grab reads from region and output to fq
    if (!all(file.exists(reads.fq)))
      {
        if (verbose)
          cat(sprintf('Fetching mapped and unmapped reads associated with region from %s\n', bam))

        reads = read.bam(bam, intervals = streduce(unlist(walks)+insert.size))
        reads.seq = .pairs2fq(reads)

        ## write read fastqs
        if (verbose)
          cat(sprintf('Writing fastqs for %s read pairs\n', length(reads[[1]])))

        bla = mapply(function(x, y) writeLines(x, y), reads.seq, reads.fq)
      }

    if (verbose)
      cat('Indexing walk fastas')

    walks.faidx = paste(walk.fa, '.bwt', sep = '')

    if (any(ix <- !file.exists(walks.faidx)))
      mclapply(walk.fa[ix], function(x) system(paste('bwa index', x)), mc.cores = mc.cores)

    if (verbose)
      cat('Running bwa mem\n')

    mclapply(1:length(walk.bam), function(x)
             {
               cmd = sprintf('bwa mem %s %s %s | samtools view -bSh -F4 - > %s; samtools sort %s %s; samtools index %s',
                     walk.fa[x], reads.fq[1], reads.fq[2], walk.bam[x], walk.bam[x], gsub('.bam$', '', walk.bam[x]), walk.bam[x])
               if (verbose)
                 cat(cmd, '\n')
               system(cmd)
               }, mc.cores = mc.cores)

    if (verbose)
      cat('Done\n')

    return(walk.bam)
  }



#' @name gr.breaks
#' @title gr.breaks
#' @description
#'
#' Break GRanges at given breakpoints into disjoint gr
#'
#' @author Xiaotong Yao
#' @import GenomicRanges
#' @param bps \code{GRanges} of width 1, locations of the bp; if any element width
#' larger than 1, both boundary will be considered individual breakpoints
#' @param query a disjoint \code{GRanges} object to be broken
#' @return \code{GRanges} disjoint object at least the same length as query,
#' with a metadata column \code{qid} indicating input index where new segment is from
#' @export
gr.breaks = function(bps=NULL, query=NULL){
   ## ALERT: big change! input parameter shuffled!

   ## if bps not provided, return back-traced disjoin wrapper
   if (is.null(bps)) {
       return(query)
   }
   else {
       ## only when bps is given do we care about what query is
       if (is.null(query)){
           message("Trying chromosomes 1-22 and X, Y.")
           query = hg_seqlengths()
           if (is.null(query)){
               message("Default BSgenome not found, let's hardcode it.")
               cs = system.file("extdata",
                                "hg19.regularChr.chrom.sizes", package = "gUtils")
               sl = read.delim(cs, header=FALSE, sep="\t")
               sl = setNames(sl$V2, sl$V1)
               query = gr.stripstrand(si2gr(sl))
           }
       }

       ## preprocess query
       if (!isDisjoint(query)){
           warning("Warning: Query GRanges not disjoint.")
           queryDj = disjoin(query)
           queryDj$qid = queryDj %N% query ## only retain the first occurence
           values(queryDj) = cbind(values(queryDj),
                                   as.data.table(values(query))[queryDj$qid])
           query = queryDj
       }
       else {
           if ("qid" %in% colnames(values(query))){
               warning("Warning: 'qid' col in query overwritten.")
           }
           query$qid = seq_along(query)
       }

       ## preprocess bps
       ## having meta fields? remove them!
       bps = bps[, c()]

       ## remove things outside of ref
       oo.seqlength = which(start(bps)<1 | end(bps)>seqlengths(bps)[as.character(seqnames(bps))])
       if (length(oo.seqlength)>0){
           warning("Warning: Some breakpoints out of chr lengths. Removing.")
           bps = bps[-oo.seqlength]
       }

       if (any(!is.null(names(bps)))){
           warning("Removing row names from bps.")
           names(bps) = NULL
       }

       ## having strand info? remove it!
       if (any(strand(bps)!="*")){
           warning("Some breakpoints have strand info. Force to '*'.")
           bps = gr.stripstrand(bps)
       }

       ## solve three edge cases
       if (any(w.0 <- (width(bps)<1))){
           warning("Some breakpoint width==0.")
           ## right bound smaller coor
           ## and there's no negative width GR allowed
           bps[which(w.0)] = gr.start(bps[which(w.0)]) %-% 1
       }
       if (any(w.2 <- (width(bps)==2))){
           warning("Some breakpoint width==2.")
           ## this is seen as breakpoint by spanning two bases
           bps[which(w.2)] = gr.start(bps[which(w.2)])
       }
       if (any(w.l <- (width(bps)>2))){
           ## some not a point? turn it into a point
           warning("Some breakpoint width>1.")
           rbps = gr.end(bps[which(w.l)])
           lbps = gr.start(bps[which(w.l)])
           start(lbps) = pmax(start(lbps)-1, 1)
           bps = c(bps[which(!w.l)], streduce(c(lbps, rbps)))
       }

       bps$inQuery = bps %^% query
       if (any(bps$inQuery==F)){
           warning("Some breakpoint not within query ranges.")
       }

       ## label and only consider breakpoints not already at the boundary of query
       bps$inner = bps$inQuery
       bps$inner[which(bps %^% gr.start(query) | bps %^% gr.end(query))]=F
       ## maybe no inner bp at all, then no need to proceed
       if (!any(bps$inner)){
           return(query)
       }
       bpsInner = bps %Q% (inner==T)
       ## map query and inner breakpoints
       qbMap = gr.findoverlaps(query, bpsInner)
       mappedQ = seq_along(query) %in% qbMap$query.id
       ## raw coors to construct ranges from
       tmpRange = data.table(qid2 = qbMap$query.id,
                             startFrom = start(query[qbMap$query.id]),
                             breakAt = start(bpsInner[qbMap$subject.id]),
                             upTo = end(query[qbMap$query.id]))
       tmpCoor = tmpRange[, .(pos=sort(unique(c(startFrom, breakAt, upTo)))), by=qid2]

       ## construct new ranges
       newRange = tmpCoor[, .(start=pos[-which.max(pos)],
                              end=pos[-which.min(pos)]), by=qid2]
       newRange[, ":="(chr = as.vector(seqnames(query)[qid2]),
                       strand = as.vector(strand(query)[qid2]))]
       newRange$start = newRange[, ifelse(start==min(start), start, start+1)]

       ## put together the mapped and broken
       newGr = GRanges(newRange, seqinfo = seqinfo(query))
       values(newGr) = values(query)[newGr$qid2, , drop=F] ## preserve the input metacol
       ## with the intact not mapped part of query
       output = sort(c(newGr, query[!mappedQ]))
       ## %Q% (order(strand, seqnames, start))
       ## browser()
       return(output)
   }
}




#' @name ra.dedup
#' @export
ra.dedup = function(grl, pad=500, ignore.strand=FALSE){

   if (!is(grl, "GRangesList")){
       stop("Error: Input must be GRangesList!")
   }

   if (any(elementNROWS(grl)!=2)){
       stop("Error: Each element must be length 2!")
   }

   if (length(grl)==0 | length(grl)==1){
       return(grl)
   }

   if (length(grl) > 1){
       ix.pair = as.data.table(
          ra.overlaps(grl, grl, pad=pad, ignore.strand = ignore.strand))[ra1.ix!=ra2.ix]
       if (nrow(ix.pair)==0){
           return(grl)
       }
       else {
           dup.ix = unique(rowMax(as.matrix(ix.pair)))
           return(grl[-dup.ix])
       }
   }
}

#' @name ra.duplicated
#' @export
ra.duplicated = function(grl, pad=500, ignore.strand=FALSE){

   if (!is(grl, "GRangesList")){
       stop("Error: Input must be GRangesList!")
   }

   if (any(elementNROWS(grl)!=2)){
       stop("Error: Each element must be length 2!")
   }

   if (length(grl)==0){
       return(logical(0))
   }

   if (length(grl)==1){
       return(FALSE)
   }

   if (length(grl)>1){

       ix.pair = as.data.table(ra.overlaps(grl, grl, pad=pad, ignore.strand = ignore.strand))[ra1.ix!=ra2.ix]

       if (nrow(ix.pair)==0){
           return(rep(FALSE, length(grl)))
       }
       else {
           dup.ix = unique(rowMax(as.matrix(ix.pair)))
           return(seq_along(grl) %in% dup.ix)
       }
   }
}

#############################
#' @name rrbind
#' @title rrbind
#'
#' @description
#'
#' like rbind, but takes the intersecting columns of the dfs
#'
#' if union flag is used then will take union of columns (and put NA's for columns of df1 not in df2 and vice versa)
#'
#' @param union logical flag whether to take union of columns (if FALSE takes intersection)
#' @author Marcin Imielinski
############################
rrbind = function(..., union = T)
{
    dfs = list(...);  # gets list of data frames
    if (any(ix <- sapply(dfs, function(x) class(x)[1])!='data.frame'))
    dfs[ix] = lapply(dfs[ix], as.data.frame)

    dfs = dfs[!sapply(dfs, is.null)]
    dfs = dfs[sapply(dfs, ncol)>0]

    ## defactorize (need to do to cat without introducing NA's in weird places)
    dfs = lapply(dfs, function(x) { for (y in names(x)) if (is.factor(x[,y])) x[, y] = as.character(x[, y]); return(x)})

    names.list = lapply(dfs, names);
    classes = unlist(lapply(dfs, function(x) sapply(names(x), function(y) class(x[, y]))))
    cols = unique(unlist(names.list));
    unshared = lapply(names.list, function(x) setdiff(cols, x));
    unshared.u = unique(unlist(unshared))
    ix = which(sapply(dfs, nrow)>0)
    expanded.dfs = lapply(ix, function(x)
    {
        dfs[[x]][, unshared[[x]]] = as.character(NA);
        return(dfs[[x]][, cols, drop = F])
    })

    out = do.call('rbind', expanded.dfs);

    if (any(uix <<- which(classes[unshared.u] != 'character')))
    {
        ix = match(unshared.u, names(out))
        for (j in uix) ### HACK to prevent stupid class mismatches leading to NA BS
        out[, ix[j]] = as(out[, ix[j]], classes[unshared.u[j]])
    }

    if (!union)
    {
        shared = setdiff(cols, unique(unlist(unshared)))
        out = out[, shared];
    }

    return(out)
}



#' @name mski_alpha
#' @title mski_alpha
#' @description
#'
#' Wrapper combining 'col2rgb()' and 'rgb()' for single colors
#'
#' Originally called 'alpha()' in 'mskilab/bamUtils'
#'
#' @param col string Any of the three kinds of R color specifications, i.e., either a color name (as listed by colors()), a hexadecimal string of the form "#rrggbb" or "#rrggbbaa" (see rgb), or a positive integer i meaning palette()[i]
#' @param alpha boolean indicating whether the alpha channel (opacity) values should be returned (default = FALSE)
#' @return string
#' @export
mski_alpha = function(col, alpha = FALSE)
{
  col.rgb = col2rgb(col, alpha = alpha)
  out = rgb(red = col.rgb['red', ]/255, green = col.rgb['green', ]/255, blue = col.rgb['blue', ]/255, alpha = alpha)
  names(out) = names(col)
  return(out)
}


#' @name get.var.col
#' @title Simple function storing default variant color scheme
#' @description
#'
#' Simple function storing default variant color scheme
#'
#' Originally in 'mskilab/bamUtils'
#'
#' @return vector of default variant colors
#' @export
get.varcol = function()
{
    VAR.COL = c('XA' = 'green', 'XG' = 'brown', 'XC' = 'blue', 'XT' = 'red', 'D' = 'white',
    'I'= 'purple', 'N' = alpha('gray', 0.2), 'XX' = 'black', 'S' = alpha('pink', 0.9))
    return(VAR.COL)
}


## #' @name counts2rpkm
## #' @title Compute rpkm counts from counts
## #' @description
## #'
## #' Takes 'Rsamtools::countBam()'' (or 'bam.cov.gr()') output "counts" and computes RPKM by aggregating across "by" variable
## #'
## #' @param counts data.table or data.frame with records, width fields
## #' @param by string Field to group counts by
## #' @note The denominator (i.e. total reads) is just the sum of counts$records
## #' @return TO BE DONE
## #' @export
counts2rpkm = function(counts, by)
{
    if (missing(counts) | missing(by)){
        stop('Error: "counts2rpkm()" requires both arguments "counts" and "by". Please see documentation for details.')
    }
    out = aggregate(1:nrow(counts), by = list(by), FUN = function(x) sum(counts$records[x])/sum(counts$width[x]/1000))
    out[,2] = out[,2] / sum(counts$records) * 1e6
    names(out) = c('by', 'rpkm')
    return(out)
}





#' @name oneoffs
#' @title Calls samtools mpileup to dump tsv of "one off" variants / sites
#' @description
#'
#' Calls samtools mpileup to dump tsv of "one off" variants / sites (i.e. that are present in exactly one read per site)
#'
#' @param out.file string Path to file in which to dump tsv 
#' @param bam string Path to BAM file 
#' @param ref tring Path to reference FASTA
#' @param min.bq integer Minimum base quality
#' @param min.mq integer Minimum mapping quality
#' @param indel boolean Flag whether to collect one off indels (default is substitution)
#' @param chunksize integer Number of mpileup lines to put into memory
#' @param verbose boolean Flag to increase verbosity (default = TRUE)
#' @note The denominator (i.e. total reads) is just the sum of counts$records
#' @export
oneoffs = function(out.file, bam, ref, min.bq = 30, min.mq = 60, indel = FALSE, chunksize = 1e4, verbose = TRUE)
{   
    if (indel){
        cmd = sprintf('samtools mpileup -x -B -Q %s -q %s -s -f %s %s | grep -P "\\w+\\s\\w+\\s\\w+\\s[\\,\\.]*[\\+\\-]\\d+[ACGTNacgtn]+[\\,\\.]*\\s"', min.bq, min.mq, ref, bam)
    }
    else{
        cmd = sprintf('samtools mpileup -x -B -Q %s -q %s -s -f %s %s | grep -P "\\w+\\s\\w+\\s\\w+\\s[\\,\\.]*[ACGTacgt][\\,\\.]*\\s"', min.bq, min.mq, ref, bam)
    }
  
    p = pipe(cmd, open = 'r')

    start = Sys.time()
    fields = c('chr', 'pos', 'ref', 'cov', 'alt', 'bq', 'mq')

    i = nv = nl = 0

    while (length(chunk <- readLines(p, n = chunksize)) > 0){

        tab = fread(paste(chunk, collapse = '\n'), sep = '\t', header = FALSE)
        setnames(tab, fields)
        tab[ ,varnum := 1:.N]

        if (indel){

            tab[, left.pad := nchar(gsub("[\\+\\-].*", '', alt))]
            tab[, wid := as.numeric(gsub('.*([\\+\\-]\\d+).*', '\\1', alt))]
            tab[, var := mapply(function(x,i) substr(x, 1, i),
                gsub('.*[\\+\\-]\\d+([ACGTNacgtn]+).*', '\\1', alt),
                abs(wid))]
            tab[wid>0, bq := mapply(function(x, i) substr(x, i, i), bq, wid)]
            tab[wid<0, bq := NA]
            varb = tab[, .(chr, pos, alt, wid, mq = NA, bq = NA)]                          
        }
        else{
            varb = tab[, .(chr = chr, pos = pos, alt = unlist(strsplit(alt, '')),
                wid = 0,
                bq = utf8ToInt(unlist(strsplit(bq, '')))-33,
                mq = utf8ToInt(unlist(strsplit(mq, '')))-33), by = varnum][!(alt %in% c(".", ",")), ]
        }

    
        fwrite(varb, out.file, append = (i>0))
        nv = nv + nrow(varb)
        nl = nl + length(chunk)
        i = i+1
        if (verbose){
            message('Wrote total of ',
                nl, ' variants to ", out.file, ". Now at coordinate ',
                varb[nrow(varb), sprintf("chr%s %s", chr, prettyNum(pos, ','))])
            print(Sys.time() - start)
        }
    }
  
    close(p)
    if (verbose){
        message('Done writing ', out.file)
    }
}


anno.gwalk = function(gwalks, hop.thresh = 1e5, bighop.thresh = 1e8)
{
  message('Computing longest hop run 1e5')
  values(gwalks)$hop.thresh = hop.thresh
  values(gwalks)$bighop.thresh = bighop.thresh
  values(gwalks)$longest.hoprun1e5 = grl.eval(gwalks,
                                            max(c(0, table(label.runs(abs(dist.nostrand)>hop.thresh & frag.pos<1e5)))),
                                            condition = dist != 1)

  message('Computing longest hop run 1e4')
  values(gwalks)$longest.hoprun1e4 = grl.eval(gwalks,
                                            max(c(0, table(label.runs(abs(dist.nostrand)>hop.thresh & frag.pos<1e4)))),
                                            condition = dist != 1)

  message('Computing longest hop run 1e3')
  values(gwalks)$longest.hoprun1e3 = grl.eval(gwalks,
                                            max(c(0, table(label.runs(abs(dist.nostrand)>hop.thresh & frag.pos<1e3)))),
                                            condition = dist != 1)

  message('Computing number 4 hop runs 1e4')
  values(gwalks)$num4hoprun1e4 = grl.eval(gwalks,
                                        sum(table(label.runs(abs(dist.nostrand)>hop.thresh & frag.pos<1e4))>3),
                                        condition = dist != 1)

  message('Computing number 4 hop runs 1e5')
  values(gwalks)$num4hoprun1e5 = grl.eval(gwalks,
                                        sum(table(label.runs(abs(dist.nostrand)>hop.thresh & frag.pos<1e5))>3),
                                        condition = dist != 1)

  message('Computing longest hop run')
  values(gwalks)$longest.hoprun = grl.eval(gwalks,
                                         max(c(0, table(label.runs(abs(dist.nostrand)>hop.thresh)))),
                                         condition = dist != 1)

  message('Computing longest big hop run 1e5')
  values(gwalks)$longest.bhoprun1e5 = grl.eval(gwalks,
                                             max(c(0, table(label.runs(abs(dist.nostrand)>bighop.thresh & frag.pos<1e5)))),
                                             condition = dist != 1)

  message('Computing longest big hop run 1e4')
  values(gwalks)$longest.bhoprun1e4 = grl.eval(gwalks,
                                             max(c(0, table(label.runs(abs(dist.nostrand)>bighop.thresh & frag.pos<1e4)))),
                                             condition = dist != 1)


  message('Computing longest big hop run 1e3')
  values(gwalks)$longest.bhoprun1e3 = grl.eval(gwalks,
                                             max(c(0, table(label.runs(abs(dist.nostrand)>bighop.thresh & frag.pos<1e3)))),
                                             condition = dist != 1)

  message('Computing number of big 4 hop runs 1e4')
  values(gwalks)$num4bhoprun1e4 = grl.eval(gwalks,
                                         sum(table(label.runs(abs(dist.nostrand)>bighop.thresh & frag.pos<1e4))>3),
                                         condition = dist != 1)


  message('Computing number of big 4 hop runs 1e5')
  values(gwalks)$num4bhoprun1e5 = grl.eval(gwalks,
                                         sum(table(label.runs(abs(dist.nostrand)>bighop.thresh & frag.pos<1e5))>3),
                                         condition = dist != 1)


  message('Computing longest big hop run')
  values(gwalks)$longest.bhoprun = grl.eval(gwalks,
                                          max(c(0, table(label.runs(abs(dist.nostrand)>bighop.thresh)))),
                                          condition = dist != 1)


  message('Computing palindromic fraction')
  values(gwalks)$palindromic.frac = grl.eval(gwalks, sum((as.numeric(width)*sign(is.dup(paste(seqnames, start, end)))))/sum(as.numeric(width)))

  message('Computing palindromic fraction for 1e5 or smaller')
  values(gwalks)$palindromic.frac1e5 = grl.eval(gwalks, sum((as.numeric(frag.pos)*sign(is.dup(paste(seqnames, start, end)))))/sum(as.numeric(frag.pos)), dist!=1 & frag.pos<1e5)


  message('Computing palindromic width')
  values(gwalks)$palindromic.width1e5 = grl.eval(gwalks, sum((as.numeric(frag.pos)*sign(is.dup(paste(seqnames, start, end))))), dist!=1 & frag.pos<1e5)

                                        #    values(gwalks)$palindromic.width1e4 = grl.eval(gwalks, sum((as.numeric(frag.pos)*sign(is.dup(paste(seqnames, start, end))))), dist!=1 & frag.pos<1e4)

                                        #   values(gwalks)$palindromic.width1e3 = grl.eval(gwalks, sum((as.numeric(frag.pos)*sign(is.dup(paste(seqnames, start, end))))), dist!=1 & frag.pos<1e3)


  message('computing 1e4 pad footprint')
  values(gwalks)$numloci1e4 = elementNROWS(grl.reduce(gwalks, 1e4))

  message('computing 1e5 pad footprint')
  values(gwalks)$numloci1e5 = elementNROWS(grl.reduce(gwalks, 1e5))

  message('computing 1e6 pad footprint')
  values(gwalks)$numloci1e6 = elementNROWS(grl.reduce(gwalks, 1e6))

  message('computing 1e7 pad footprint')
  values(gwalks)$numloci1e7 = elementNROWS(grl.reduce(gwalks, 1e7))

  message('computing 1e8 pad footprint')
  values(gwalks)$numloci1e8 = elementNROWS(grl.reduce(gwalks, 1e8))

  message('computing 1e9 pad footprint')
  values(gwalks)$numloci1e9 = elementNROWS(grl.reduce(gwalks, 1e9, clip = TRUE))

  return(gwalks)
}

#' @name anno.hops
#' @title anno.hops
#' @description
#'
#' Adds simple annotations to GRangesList of walks including
#' distance along each reference fragment and distance
#' between "hops"
#' 
#' @export
#' @param walks walks to annotate
#' 
#' @author Marcin Imielinski
anno.hop = function(walks)
{
  gw = gr2dt(grl.unlist(walks))

  if (is.null(gw$ab.id))
  {
    ## mark nodes that precede a reference junction
    gw$ab.id = as.integer(NA)
    gw[, d.to.next := c((start-data.table::shift(end))[-1], NA), by = pid]
    gw[, d.to.next.neg := c((data.table::shift(start)-end)[-1], NA), by = pid]
    gw[, same.strand := c((strand==data.table::shift(strand))[-1], NA), by = pid]
    gw[, same.chrom := c((as.character(seqnames)==data.table::shift(as.character(seqnames)))[-1], NA), by = pid]
    gw[, last.node := 1:.N == .N, by = pid]
    gw[, before.ref :=
                 (((d.to.next<=1 & d.to.next>=0 & strand == '+') |
                   (d.to.next.neg<=1 & d.to.next.neg>=0 & strand == '-')
                 ) & same.strand & same.chrom)]
    gw[is.na(before.ref), before.ref := FALSE]
    
    ## label reference runs of nodes then collapse
    .labrun = function(x) ifelse(x, cumsum(diff(as.numeric(c(FALSE, x)))>0), as.integer(NA))
    gw[, ref.run := .labrun(before.ref), by = pid]
    gw[, ref.run.last := data.table::shift(ref.run), by = pid]
    gw[is.na(ref.run) & !is.na(ref.run.last), ref.run := ref.run.last]
    gw[!is.na(ref.run), ref.run.id := paste(pid, ref.run)]
    gw[is.na(ref.run), ab.id := 1:.N]
  }

  gw[, ab.chunk := cumsum(!is.na(ab.id)),  by = grl.ix]
  gw[, dist := c(ifelse((seqnames[-1] != seqnames[-length(seqnames)]) |
                        (strand[-1] != strand[-length(strand)]), Inf,
                 ifelse(strand[-1]=='+',
                        start[-1]-end[-length(end)],
                        start[-length(end)]-end[-1]))
               , Inf), by = grl.ix]

  gw[, dist.nostrand := c(ifelse((seqnames[-1] != seqnames[-length(seqnames)]), Inf,
                          ifelse(strand[-1]=='+',
                                 start[-1]-end[-length(end)],
                                 start[-length(end)]-end[-1]))
                        , Inf), by = grl.ix]

  gw[, dist := ifelse((1:length(grl.iix) %in% length(grl.iix)), as.numeric(NA), dist), by = grl.ix]
  gw[, dist.nostrand := ifelse((1:length(grl.iix) %in% length(grl.iix)), as.numeric(NA), dist),
     by = grl.ix]
  
  ## gw[, frag.id := as.character(frag.id)]
  gw[, ":="(frag.id = paste(grl.ix, ab.chunk), 
            frag.iid = 1:length(grl.ix),
            frag.pos = cumsum(width)
            ), by = .(ab.chunk, grl.ix)]

  gr.out = dt2gr(gw)
  gr.out$width = NULL

  setkey(gw, frag.id)

  grl.out = split(gr.out[, c('ab.id','grl.iix', 'cn', 'cn.1', 'frag.id', 'frag.iid', 'frag.pos', 'dist', 'dist.nostrand')], gr.out$grl.ix)[as.character(1:length(walks))]

  names(grl.out) = names(walks)
  values(grl.out) = values(walks)

    return(grl.out)
}



####################################################################
#' ppgrid
#'
#' least squares grid search for purity and ploidy modes
#'
#' @param segstats GRanges object of intervals with meta data fields "mean" and "sd" (i.e. output of segstats function)
#' @param allelic logical flag, if TRUE will also look for mean_high, sd_high, mean_low, sd_low variables and choose among top solutions from top copy number according to the best allelic fit
#' @param mc.cores integer number of cores to use (default 1)
#' @return data.frame with top purity and ploidy solutions and associated gamma and beta values, for use in downstream jbaMI
############################################
ppgrid = function(
  segstats, # n x 1 GRanges object with "mean" and "sd" value fields, optional field $ncn for "normal tissue" cn (default = 2)

  ########### optional args to describe the "valid modes"
    allelic = FALSE, ## if TRUE will also look for mean_high, sd_high, mean_low, sd_low variables and choose among top solutions from top copy number according to the best allelic fit
    purity.min = 0.01,
    purity.max = 1.0,
    ploidy.step = 0.01,
    purity.step = 0.01,
    ploidy.min = 1.2, # ploidy bounds (can be generous)
    ploidy.max = 6,
    plot = F,
    verbose = F,
    mc.cores = 10
  )
{

  if (verbose)
      message('setting up ppgrid matrices .. \n')

    if (is.na(ploidy.min)) ploidy.min = 1.2
    if (is.na(ploidy.max)) ploidy.max = 6
    if (is.na(purity.min)) purity.min = 0.01
    if (is.na(purity.max)) purity.max = 1

    ##  purity.guesses = seq(0, 1, purity.step)
    purity.guesses = seq(pmax(0, purity.min), pmin(1.00, purity.max), purity.step)
    ## ploidy.guesses = seq(pmin(0.5, ploidy.min), pmax(10, ploidy.max), ploidy.step)
    ploidy.guesses = seq(pmax(0.5, ploidy.min), pmax(0.5, ploidy.max), ploidy.step)

    if (allelic)
        if (!all(c('mean_high', 'mean_low', 'sd_high', 'sd_low') %in% names(values(segstats))))
        {
            warning('If allelic = TRUE then must have meta data fields mean_high, mean_low, sd_high, sd_low in input segstats')
            allelic = FALSE
        }

    if (is.null(segstats$mean))
        stop('segstats must have field $mean')

    segstats = segstats[!is.na(segstats$mean) & !is.na(segstats$sd)]

    if (!is.null(segstats$ncn))
        segstats = segstats[segstats$ncn==2, ]

    ## if (is.null(segstats$ncn))
    ##     ncn = rep(2, length(mu))
    ## else
    ##     ncn = segstats$ncn

    mu = segstats$mean
    w = as.numeric(width(segstats))
    Sw = sum(as.numeric(width(segstats)))
    sd = segstats$sd
    m0 = sum(as.numeric(mu*w))/Sw

    if (verbose)
        cat(paste(c(rep('.', length(purity.guesses)), '\n'), collapse = ''))

    NLL = matrix(unlist(mclapply(1:length(purity.guesses), function(i)
    {
        if (verbose)
            cat('.')
        nll = rep(NA, length(ploidy.guesses))
        for (j in 1:length(ploidy.guesses))
        {
            alpha = purity.guesses[i]
            tau = ploidy.guesses[j]
            gamma = 2/alpha - 2
            beta = (tau + gamma)/m0 ## replaced with below 9/10/14
                                        #          beta = ( tau + tau_normal * gamma /2 ) / m0
                                        #          v = pmax(0, round(beta*mu-ncn*gamma/2))
            v = pmax(0, round(beta*mu-gamma))

            ## using normal cn
                                        #          nll[j] = sum((v-beta*mu+ncn*gamma/2)^2/((sd)^2))

                                        # REVERT
            nll[j] = sum((v-beta*mu+gamma)^2/((sd)^2))

                                        #  mv = pmax(20, pmin(20, max(v, na.rm = TRUE)))
                                        #  mv = 2
                                        # log likelihood matrix across approximately "all" integer copy states
                                        # we are obviously ignoring very high states in this estimate
                                        # but they are likely to have high sd and thus contribute less to the overall log likelihood
            ##           ll = -sapply(0:mv, function(x) (x-beta*mu+gamma)^2/((sd)^2)) ## OG VERSION
                                        #    ll = -sapply(0:mv, function(x) ((x+gamma)/beta-mu)^2 / sd^2)  ## NEWFANGLED VERSION
            ##        ml = apply(ll, 1, max) ##  get maximum likelihood
            ##       probs  = 1/rowSums(exp(ll-ml)) ## normalize to get posterior probabilities (assuming uniform prior)
            ##      nll[j] = sum(-log(probs))
        }
        return(nll)
    }, mc.cores = mc.cores)), nrow = length(purity.guesses), byrow = T)

    dimnames(NLL) = list(as.character(purity.guesses), as.character(ploidy.guesses))

  if (verbose)
    cat('\n')

    ## rix = as.numeric(rownames(NLL))>=purity.min & as.numeric(rownames(NLL))<=purity.max
    ## cix = as.numeric(colnames(NLL))>=ploidy.min & as.numeric(colnames(NLL))<=ploidy.max
    ## NLL = NLL[rix, cix, drop = FALSE]

    a = rep(NA, nrow(NLL));
    b = rep(NA, ncol(NLL)+2)
    b.inf = rep(Inf, ncol(NLL)+2)
    #  a = rep(Inf, nrow(NLL));
    #  b = rep(Inf, ncol(NLL)+2)
    NLLc = rbind(b, cbind(a, NLL, a), b) ## padded NLL and all of its shifts
    NLLul = rbind(cbind(NLL, a, a), b.inf, b)
    NLLuc = rbind(cbind(a, NLL, a), b.inf, b)
    NLLur = rbind(cbind(a, a, NLL), b.inf, b)
    NLLcl = rbind(b, cbind(NLL, a, a), b)
    NLLcr = rbind(b, cbind(a, a, NLL), b)
    NLLll = rbind(b, b, cbind(NLL, a, a))
    NLLlc = rbind(b, b, cbind(a, NLL, a))
  NLLlr = rbind(b, b, cbind(a, a, NLL))

  NLLm = melt(NLL) %>% as.data.table %>% setnames(c('purity', 'ploidy', 'NLL'))

  NLLm[, purity.id := factor(purity, purity.guesses) %>% as.integer]
  NLLm[, ploidy.id := factor(ploidy, ploidy.guesses) %>% as.integer]

  setkeyv(NLLm, c('purity.id', 'ploidy.id'))
  NLLm$NLL.lc = NLLm[.(NLLm$purity.id-1, NLLm$ploidy.id), NLL]
  NLLm$NLL.rc = NLLm[.(NLLm$purity.id+1, NLLm$ploidy.id), NLL]
  NLLm$NLL.lu = NLLm[.(NLLm$purity.id-1, NLLm$ploidy.id+1), NLL]
  NLLm$NLL.ru = NLLm[.(NLLm$purity.id+1, NLLm$ploidy.id+1), NLL]
  NLLm$NLL.ld = NLLm[.(NLLm$purity.id-1, NLLm$ploidy.id-1), NLL]
  NLLm$NLL.rd = NLLm[.(NLLm$purity.id+1, NLLm$ploidy.id-1), NLL]
  NLLm$NLL.cd = NLLm[.(NLLm$purity.id, NLLm$ploidy.id-1), NLL]
  NLLm$NLL.cu = NLLm[.(NLLm$purity.id, NLLm$ploidy.id+1), NLL]
  NLLm[, minima := NLL <= pmin(NLL.lc, NLL.rc, NLL.lu, NLL.ru, NLL.ld, NLL.rd, NLL.cd, NLL.cu, na.rm = TRUE)]
  out = NLLm[minima == TRUE, .(purity, ploidy, NLL, i = purity.id, j = ploidy.id)]
  

    ## if (min(c(ncol(NLL), nrow(NLL)))>1) ## up up down down left right left right ba ba start
    ##     M = (NLLc < NLLul & NLLc < NLLuc & NLLc < NLLur & NLLc < NLLcl & NLLc < NLLcr & NLLc < NLLll & NLLc < NLLlc & NLLc < NLLlr)[-c(1, nrow(NLLc)), -c(1, ncol(NLLc)), drop = FALSE]
    ## else if (ncol(NLL)==1) ## one column, only go up and down
    ##     M = (NLLc < NLLuc & NLLc < NLLlc)[-c(1, nrow(NLLc)), -c(1, ncol(NLLc)), drop = FALSE]
    ## else ## only row, only go left right
    ##     M = (NLLc < NLLcl & NLLc < NLLcr)[-c(1, nrow(NLLc)), -c(1, ncol(NLLc)), drop = FALSE]

    ## if (length(M)>1)
    ## {
    ##     ix = Matrix::which(M, arr.ind= T);
    ##     if (nrow(ix)>1)
    ##     {
    ##         C = hclust(d = dist(ix), method = 'single')
    ##         cl = cutree(C, h = min(c(nrow(NLL), ncol(NLL), 2)))
    ##         minima = ix[vaggregate(1:nrow(ix), by = list(cl), function(x) x[which.min(NLL[ix[x, drop = FALSE]])]), , drop = FALSE]
    ##     }
    ##     else
    ##         minima = ix[1,, drop = FALSE]
    ## }
    ## else
    ##     minima = cbind(1,1)

    ## out = data.frame(purity = as.numeric(rownames(NLL)[minima[,1]]), ploidy = as.numeric(colnames(NLL)[minima[,2]]), NLL = NLL[minima],
    ##                  i = minima[,1], j = minima[,2])

    out = out[order(out$NLL), , drop = FALSE]
    rownames(out) = 1:nrow(out)
    ## Saturday, Sep 02, 2017 10:33:26 PM
    ## Noted floating point error, use the epsilon trick to replace '>='
    ## out = out[out$purity>=purity.min & out$purity<=purity.max & out$ploidy>=ploidy.min & out$ploidy<=ploidy.max, ]
    eps = 1e9
    out = out[out$purity - purity.min >= -eps &
              out$purity - purity.max <= eps &
              out$ploidy - ploidy.min >= -eps &
              out$ploidy - ploidy.max <= eps, ]
    out$gamma = 2/out$purity -2
    out$beta = (out$ploidy + out$gamma)/m0
    out$mincn = mapply(function(gamma, beta) min(round(beta*mu-gamma)), out$gamma, out$beta)
    out$maxcn = mapply(function(gamma, beta) max(round(beta*mu-gamma)), out$gamma, out$beta)

    ## group solutions with (nearly the same) slope (i.e. 1/beta), these should have almost identical
    ## NLL (also take into account in-list distance just be safe)
    if (nrow(out)>1)
        out$group = cutree(hclust(d = dist(cbind(100/out$beta, 1:nrow(out)), method = 'manhattan'), method = 'single'), h = 2)
    else
        out$group = 1
    out = out[out$group<=3, ,drop = FALSE] ## only pick top 3 groups

    if (allelic) ## if allelic then use allelic distance to rank best solution in group
    {
        ## remove all NA allelic samples
        segstats = segstats[!is.na(segstats$mean_high) & !is.na(segstats$sd_high) & !is.na(segstats$mean_low) & !is.na(segstats$sd_low)]
        out$NLL.allelic = NA
        mu = cbind(segstats$mean_high, segstats$mean_low)
        w = matrix(rep(as.numeric(width(segstats)), 2), ncol = 2, byrow = TRUE)
        Sw = sum(as.numeric(width(segstats)))*2
        sd = cbind(segstats$sd_high, segstats$sd_low)
        m0 = sum(as.numeric(mu*w))/Sw

        if (verbose)
            cat(paste(c(rep('.', length(purity.guesses)), '\n'), collapse = ''))

        for (i in 1:nrow(out))
        {
          if (verbose)
            {
              message(sprintf('Evaluating alleles for solution %s of %s\n', i, nrow(out)))
            }
            alpha = out$purity[i]
            tau = out$ploidy[i]
                                        #                  gamma = 2/alpha - 2
            gamma = 1/alpha - 1 ## 1 since we are looking at hets
            beta = (tau + gamma)/m0 ## replaced with below 9/10/14
                                        #          beta = ( tau + tau_normal * gamma /2 ) / m0
                                        #          v = pmax(0, round(beta*mu-ncn*gamma/2))
            v = pmax(0, round(beta*mu-gamma))

            vtot = round(out$beta[i]*segstats$mean-out$gamma[i])
            vlow.mle = rep(NA, length(vtot))

            for (j in 1:length(vlow.mle))
            {
                if (vtot[j]==0)
                    vlow.mle[j] = 0
                else
                {
                    vlow = 0:floor(vtot[j]/2)
                    vhigh = vtot[j]-vlow
                    tmp.nll = cbind((vlow-beta*mu[j,2]+gamma)^2/(sd[j,2])^2, (vhigh-beta*mu[j, 1]+gamma)^2/((sd[j,1])^2))
                    vlow.mle[j] = vlow[which.min(rowSums(tmp.nll))]
                }
            }

            vlow.mle = apply(cbind(mu, sd, vtot), 1, function(x) {
                tot = x[5]
                if (tot == 0)
                    return(0)
                else
                {
                    vlow = 0:floor(tot/2)
                    vhigh = tot-vlow
                    muh = x[1]
                    mul = x[2]
                    sdh = x[3]
                    sdl = x[4]
                    tmp.nll = cbind((vlow-beta*mul+gamma)^2/(sdl)^2, (vhigh-beta*muh+gamma)^2/((sdh)^2))
                    return(vlow[which.min(rowSums(tmp.nll))])
                }
            })

            out$NLL.allelic[i] = sum((cbind(vtot-vlow.mle, vlow.mle)-beta*mu+gamma)^2/sd^2)
        }

        out$NLL.tot = out$NLL
        out$NLL = out$NLL.tot + out$NLL.allelic
        out.all = out
        ix = vaggregate(1:nrow(out), by = list(out$group), FUN = function(x) x[order(abs(out$NLL[x]))][1])
    }
    else ## otherwise choose the one that gives the lowest magnitude copy number
    {
        out.all = out
        ix = vaggregate(1:nrow(out), by = list(out$group), FUN = function(x) x[order(abs(out$mincn[x]), out$mincn[x]<0)][1])
    }

    out = out[ix, , drop = FALSE]
    out$NLL = vaggregate(out$NLL, by = list(out$group), FUN = min)

    out.all$keep = 1:nrow(out.all) %in% ix ## keep track of other ploidy group peaks for drawing purposes
    out.all = out.all[out.all$group %in% out$group, ] ## only draw the groups in the top solution

    if (plot)
    {
        library(ellipse)
        library(numDeriv)

        ## xval = as.numeric(rownames(NLL))
        ## yval = as.numeric(colnames(NLL))
        ## f = function(x) {
        ##     i = x[1]; ## interpolate between closest values of NLL
        ##     im = which(i<=xval)[1]
        ##     ip = (i-xval[im-1])/diff(xval[c(im-1, im)]);  ## proportion of lower value to consider
        ##     j = x[2];
        ##     jm = which(j<=yval)[1]
        ##     jp = (j-yval[jm-1])/diff(yval[c(jm-1, jm)]);  ## proportion of lower value to consider
        ##     nllm = NLL[c(im-1, im), c(jm-1, jm)] ## piece of NLL matrix containing the low and high i and j matches
        ##     nllp = cbind(c(ip, 1-ip)) %*% rbind(c(jp, 1-jp)) ## proportion of values to input into interpolation
        ##     return(sum(-nllm*nllp))
        ## }

        plot(out.all$purity, out.all$ploidy, pch = 19,
             xlim = c(purity.min, purity.max), ylim = c(ploidy.min, ploidy.max), xlab = 'purity', ylab = 'ploidy', cex = 0.2, col = alpha('white', out.all$intensity))


        f = function(x) -NLL[((x[1]-1) %% nrow(NLL))+1, ((x[2]-1) %% ncol(NLL))+1]
          ir = range(as.numeric(rownames(NLL)))
          jr = range(as.numeric(colnames(NLL)))
          txf = function(z) cbind(affine.map(z[,1], ir, c(1, nrow(NLL))), affine.map(z[,2], jr, c(1, ncol(NLL))))

          levs = c(0.95, 0.99)
                                        #          levs = c(1-1e-7)
          tmp.out = out.all
          tmp.out$NLL = as.data.table(tmp.out)[, NLL := min(NLL), by = group][, NLL]
          tmp.out$intensity = affine.map(tmp.out$NLL, c(1, 0.1))
          tmp.out = tmp.out[rev(1:nrow(tmp.out)), ]
                                        #          tmp.out$col = brewer.master(length(levels(tmp.out$group)), 'PuRd')[as.integer(tmp.out$group)]
          tmp.out$col = brewer.pal(length(unique(out$group))+2, 'PuRd')[match(tmp.out$group, unique(tmp.out$group))]
          tmp.out$rank = ''; ## hacky stuff to just plot ranks for the top per group solutions
          tmp.out$rank[tmp.out$keep] = match(tmp.out$group[tmp.out$keep], out$group)

          require(DiceKriging)
          bla = mapply(function(x, y, c, a)
          {
              ## grab k square from computed NLL values to krig around
              k = 4
              ir = pmin(nrow(NLL), pmax(1, (x-k):(x+k)))
              jr = pmin(ncol(NLL), pmax(1, (y-k):(y+k)))
              ij = expand.grid(ir, jr)
              xy = expand.grid(purity.guesses[ir], ploidy.guesses[jr])
              m = tryCatch(km(design = xy[, 1:2], response = -NLL[as.matrix(ij)]), error = function(e) NULL)
              ## custom function gives best krigged interpolation in
              ## region will be used for hessian computaiton
              f = function(x) predict(m, newdata = data.frame(Var1 = x[1], Var2 = x[2]), type = 'UK')$mean
              for (lev in levs) ## TOFIX: MESS
              {
                  if (!is.null(m))
                      F = tryCatch(-hessian(f, c(purity.guesses[x], ploidy.guesses[y])), error = function(e) matrix(NA, ncol = 2, nrow = 2))
                  else
                      F = NA

                  if (all(!is.na(F)))
                      V = solve(F)
                  else
                      V = NA
                  if (any(is.na(V)))
                      V = diag(rep(1,2))
                  M = cov2cor(V)
                  XY = ellipse(M, scale = .01*c(diff(par('usr')[1:2]),diff(par('usr')[3:4])), centre = c(purity.guesses[x], ploidy.guesses[y]), level = lev)
                  polygon(XY[,1], XY[,2], border = NA, col = alpha(c, a*affine.map(lev, c(1, 0.3))));
              }
          }, tmp.out$i, tmp.out$j, tmp.out$col, tmp.out$intensity, SIMPLIFY = FALSE)

          tmp.out = tmp.out[tmp.out$keep, ]
          text(tmp.out$purity, tmp.out$ploidy, tmp.out$rank, col = alpha('black', 0.5))

          tmp.out = tmp.out[rev(1:nrow(tmp.out)), ]
          legend(0, par('usr')[4]*0.98, legend = paste(tmp.out$rank, ') ', sapply(tmp.out$purity, function(x) sprintf('%0.2f',x)), ', ',
                                                       sapply(tmp.out$ploidy, function(x) sprintf('%0.2f',x)),
                                                       ' (gamma = ',sapply(tmp.out$gamma, function(x) sprintf('%0.3f',x)),
                                                       ', beta = ',sapply(tmp.out$beta, function(x) sprintf('%0.3f',x)),
                                                       ')', sep = ''), fill = tmp.out$col, cex = 0.8, yjust = 1, ncol = 1)
      }

    out = out.all;
    out = out[order(out$group, !out$keep, out$NLL), ]
    out$rank = NA
    out$rank[out$keep] = 1:sum(out$keep)
    out$keep = out$i = out$j = NULL
    rownames(out) = NULL
    return(out)
}

#' Identify matches between query and dictionary
#'
#' Wrapper around matchPdict to identify matches between a query
#' string query and dictionary dict (both BString objects or subclasses)
#'
#' @param query Query XStringSet / DNAStringSet 
#' @param dict Dictionary
#' @param midpoint Flag for output the coordinates of the match as the location,
#'   where the midpoint of the dict string matches the given query. Default FALSE
#' @return a vector of indices of length width(query) that contains
#' indices of the (starting) dictionary in the query string
match.bs = function(query, dict, midpoint = FALSE)
{
  names(dict) = as.character(1:length(dict))

  tmp = sort(unlist(Biostrings::matchPDict(dict, query)))
  out = rep(NA, length(query))

  if (!midpoint)
    out[start(tmp)] = as.numeric(names(tmp))
  else
    out[floor((start(tmp)+end(tmp))/2)] = as.numeric(names(tmp))

  return(out)
}

## not exported in dev branch, gUtils
## #' @name grl.stripnames
## # ' @title Remove \code{GRanges} names inside a \code{GRangesList}
## #' @description
## #'
## #' Remove \code{GRanges} names inside a \code{GRangesList}
## #'
## #' @param grl \code{GRangesList} with names elements
## #' @return \code{GRangesList} where \code{GRanges} have no names
grl.stripnames = function(grl)
{
    ele = tryCatch(as.data.frame(grl)$element, error = function(e) NULL)
    if (is.null(ele)){
        ele = unlist(lapply(1:length(grl), function(x) rep(x, length(grl[[x]]))))
    }

    gr = unlist(grl);
    names(gr) = NULL;

    out = split(gr, ele);
    values(out) = values(grl)
    names(out) = names(grl)

    return(out)
}




#' @name hets
#' @export
#' @title Simple het "caller" meant to be used at validated het SNP sites for tumor / normal pairs
#' @description 
#'
#' hets() outputs a tsv file of ALT ($alt.count.t, $alt.count.n) and REF ($ref.count.t,, $ref.count.n) read counts to out.file
#' for a tumor / normal pair across a set of sites specified by an input VCF
#'
#' @param tum.bam string path to tumor sample, input to Bamfile()
#' @param norm.bam string path to normal sample, input to Bamfile()(optional) (default = NULL)
#' @param out.file string path to TSV output file to be generated 
#' @param vcf.file string path to VCF file of sites (eg hapmap or 1000G) at which to compute read counts
#' @param chunk.size1 integer Number of variants to process from VCF file at a time (default = 1e3)
#' @param chunk.size2 integer Number of variants to access from BAM file in a single iteration (default = 1e2)
#' @param mc.cores integer Number of cores in mclapply (default = 1)
#' @param verbose boolean Flag to increase verbosity (default = TRUE)
#' @param na.rm logical Flag to remove rows with NA counts (default = TRUE)
#' @param filt.norm logical Flag to remove any sites that have allele fraction of 0 or 1 or NA in MAF; if TRUE will remove any sites that have allele fraction 0 or 1 or NA in MAF 
#' @return nil
#' @author Marcin Imielinski
#' @export
hets = function(tum.bam, norm.bam = NULL, out.file, vcf.file, chunk.size1 = 1e3, chunk.size2 = 1e2, mc.cores = 1, verbose = TRUE, na.rm = TRUE, filt.norm = TRUE)
{    
    f = file(vcf.file, 'r')
      
    if (grepl('VCF', readLines(f, 1))){
        vcf = TRUE
    }
    else{
        vcf = FALSE
    }

    sl = hg_seqlengths()

    if (verbose){
        st = Sys.time()
    }

    nprocessed = 0
    nhets = 0
    first = TRUE
    ## get past headers

    ## while (grepl('^#', last.line <<- readLines(f, n=1))){}

    if (verbose){
        cat('Opened vcf, writing hets to text file', out.file, '\n')
    }

    out.cols = c('seqnames', 'start', 'end', 'Tumor_Seq_Allele1', 'Reference_Allele', 'ref.count.t', 'alt.count.t', 'ref.count.n', 'alt.count.n', 'alt.frac.t', 'ref.frac.t', 'alt.frac.n', 'ref.frac.n')

    if (vcf){
        col.ix = 1:5
    }
    else{
        col.ix = match(c("Chromosome", "Start_position", "End_position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"), strsplit(last.line, '\t')[[1]])
        if (any(is.na(col.ix))){
            stop('Error: failure processing variant file: must be valid VCF or MAF')
        }
    }
      
    while (!is.null(tmp <- tryCatch(read.delim(file = f, as.is = T, header = F, nrows = chunk.size1)[, col.ix], error = function(x) NULL))){
        
        if (vcf){
            names(tmp) = c('chr', 'start', 'name', 'ref', 'alt')
        }
        else{
            names(tmp) = c('chr', 'start', 'name', 'ref', 'alt', 'alt2')
            ## just in case the first tumor seq allele is equal to reference .. which happens in mafs
            tmp$alt = ifelse(tmp$alt==tmp$ref, tmp$alt2, tmp$alt)
        }
              
        loc = seg2gr(tmp, seqlengths = sl)    
        clock({loc.count = mafcount(tum.bam, norm.bam, loc, indel = T, chunk.size = chunk.size2, mc.cores = mc.cores)})
        nprocessed = nprocessed + length(loc.count)
              
        if (filt.norm & !is.null(loc.count$alt.frac.n)){
            loc.count = loc.count[which(loc.count$alt.frac.n != 1 & loc.count$alt.frac.n != 0)]
        }
              
        nhets = nhets + length(loc.count)
        if (length(loc.count)>0){

            df = as.data.frame(loc.count)
            ## remove any entries with 0 ref or alt reads in tumor or normal
            if (na.rm){
                if (!is.null(norm.bam)){
                    naix = apply(df[, c('alt.count.t', 'ref.count.t', 'alt.count.n', 'ref.count.n')], 1, function(x) all(is.na(x)))
                }
                else{
                    naix = apply(df[, c('alt.count.t', 'ref.count.t')], 1, function(x) all(is.na(x)))
                }
                df = df[which(!naix), ]
            }

            out.cols = intersect(out.cols, names(df))

            if (first){
                write.tab(df[, out.cols], out.file, append = F, col.names = T)
                first = F
            }
            else{
                write.tab(df[, out.cols], out.file, append = T, col.names = F)
            }                     
              
            if (verbose){
                cat(sprintf('Processed %s sites, wrote %s candidate hets\n', nprocessed, nhets))
            }

            if (verbose){
                cat('Time elapsed:\n')
                print(Sys.time() - st)
            }              
        }
    }
      
    close(f)
     
    if (verbose){
        cat('Finished het processing wrote to file', out.file, '\n')
    }
}











#' Merges rearrangements represented by \code{GRangesList} objects
#'
#' Determines overlaps between two or more piles of rearrangement junctions (as named or numbered arguments) +/- padding
#' and will merge those that overlap into single junctions in the output, and then keep track for each output junction which
#' of the input junctions it was "seen in" using logical flag  meta data fields prefixed by "seen.by." and then the argument name
#' (or "seen.by.ra" and the argument number)
#'
#' @param ... GRangesList representing rearrangements to be merged
#' @param pad non-negative integer specifying padding
#' @param ind  logical flag (default FALSE) specifying whether the "seen.by" fields should contain indices of inputs (rather than logical flags) and NA if the given junction is missing
#' @param ignore.strand whether to ignore strand (implies all strand information will be ignored, use at your own risk)
#' @return \code{GRangesList} of merged junctions with meta data fields specifying which of the inputs each outputted junction was "seen.by"
#' @name ra.merge
#' @export
#' @examples
#'
#' # generate some junctions
#' gr1 <- GRanges(1, IRanges(1:10, width = 1), strand = rep(c('+', '-'), 5))
#' gr2 <- GRanges(1, IRanges(4 + 1:10, width = 1), strand = rep(c('+', '-'), 5))
#' ra1 = split(gr1, rep(1:5, each = 2))
#' ra2 = split(gr2, rep(1:5, each = 2))
#'
#' ram = ra.merge(ra1, ra2)
#' values(ram) # shows the metadata with TRUE / FALSE flags
#'
#' ram2 = ra.merge(ra1, ra2, pad = 5) # more inexact matching results in more merging
#' values(ram2)
#'
#' ram3 = ra.merge(ra1, ra2, ind = TRUE) #indices instead of flags
#' values(ram3)
#' @export
ra.merge = function(..., pad = 0, ind = FALSE, ignore.strand = FALSE){
    ra = list(...)
    ra = ra[which(!sapply(ra, is.null))]
    nm = names(ra)
    if (is.null(nm)){
        nm = paste('ra', 1:length(ra), sep = '')
    }
    nm = paste('seen.by', nm, sep = '.')
    if (length(nm)==0){
        return(NULL)
    }
    out = ra[[1]]
    values(out) = cbind(as.data.frame(matrix(FALSE, nrow = length(out), ncol = length(nm), dimnames = list(NULL, nm))), values(out))

    if (!ind){
        values(out)[, nm[1]] = TRUE
    } else{
        values(out)[, nm[1]] = 1:length(out)
    }

    if (length(ra)>1){
        for (i in 2:length(ra)){
            this.ra = ra[[i]]
            if (length(this.ra)>0){
                values(this.ra) = cbind(as.data.frame(matrix(FALSE, nrow = length(this.ra), ncol = length(nm), dimnames = list(NULL, nm))), values(this.ra))
                ovix = ra.overlaps(out, this.ra, pad = pad, ignore.strand = ignore.strand)

                if (!ind){
                    values(this.ra)[[nm[i]]] = TRUE
                } else{
                    values(this.ra)[[nm[i]]] = 1:length(this.ra)
                }

                if (!ind){
                    if (!all(is.na(ovix))){
                        values(out)[, nm[i]][ovix[,1]] = TRUE
                    }
                } else{
                    values(out)[, nm[i]] = NA
                    if (!all(is.na(ovix))){
                        values(out)[, nm[i]][ovix[,1]] = ovix[,1]
                    }
                }
                ## which are new ranges not already present in out, we will add these
                if (!all(is.na(ovix))){
                    nix = setdiff(1:length(this.ra), ovix[,2])
                } else{
                    nix = 1:length(this.ra)
                }

                if (length(nix)>0){
                    val1 = values(out)
                    val2 = values(this.ra)
                    if (ind){
                        val2[, nm[1:(i-1)]] = NA
                    }
                    else{
                        val2[, nm[1:(i-1)]] = FALSE
                    }
                    values(out) = NULL
                    values(this.ra) = NULL
                    out = grl.bind(out, this.ra[nix])
                    values(out) = rrbind(val1, val2[nix, ])
                } else if (length(setxor(1:length(this.ra), ovix[, 2])) == 0) { ## fix if there is perfect overlap... this case was previously not dealt with
                    merge_cols = setdiff(colnames(values(this.ra)), nm)
                    if (length(merge_cols) > 0) {
                        tmp_field = "out_ix_50169346127375946273"
                        val1 = values(out)
                        val2 = values(this.ra)
                        val1[[tmp_field]] = seq_along(out)
                        val2[ovix[,2],tmp_field] = ovix[,1]
                        new_val = merge(val1, val2[,c(tmp_field, merge_cols)], by = tmp_field, all.x = TRUE)
                        new_val = new_val[order(new_val[[tmp_field]]),]
                        new_val[[tmp_field]] = NULL
                        values(out) = new_val
                    }
                }
            }
        }
    }
    return(out)
}




#' @name jabba.gwalk
#' @export
#' @rdname internal
#' @title jabba.gwalk
#' @description
#'
#' Computes greedy collection (i.e. assembly) of genome-wide walks (graphs and cycles) by finding shortest paths in JaBbA graph.
#'
#' @param jab JaBbA object
#' #
#' @return GRangesList of walks with copy number as field $cn, cyclic walks denoted as field $is.cycle == TRUE, and $wid (width) and $len (segment length) of walks as additional metadata#' 
#' @import igraph
#' @author Marcin Imielinski
#' @author Xiaotong Yao
jabba.gwalk = function(jab, verbose = FALSE, return.grl = TRUE)
{
    cn.adj = jab$adj
    adj = as.matrix(cn.adj)
    adj.new = adj*0
    ## ALERT!!! see below
    adj[Matrix::which(adj!=0, arr.ind = TRUE)] = width(jab$segstats)[Matrix::which(adj!=0, arr.ind = TRUE)[,2]] ## make all edges a large number by default
    ## adj[which(adj!=0, arr.ind = TRUE)] = width(jab$segstats)[which(adj!=0, arr.ind = TRUE)[,1]] ## make all edges a large number by default
    if (verbose){
        ## ALERT!!! I'm gonna switch to source node width for default weight of edges
        jmessage('Setting edge weights to destination widths for reference edges and 1 for aberrant edges')
        ## jmessage('Setting default edge weights to SOURCE widths for edges and 1% less for aberrant edges')
    }

    ab.edges = rbind(jab$ab.edges[,1:2, 1], jab$ab.edges[,1:2, 2])
    ab.edges = ab.edges[Matrix::rowSums(is.na(ab.edges))==0, ]
    ## ALERT!!!
    adj[ab.edges] = sign(cn.adj[ab.edges]) ## make ab.edges = 1
    ## adj[ab.edges] = adj[ab.edges] * 0.99 ## make ab.edges 1 bp shorter than ref!
    adj[is.na(adj)] = 0
    cn.adj[which(is.na(cn.adj))] = 0

    ## ALERT!!! major change
    ## adjj = adj/as.matrix(cn.adj)
    ## adjj[which(is.nan(adjj))] = 0
    ## adjj[which(adjj<0)] = 0
    ## G = graph.adjacency(adjj, weighted = 'weight')
    ## esl = which(adj != 0, arr.ind=T)
    ## eids = paste(esl[,1], esl[,2])
    ## weights = adj[esl]
    ## eclasses = ed[.(eids), eclass]
    G = graph.adjacency(adj, weighted = 'weight')
    ## G = make_graph(t(esl), )

    ## DD = shortest.paths(G, mode="out")
    ## IJ = which(!is.infinite(DD), arr.ind=T)

    ## define ends not using degree (old method) but using either telomeres or loose ends
    ## (otherwise lots of fake ends at homozygous deleted segments)
    ss = gr2dt(jab$segstats)[ , vid:= 1:length(seqnames)]
    ss[loose == TRUE, is.end := TRUE]
    ss[loose == FALSE, is.end := 1:length(loose) %in% c(which.min(start), which.max(end)), by = list(seqnames, strand)]
    ends = which(ss$is.end)
    src = (Matrix::colSums(adj)[ends]==0) ## indicate which are sources

    ## sanity check
    unb = which(!ss$is.end & Matrix::rowSums(jab$adj, na.rm = TRUE) != Matrix::colSums(jab$adj, na.rm = TRUE))

    if (length(unb)>0)
    {
        jmessage(sprintf('JaBbA model not junction balanced at %s non-ends! Adding these to "ends"', length(unb)))
        ends = c(ends, unb)         ## shameless HACK ... TOFIX
    }

    ## ends = which(degree(G, mode = 'out')==0 | degree(G, mode = 'in')==0)
    i = 0
    ## adjust weight just before creating D
    ## assign lighter weight to higher copy
    ## D = shortest.paths(G, v = ends, mode = 'out', weight = E(G)$weight)[, ends]

    ## D records distance from ends to every node
    D = shortest.paths(G, v = ends, mode = 'out', weight = E(G)$weight)[, ends]

    ## sort shortest paths
    ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row != col, ][order(dist), ][, row := ends[row]][, col := ends[col]]

    ## ij only record end to end distance
    ## ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[col %in% ends, ][, dist := D[cbind(row, col)]][, row := ends[row]][row != col, ][order(dist), ]

    maxrow = length(ends)*max(cn.adj[ends, ends], na.rm = TRUE)
    vpaths = rep(list(NA), maxrow)
    epaths = rep(list(NA), maxrow)
    cns = rep(NA, maxrow)
    palindromic.path = rep(FALSE, maxrow)
    palindromic.cycle = rep(FALSE, maxrow)

    nb.all = which(Matrix::rowSums(cn.adj) != Matrix::colSums(cn.adj))
    cn.adj0 = cn.adj
    G0 = G
    D0 = D

    #' first peel off "simple" paths i.e. zero degree
    #' ends with >0 copy number
    psimp =  which(degree(G, mode = 'out')==0 & degree(G, mode = 'in')==0 & jab$segstats$cn>0)
    i = 0
    if (length(psimp)>0)
    {
        vpaths[1:length(psimp)] = split(psimp, 1:length(psimp))
        epaths[1:length(psimp)] = lapply(psimp, function(x) cbind(NA, NA)) ## there is no "edge" associated with a zero total degree node
        cns[1:length(psimp)] = jab$segstats$cn[psimp]
        i = length(psimp)
    }

    ## now iterate from shortest to longest path
    ## peel that path off and see if it is still there ..
    ## and see if it is still there
    ## peel off top path and add to stack, then update cn.adj

    jab$segstats$tile.id = jab$segstats$tile.id + as.numeric(jab$segstats$loose)*0.5

    tile.map =
        gr2dt(jab$segstats)[, .(id = 1:length(tile.id),
                                tile.id = ifelse(strand == '+', 1, -1)*tile.id)]
    rtile.map =
        gr2dt(jab$segstats)[, .(id = 1:length(tile.id),
                                tile.id = ifelse(strand == '+', 1, -1)*tile.id)]
    setkey(tile.map, id)
    setkey(rtile.map, tile.id)

    ## unique pair of edge ids: rev comp of a foldback edge will be identical to itself!!!
    ed = data.table(jab$edges)[cn>0, .(from, to , cn)]

    if (nrow(ed)==0)
        return(GRangesList())

    ed[, ":="(fromss = tile.map[ .(from), tile.id],
              toss = tile.map[ .(to), tile.id]),
       by = 1:nrow(ed)]
    ed[, weight :=  adj[cbind(from, to)]]
    print(ed)
    ed[fromss*toss > 0, eclass := ifelse(fromss>0, paste(fromss, toss), paste(-toss, -fromss))]
    ed[fromss*toss < 0, eclass := ifelse(abs(fromss)<=abs(toss),
                                         paste(fromss, toss), paste(-toss, -fromss))]
    ed[, eclass := as.numeric(as.factor(eclass))]
    ed[, eid := paste(from, to)]
    setkey(ed, "eid")
    eclass.cn = ed[!duplicated(eclass), setNames(cn, eclass)]

    cleanup_mode = FALSE


    while (nrow(ij)>0)
    {
        if (verbose)
            jmessage('Path peeling iteration ', i, ' with ', sum(adj!=0, na.rm = TRUE), ' edges left and ', nrow(ij), ' end-pairs to resolve' )
        i = i+1
        ## swap this
        ##        vpaths[[i]] = p = as.numeric(get.shortest.paths(G, ij[1, 1], ij[1, 2], mode = 'out', weight = E(G)$weight)$vpath[[1]])

        p = get.constrained.shortest.path(cn.adj, G, v = ij[1, 1], to = ij[1, 2], weight = E(G)$weight, edges = ed, verbose = TRUE, mip = cleanup_mode)

        if (is.null(p)){
            jmessage('Came up empty!')
            i = i -1
            ij = ij[-1, , drop = FALSE]
        }
        else
        {
            ## Don't forget to update ed here
            ed$cn = cn.adj[cbind(ed$from, ed$to)]

            vpaths[[i]] = p
            epaths[[i]] = cbind(p[-length(p)], p[-1])
            eids = paste(epaths[[i]][,1], epaths[[i]][,2])
            cns[i] = ed[.(eids), if (length(cn)>1) cn/2 else cn, by = eclass][, floor(min(V1))] ## update cn correctly, adjusting constraints for palinrdomic edges by 1/2

            rvpath = rtile.map[list(tile.map[list(vpaths[[i]]), -rev(tile.id)]), id]
            repath = cbind(rvpath[-length(rvpath)], rvpath[-1])
            plen = length(rvpath)
            hplen = floor(length(rvpath)/2)

            ## (awkward) check for palindromicity for odd and even length palindromes
            ## if (all((vpaths[[i]]==rvpath)[c(1:hplen,(plen-hplen+1):plen)]))
            if (ed[eids, any(table(eclass)>1)])
                palindromic.path[i] = TRUE
            ## else
            ## {
            vpaths[[i+1]] = rvpath
            epaths[[i+1]] = repath
            cns[i+1] = cns[i]
            palindromic.path[i+1] = TRUE
            ## }
            ##        palindromic = TRUE ## set to true while we "figure things out"


            #' so now we want to subtract that cn units of that path from the graph
            #' so we want to update the current adjacency matrix to remove that path
            #' while keeping track of of the paths on the stack
            cn.adj[epaths[[i]]] = cn.adj[epaths[[i]]]-cns[i]

            ## if (!palindromic) ## update reverse complement unless palindromic
            cn.adj[epaths[[i+1]]] = cn.adj[epaths[[i+1]]]-cns[i+1]

            if (!all(cn.adj[epaths[[i]]]>=0)) ## something wrong, backtrack
            {
                jmessage('backtracking ...') ## maybe we got stuck in a quasi-palindrome and need to backtrack
                                        #            browser()
                cn.adj[epaths[[i]]] = cn.adj[epaths[[i]]]+cns[i]
                ## if (!palindromic) ## update reverse complement unless palindromic
                cn.adj[epaths[[i+1]]] = cn.adj[epaths[[i+1]]]+cns[i+1]
                i = i-1
                ij = ij[-1, , drop = FALSE]
            }
            else ## continue, reduce
            {
                adj.new[epaths[[i]]] = adj.new[epaths[[i]]] + cns[i]
                ## if (!palindromic)
                adj.new[epaths[[i+1]]] = adj.new[epaths[[i+1]]] + cns[i]

                ## ## make sure I didn't overuse any edge
                ## if (nrow(overdue <- which((as.matrix(jab$adj)-adj.new)<0, arr.ind=T))>0) {
                ##     print("Edge copy deficit!")
                ##     browser()
                ## }

                ## intermediate check
                ## if (length(which(((adj.new + cn.adj) - jab$adj)!=0, arr.ind = TRUE)))
                ##     browser()

                to.rm = epaths[[i]][which(cn.adj[epaths[[i]]]==0), ,drop = FALSE]
                ## if (!palindromic) ## update reverse complement
                to.rm = rbind(to.rm, epaths[[i+1]][which(cn.adj[epaths[[i+1]]]==0), ,drop = FALSE])

                if (nrow(to.rm)>0)
                {
                    adj[to.rm] = 0
                    ## ALERT!!! major change
                    ## adjj = adj/as.matrix(cn.adj)
                    ## adjj[which(is.nan(adjj))] = 0
                    ## adjj[which(adjj<0)] = 0
                    G = graph.adjacency(adj, weighted = 'weight')
                    ## G = graph.adjacency(adjj, weighted = 'weight')
                    new.ends = setdiff(which(
                    (degree(G, mode = 'out')==0 | degree(G, mode = 'in')==0)
                    & degree(G)>0), ends)

                    ## ## check if cn.adj out of balance
                    ## if (any((Matrix::colSums(cn.adj)*Matrix::rowSums(cn.adj) != 0) & (Matrix::colSums(cn.adj) != Matrix::rowSums(cn.adj)))){
                    ##     print("Junction OUT OF BALANCE!")
                    ##     browser()
                    ## }

                    ## ## should be no new ends
                    ## if (length(new.ends)>0){
                    ##     print("Please, no new ends!")
                    ##     browser()
                    ## }

                    ## remain = as.matrix(jab$adj) - adj.new
                    ## nb <- which(Matrix::colSums(remain) != Matrix::rowSums(remain))
                    ## if (any(!is.element(nb, nb.all)))
                    ##     browser()

                    D = shortest.paths(G, v = ends, mode = 'out', weight = E(G)$weight)[, ends]
                    ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row != col, ][order(dist), ][, row := ends[row]][, col := ends[col]]
                }
                else
                    ij = ij[-1, , drop = FALSE]

                ## if (!palindromic) ## increase extra counter to account for reverse complement
                ## TOFIX: just update counter by 2 above, since we are just doing every path and its rc
                i = i+1
            }
        }


        ## DEBUG DEBUG DEBUG
        seg.ix = which(as.character(strand(jab$segstats))=='+'); seg.rix = which(as.character(strand(jab$segstats))=='-');


        if (nrow(ij)==0 & cleanup_mode == FALSE)
        {
            jmessage('!!!!!!!!!!!!!!!!!!!!!!!!!!STARTING CLEANUP MODE FOR PATHS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row != col, ][order(dist), ][, row := ends[row]][, col := ends[col]]
            cleanup_mode = TRUE
        }
    }
    if (verbose)
        jmessage('Path peeling iteration ', i, ' with ', sum(adj!=0, na.rm = TRUE), ' edges left ', nrow(ij) )

    ## ## record G, D, remaining edges at the end of path peeling
    ## G1 = G
    ## D1 = D
    ## remain1 = remain

    vpaths = vpaths[1:i]
    epaths = epaths[1:i]
    cns = cns[1:i]
    palindromic.path = palindromic.path[1:i]

    vcycles = rep(list(NA), maxrow)
    ecycles = rep(list(NA), maxrow)
    ccns = rep(NA, maxrow)

    csimp = which(diag(cn.adj)!=0)
    ipath = i
    i = 0
    if (length(csimp)>0)
    {
        vcycles[1:length(csimp)] = split(csimp, 1:length(csimp))
        ecycles[1:length(csimp)] = lapply(csimp, function(x) cbind(x, x))
        ccns[1:length(csimp)] = diag(cn.adj)[csimp]
        cn.adj[cbind(csimp, csimp)] = 0
        adj[cbind(csimp, csimp)] = 0
        i = length(csimp)

        for (j in 1:length(csimp))
            adj.new[ecycles[[j]]] = adj.new[ecycles[[j]]] + ccns[j]
    }

    ## sort shortest paths and find which connect a node to its ancestor (i.e. is a cycle)
    .parents = function(adj)
    {
        tmp = apply(adj, 2, function(x) which(x!=0))
        ix = which(sapply(tmp, length)>0)
        if (length(ix)>0)
        {
            parents = rbindlist(lapply(ix, function(x) data.table(x, tmp[[x]])))
            setnames(parents, c('node', 'parent'))
            setkey(parents, node)
        } else {
            parents = data.table(node = 0, parent = NA)
            setkey(parents, node)
        }
    }

    parents = .parents(adj)

    #' then find paths that begin at a node and end at (one of its) immediate upstream neighbors
    #' this will be a path for whom col index is = parent(row) for one of the rows
    ## ALERT!!! major change
    ## adjj = adj/as.matrix(cn.adj)
    ## adjj[which(is.nan(adjj))] = 0
    ## adjj[which(adjj<0)] = 0
    G = graph.adjacency(adj, weighted = 'weight')
    ## G = graph.adjacency(adjj, weighted = 'weight')
    D = shortest.paths(G, mode = 'out', weight = E(G)$weight)

    ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row %in% parents$parent & row != col, ][order(dist), ][, is.cycle := parents[list(row), col %in% parent], by = row][is.cycle == TRUE, ]


    ## now iterate from shortest to longest path
    ## peel that path off and see if it is still there ..
    ## and see if it is still there

    ## peel off top path and add to stack, then update cn.adj

    cleanup_mode = FALSE
    while (nrow(ij)>0)
    {
        if (verbose)
            jmessage('Cycle-peeling iteration ', i, ' with ', sum(adj!=0, na.rm = TRUE), ' edges left ', nrow(ij) )
        i = i+1
                                        #        p = as.numeric(get.shortest.paths(G, ij[1, 1], ij[1, 2], mode = 'out', weight = E(G)$weight)$vpath[[1]])

        p = get.constrained.shortest.path(cn.adj, G, allD = D, v = ij[1, 1], to = ij[1, 2], weight = E(G)$weight, edges = ed, verbose = TRUE, mip = cleanup_mode)

        if (is.null(p)){
            jmessage('Came up empty!')
            i = i -1
            ij = ij[-1, , drop = FALSE]
        } else
        {

            ed$cn = cn.adj[cbind(ed$from, ed$to)]
            vcycles[[i]] = p
            ecycles[[i]] = cbind(p, c(p[-1], p[1]))
            eids = paste(ecycles[[i]][,1], ecycles[[i]][,2])
            ccns[i] = ed[.(eids), if (length(cn)>1) cn/2 else cn, by = eclass][, floor(min(V1))] ## update cn correctly, adjusting constraints for palindromic edges by 1/2

            rvcycle = rtile.map[list(tile.map[list(vcycles[[i]]), -rev(tile.id)]), id]
            recycle = cbind(rvcycle, c(rvcycle[-1], rvcycle[1]))
            clen = length(rvcycle)
            hclen = floor(length(rvcycle)/2)
            ## (awkward) check for palindromicity for odd and even length palindromes

            ## if (all((vcycles[[i]]==rvcycle)[c(1:hclen,(clen-hclen+1):clen)]))
            if (ed[eids, any(table(eclass)>1)])
                palindromic.cycle[i] = TRUE
            ## else
            ## {
            vcycles[[i+1]] = rvcycle
            ecycles[[i+1]] = recycle
            ccns[i+1] = ccns[i]
            palindromic.cycle[i+1] = TRUE
            ##     palindromic = FALSE
            ## }
            ##        palindromic = TRUE ## set to true while we "figure things out"

            #' so now we want to subtract that cn units of that path from the graph
            #' so we want to update the current adjacency matrix to remove that path
            #' while keeping track of of the cycles on the stack
            cn.adj[ecycles[[i]]] = cn.adj[ecycles[[i]]]-ccns[i]
            ## if (!palindromic) ## update reverse complement unless palindromic
            cn.adj[ecycles[[i+1]]] = cn.adj[ecycles[[i+1]]]-ccns[i+1]

            if (!all(cn.adj[ecycles[[i]]]>=0))
            {
                jmessage('backtracking')
                ## browser()
                cn.adj[ecycles[[i]]] = cn.adj[ecycles[[i]]]+ccns[i]
                ## if (!palindromic) ## update reverse complement unless palindromic
                cn.adj[ecycles[[i+1]]] = cn.adj[ecycles[[i+1]]]+ccns[i+1]
                i = i-1
                ij = ij[-1, , drop = FALSE]
            }
            else
            {
                adj.new[ecycles[[i]]] = adj.new[ecycles[[i]]] + ccns[i]

                ## ## if (!palindromic)
                ##     adj.new[ecycles[[i+1]]] = adj.new[ecycles[[i+1]]] + ccns[i]

                ## ## ## make sure I didn't overuse any edge
                ## ## if (length(overdue <- which((as.matrix(jab$adj)-adj.new)<0))) {
                ## ##     print("Edge copy deficit!")
                ## ##     browser()
                ## ## }

                ## ## ## intermediate cross check
                ## ## if (length(which(((adj.new + cn.adj) - jab$adj)!=0, arr.ind = TRUE)))
                ## ##     browser()

                to.rm = ecycles[[i]][which(cn.adj[ecycles[[i]]]==0), ,drop = FALSE]

                ## if (!palindromic) ## update reverse complement
                to.rm = rbind(to.rm, ecycles[[i+1]][which(cn.adj[ecycles[[i+1]]]==0), ,drop = FALSE])

                if (nrow(to.rm)>0)
                {
                    adj[to.rm] = 0
                    parents = .parents(adj)
                    ## G = graph.adjacency(adj, weighted = 'weight')

                    ## ALERT!!! major change
                    ## adjj = adj/as.matrix(cn.adj)
                    ## adjj[which(is.nan(adjj))] = 0
                    ## adjj[which(adjj<0)] = 0
                    G = graph.adjacency(adj, weighted = 'weight')
                    ## G = graph.adjacency(adjj, weighted = 'weight')

                    ## if (any((Matrix::colSums(cn.adj)*Matrix::rowSums(cn.adj) != 0) & (Matrix::colSums(cn.adj) != Matrix::rowSums(cn.adj)))){
                    ##     print("Junction OUT OF BALANCE!")
                    ##     browser()
                    ## }

                    ## remain = as.matrix(jab$adj) - adj.new
                    ## nb <- which(Matrix::colSums(remain) != Matrix::rowSums(remain))
                    ## if (any(!is.element(nb, nb.all)))
                    ##     browser()

                    D = shortest.paths(G, mode = 'out', weight = E(G)$weight)
                    ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row %in% parents$parent & row != col, ][order(dist), ][, is.cycle := parents[list(row), col %in% parent], by = row][is.cycle == TRUE, ]
                }
                else
                    ij = ij[-1, ,drop = FALSE]

                ## if (!palindromic) ## increase extra counter to account for reverse complement
                i = i+1
            }
        }

        if (nrow(ij)==0 & cleanup_mode == FALSE)
        {
            jmessage('!!!!!!!!!!!!!!!!!!!!!!!!!!STARTING CLEANUP MODE FOR CYCLES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row %in% parents$parent & row != col, ][order(dist), ][, is.cycle := parents[list(row), col %in% parent], by = row][is.cycle == TRUE, ]

            cleanup_mode = TRUE
        }
    }

    if (verbose)
        jmessage('Cycle peeling iteration ', i, ' with ', sum(adj!=0, na.rm = TRUE), ' edges left ', nrow(ij) )


    if (i>0)
    {
        vcycles = vcycles[1:i]
        ecycles = ecycles[1:i]
        ccns = ccns[1:i]
    }
    else
    {
        vcycles = NULL
        ecycles = NULL
        ccns = NULL
    }

    vall = c(vpaths, vcycles)
    eall = c(epaths, ecycles)
    ecn = c(cns, ccns)

    ## ## record G, D, remaining edges at the end of cycle peeling
    ## G2 = G
    ## D2 = D
    ## remain2 = remain
    remain = as.matrix(jab$adj) - adj.new
    remain.ends = which(Matrix::colSums(remain)*Matrix::rowSums(remain)==0 & Matrix::colSums(remain)-Matrix::rowSums(remain)!=0)
    if (length(remain.ends)>0){
        if (verbose)
            jmessage(length(remain.ends), "ends were not properly assigned a path. Do them.")
    }

    tmp = cbind(do.call(rbind, eall), rep(ecn, sapply(eall, nrow)), munlist(eall))
    ix = which(Matrix::rowSums(is.na(tmp[, 1:2]))==0)

    if (length(ix)>0)
        adj.new = sparseMatrix(tmp[ix,1], tmp[ix,2], x = tmp[ix,3], dims = dim(adj))
    else
        adj.new = sparseMatrix(1, 1, x = 0, dims = dim(adj))
    vix = munlist(vall)

    jab$segstats$node.id = 1:length(jab$segstats)
    pathsegs = jab$segstats[vix[,3]]
    pathsegs$grl.ix = vix[,1]
    abjuncs =  as.data.table(rbind(jab$ab.edges[, 1:2, '+'], jab$ab.edges[, 1:2, '-']))[,
                                   id := rep(1:nrow(jab$ab.edges),2)*
                                       rep(c(1, -1), each = nrow(jab$ab.edges))][!is.na(from), ]
    abjuncs = abjuncs[, tag := structure(paste(from, to), names = id)]
    setkey(abjuncs, tag)

    ## annotate ab.id (if any) following each segment in each path
    pathsegs$ab.id = gr2dt(pathsegs)[ , ab.id := c(abjuncs[paste(node.id[-length(node.id)], node.id[-1]), id], NA), by = grl.ix][, ab.id]

    paths = split(pathsegs, vix[,1])
    values(paths)$ogid = 1:length(paths)
    values(paths)$cn = ecn[as.numeric(names(paths))]
    values(paths)$label = paste('CN=', ecn[as.numeric(names(paths))], sep = '')
    values(paths)$is.cycle = !(as.numeric(names(paths)) %in% 1:length(vpaths))
    values(paths)$numsegs = elementNROWS(paths)
    values(paths)$num.ab = sapply(paths, function(x) sum(!is.na(x$ab.id)))
    values(paths)$wid = sapply(lapply(paths, width), sum)

    check = which((adj.new - jab$adj) !=0, arr.ind = TRUE)

    if (length(check)>0)
        stop('Alleles do not add up to marginal copy number profile!')
    else if (verbose)
        jmessage('Cross check successful: sum of walk copy numbers = marginal JaBbA edge set!')

    ## match up paths and their reverse complements
    psig = lapply(paths, function(x) ifelse(as.logical(strand(x)=='+'), 1, -1)*x$tile.id)
    psig.flip = sapply(psig, function(x) -rev(x))

    unmix = data.table(
        ix = 1:length(paths),
        mix = match(sapply(psig, paste, collapse = ','), sapply(psig.flip, paste, collapse = ',')))[, pos := 1:length(mix)<mix][order(!pos), ]
    setkey(unmix, ix)
    unmix[is.na(mix), pos := TRUE] ## if we have paths with no reverse complement i.e. NA mix, then call "+" for now

    remix = rbind(
        unmix[pos == TRUE, ][, id := 1:length(ix)],
        unmix[list(unmix[pos == TRUE, mix]), ][, id := 1:length(ix)][!is.na(ix), ]
    )

    paths = paths[remix$ix]
    names(paths) = paste(remix$id, ifelse(remix$pos, '+', '-'), sep = '')
    values(paths)$id = remix$id
    values(paths)$str = ifelse(remix$pos, '+', '-')

    if (length(setdiff(values(paths)$ogid, 1:length(paths))))
        jmessage('Warning!!! Some paths missing!')

    ## for gGnome compatibiliity
    if (!return.grl)
    {
      tmp.dt = as.data.table(copy(paths))[, pid := group_name][, nix := 1:.N, by =pid]
      setkeyv(tmp.dt, c('pid', 'nix'))

      ## mark nodes that precede a reference junction
      tmp.dt[, d.to.next := c((start-data.table::shift(end))[-1], NA), by = pid]
      tmp.dt[, d.to.next.neg := c((data.table::shift(start)-end)[-1], NA), by = pid]
      tmp.dt[, same.strand := c((strand==data.table::shift(strand))[-1], NA), by = pid]
      tmp.dt[, same.chrom := c((as.character(seqnames)==data.table::shift(as.character(seqnames)))[-1], NA), by = pid]
      tmp.dt[, last.node := 1:.N == .N, by = pid]
      tmp.dt[, before.ref :=
                 (((d.to.next<=1 & d.to.next>=0 & strand == '+') |
                   (d.to.next.neg<=1 & d.to.next.neg>=0 & strand == '-')
                 ) & same.strand & same.chrom)]
      tmp.dt[is.na(before.ref), before.ref := FALSE]

      ## label reference runs of nodes then collapse
      .labrun = function(x) ifelse(x, cumsum(diff(as.numeric(c(FALSE, x)))>0), as.integer(NA))
      tmp.dt[, ref.run := .labrun(before.ref), by = pid]
      tmp.dt[, ref.run.last := data.table::shift(ref.run), by = pid]
      tmp.dt[is.na(ref.run) & !is.na(ref.run.last), ref.run := ref.run.last]
      tmp.dt[!is.na(ref.run), ref.run.id := paste(pid, ref.run)]
      #tmp.dt[loose == TRUE, ref.run.id := NA] ## make sure loose ends stay ungrouped
      collapsed.dt = tmp.dt[!is.na(ref.run.id), .(
                                                  nix = nix[1],
                                                  pid = pid[1],
                                                  seqnames = seqnames[1],
                                                  start = min(start),
                                                  end = max(end),
                                                  loose = FALSE,
                                                  strand = strand[1]
                                                ), by = ref.run.id]

      ## concatenate back with nodes that precede a non reference junction
      tmp.dt = rbind(tmp.dt[is.na(ref.run.id), .(pid, nix, seqnames, start, end, strand, loose)],
                     collapsed.dt[, .(pid, nix, seqnames, start, end, strand, loose)])
      setkeyv(tmp.dt, c('pid', 'nix'))

      tmp.gr = dt2gr(tmp.dt)
      tmp.segs = unique(tmp.gr)
      tmp.gr$seg.id = match(tmp.gr, tmp.segs)
      tmp.paths = split(tmp.gr$seg.id, tmp.gr$pid)
      tmp.vals = as.data.frame(values(paths[names(tmp.paths)]))

      names(tmp.paths) = ifelse(grepl('\\-', names(tmp.paths)), -1, 1)*as.numeric(gsub('\\D', '', names(tmp.paths)))

      ##      gw = gGnome::gWalks$new(segs=tmp.segs,
      gw = gWalks$new(segs=tmp.segs,
                      paths=tmp.paths,
                      metacols=tmp.vals)
      return(gw)
    }


    return(paths)
}



get.constrained.shortest.path = function(cn.adj, ## copy number matrix
                                         G, ## graph with distances as weights
                                         allD=NULL, ## shortest path between all nodes in graph
                                         v,
                                         to,
                                         weight,
                                         edges,
                                         verbose = TRUE,
                                         mip = TRUE
                                         )
{

    if (is.null(allD)) allD = shortest.paths(G, mode="out", weights = weight)

    v = as.numeric(v)
    to = as.numeric(to)

    if (is.infinite(allD[v, to]) | allD[v, to]==0) return(NULL)

    edges$cn = cn.adj[cbind(edges$from, edges$to)]

  ## ASSUME: from, to are scalars, within node range, to is reachable from from
  ## ASSUME edges contains eid key and eclass mapping

  ## getting shortest path and associated nodes and edges
    tmp.p = as.numeric(get.shortest.paths(G, from=v, to=to, "out", weights=weight)$vpath[[1]])
    tmp.e = cbind(tmp.p[-length(tmp.p)], tmp.p[-1])
    tmp.eid = paste(tmp.e[, 1], tmp.e[, 2])
    tmp.eclass = edges[.(tmp.eid), eclass] ## grouping edges by reverse complement

    ## the cn of this path is the max number of copies that the network will allow
    ## here we have to group by eclass, i.e. so if there are two edges from an eclass
    ## in a given path then we need to halve the "remaining copies" constraint
    tmp.pcn = edges[.(tmp.eid), if (length(cn)>1) cn/2 else cn, by = eclass][, floor(min(V1))]

    edges[, rationed := cn<(tmp.pcn*2)]

    D.totarget = allD[, as.numeric(to)]
    edges[, distance_to_target :=  D.totarget[to]]
    edges = edges[!is.infinite(distance_to_target) & cn>0, ]

    rationed.edges = edges[rationed == TRUE, ]

    ## find overdrafted eclasses - meaning two instances in this path but only one remaining copy
    overdrafts.eclass = intersect(names(which(table(tmp.eclass)==2)), rationed.edges$eclass)

    first.overdraft = which(tmp.eclass %in% overdrafts.eclass & duplicated(tmp.eclass))[1]

    ## no overdrafts?, then return
    if (is.na(first.overdraft) & tmp.pcn>0)
    {
        if (verbose)
            jmessage('Shortest path is good enough!')
        return(tmp.p)
    }

    if (!mip)
        return(NULL)

    ## use MIP to find constrained path
    edges[, enum := 1:length(eid)]

    ## incidence matrix constraints + 1 for tmp.pcn
    A = sparseMatrix(edges$to, edges$enum, x = 1, dims = c(nrow(cn.adj), nrow(edges))) -
        sparseMatrix(edges$from, edges$enum, x = 1, dims = c(nrow(cn.adj), nrow(edges)))
    b = rep(0, nrow(A))
    b[v] = -1
    b[to] = 1

    ix = which(Matrix::rowSums(A!=0)!=0) ## remove zero constraints

    ## "ration" or reverse complementarity constraints
    tmp.constraints = edges[, list(e1 = enum[1], e2 = enum[2], ub = cn[1]), by = eclass]
    tmp.constraints = tmp.constraints[!is.na(e1) & !is.na(e2), ]

    R = sparseMatrix(rep(1:nrow(tmp.constraints), 2),
                     c(tmp.constraints$e1, tmp.constraints$e2),
                     x = 1, dims = c(nrow(tmp.constraints), nrow(edges)))
    Rb = tmp.constraints$ub

    ## minimize weight of path
    c = edges$weight

    res = Rcplex(c, rbind(A[ix,], R), c(b[ix], Rb), sense = c(rep('E', length(ix)), rep('L', length(Rb))),
                 lb = 0, vtype = "B",
                 objsense = 'min')


    if (verbose)
        jmessage('YES WE ARE DOING PROPER MIP!!!!')

    if (!(res$status %in% c(101, 102)))
    {
        if (verbose)
            jmessage('No solution to MIP!')

#        browser()
        return(NULL)
    }

    ## use igraph to sort these edges into a path, i.e. make simple graph with one path and extract it using igraph (lazy :)
    tmp.p = as.numeric(get.shortest.paths(graph_from_edgelist(edges[res$xopt!=0, cbind(from, to)]), v, to)$vpath[[1]])

    ## check if overdrafted
    if (verbose)
    {
        tmp.e = cbind(tmp.p[-length(tmp.p)], tmp.p[-1])
        tmp.eid = paste(tmp.e[, 1], tmp.e[, 2])
        tmp.eclass = edges[.(tmp.eid), eclass]
        tmp.pcn = edges[.(tmp.eid), if (length(cn)>1) cn/2 else cn, by = eclass][, min(V1)]
        overdrafts.eclass = intersect(names(which(table(tmp.eclass)==2)), rationed.edges$eclass)
        if (length(overdrafts.eclass)==0)
            jmessage('No overdrafts after MIP')
        else
        {
            jmessage('Still overdraft!')
            browser()
        }
    }

    return(tmp.p)
}



####################################################
#' @name jabba.gwalk
#' @export
#' @rdname internal
#' jabba.walk
#'
#' Computes walks around all aberrant edges in JABbA object
#'
#' Takes in JaBbA solution and computes local
#' reconstructions around all aberrant edges (default).  Reconstructions (i.e. Huts) consists
#' of collections of walks, each walk associated with a copy number, and a given
#' region (collection of genomic windows).  The interval sum of walks in a given region, weighted
#' by copy numbers will recapitulate the marginal copy profile (as estimated by JaBbA).
#' The reconstruction is chosen to maximize parsimony.
#'
#' Optional flags allow making huts around specific junctions or specified loci (GRangesList)
#'
#' Walks are reconstructed locally within "clustersize" nodes of each aberrant edge, where
#' clustersize is measured by the number of total edges.  Larger cluster sizes may fail to be
#' computationally tractable, i.e. with a highly rearranged genome in an area of dense interconnectivity.
#'
#' @param sol JaBbA object
#' @param outdir output directory
#' @param junction.ix junction indices around which to build walks (default is all junctions)
#' @param loci  loci around which to build walks (over-rides junction.ix), alternatively can be a list of  "all.paths" objects (i.e. each a list utput of initial all.paths = TRUE run  +/- field $prior for walk to re-eval a given all.paths combo
#' @param clustersize size of the cluster to output around the locus or junction of interest
#' @param trim logical flag whether trim in neighborhood of junction (only applicable if loci = NULL, default = TRUE)
#' @param trim.w integer width to which to trim to
#' @param prune flag whether to prune trivial walks for whom a path can be drawn from first to last interval in a graph linking intervals with pairwise distance < d1 on the walk or distance < d2 on the reference
#' @param prune.d1 local distance threshold for walk pruning
#' @param prune.d2 referenc distance threshold for walk pruning
#' @param mc.cores number of cores to use, default 1
#' @param genes character vector of gene symbols with which to annotate walk (eg cancer genes)
#' @param verbose logical flag
#' @return list of walk set around each locus or junction that is inputted to analysis, each list item is a list with the following fields
#' $win = input locus of interest, $grl = GRangesList of walks, $grs is a collapsed footprint of all walks in the walk list for this locu
#' $gtrack gTrack of of the output, additional outputs for debugging: $sol, $K, $Bc, $eix, $vix, $h
####################################################
jabba.walk = function(sol, kag = NULL, digested = T, outdir = 'temp.walk', junction.ix = NULL, loci = NULL, clustersize = 100,
  trim = FALSE, ## whether to trim around junction (only applicable when loci = NULL)
  trim.w = 1e6, ## how far to trim in neighborhood of junction (only applicable when loci = NULL
  prune = FALSE, ## whether to prune trivial walks i.e. those for whom a path can be drawn from first to last interval in a graph linking intervals with pairwise distance < d1 on the walk or distance < d2 on the reference
  prune.d1 = 1e5, ## local distance threshold for walk pruning
  prune.d2 = 1e5, ## reference distance threshold for walk pruning
  maxiterations = Inf, mc.cores = 1, genes = read.delim('~/DB/COSMIC/cancer_gene_census.tsv', strings = F)$Symbol, verbose = T, max.threads = 4, customparams = T, mem = 6, all.paths = FALSE, nomip = F, tilim = 100, nsolutions = 100, cb.interval = 1e4, cb.chunksize = 1e4, cb.maxchunks = 1e10)
{
  system(paste('mkdir -p', outdir))
  ## awkward workaround to limit the number of processors Cplex will gobble up
  ##

  if (customparams)
    {
      out.file = paste(outdir, 'tmp.prm', sep = '/')
      max.threads = Sys.getenv("LSB_DJOB_NUMPROC")
      if (nchar(max.threads) == 0)
        max.threads = Inf
      else
        max.threads = as.numeric(max.threads)
      max.threads = min(max.threads, mc.cores)
      if (is.infinite(max.threads))
        max.threads = 0

      param.file = paste(out.file, '.prm', sep = '')
      .cplex_customparams(param.file, max.threads, treememlim = mem * 1e3)

      Sys.setenv(ILOG_CPLEX_PARAMETER_FILE = normalizePath(param.file))
      print(Sys.getenv('ILOG_CPLEX_PARAMETER_FILE'))
    }


   if (is.null(sol))
      sol = kag

  if (is.null(sol$segstats))
      {
          sol$segstats = sol$tile
          sol$segstats$cn = 2
          sol$segstats$eslack.out = 0
          sol$segstats$eslack.in = 0
      }

  if (is.null(kag))
      kag = sol


  out = list()
  tmp.adj = sol$adj
  if (digested)  ## if input is already "digested", then don't need to bother with slacks
      {
      sol$segstats$eslack.in = 0
      sol$segstats$eslack.out = 0
      G = sol$G
    }
  else ## soon to be deprecated
    {
      ix = which(sol$segstats$eslack.in!=0 | sol$segstats$eslack.out!=0)
      tmp.adj[ix, ix] = 0
      pos.ix = which(as.character(strand(sol$segstats))=='+')
      sol$segstats$tile.id = match(gr.stripstrand(sol$segstats), gr.stripstrand(sol$segstats[pos.ix]))
      G = graph.adjacency(tmp.adj!=0)
    }

  if (verbose)
    jmessage(paste('Processing JaBbA'))

  h = jbaMIP.process(sol)

  if (!is.null(genes))
    td.rg = track.gencode(genes = genes, height = 3)

  if (is.null(junction.ix) & is.null(loci))
    junction.ix = 1:nrow(kag$ab.edges)

  if (!is.null(junction.ix))
    if (is.null(names(junction.ix)))
      names(junction.ix) = 1:length(junction.ix)

  if (is.null(loci)) ## junction.ix should be not null here
    {
      loci = do.call('GRangesList', mclapply(junction.ix, function(i)
        {
             if (verbose)
               cat(paste('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\nDefining subgraph around junction', i, '\n'))
             vix = vix.i = setdiff(kag$ab.edges[i, 1:2, ], NA)
             if (length(vix)==0)
                 return(GRanges())
             k = 0
             last.clustersize = 0
             while (length(vix)<clustersize & k < maxiterations & length(vix)>last.clustersize)
               {
                 k = k + 1
                 last.clustersize = length(vix)
                 vix = unique(unlist(neighborhood(G, vix.i, order = k)))
               }
             if (verbose)
               cat(paste('Outputting', length(vix), 'vertices around junction', i, '\n'))

             return(kag$segstats[vix])
           }
        , mc.cores = mc.cores))

      names(loci) = names(junction.ix)
      loci = loci[sapply(loci, length)>0]
    }
  else ## if loci are provided (i.e. not junction centric) then we will not trim or prune
    {
      trim = F
      prune = F
    }

  if (verbose)
    cat(paste('Finished defining subgraphs\n'))

  starts = gr.start(sol$segstats, ignore.strand = F)
  ends = gr.end(sol$segstats, ignore.strand = F)

  names(sol$segstats) = 1:length(sol$segstats)

  if (is.null(names(loci)))
    lnames =  paste('locus', 1:length(loci), sep = '')
  else
    lnames = names(loci)

  all.junc.pair = c(paste(sol$ab.edges[, 1, 1], sol$ab.edges[, 2, 1], sep = ','), paste(sol$ab.edges[, 1, 2], sol$ab.edges[, 2, 2], sep = ','))
  names(all.junc.pair) = c(1:nrow(sol$ab.edges), -c(1:nrow(sol$ab.edges)))

  if (length(loci)>0)
    {
      out = mclapply(1:length(loci), function(i)
          {
            label = lnames[i]
            mprior = NULL
              outfile.rds = sprintf('%s/%s.rds', outdir, label)
              outfile.pdf = sprintf('%s/%s.pdf', outdir, label)
              outfile.txt = sprintf('%s/%s.txt', outdir, label)
              outfile.allpaths.txt = sprintf('%s/%s.allpaths.txt', outdir, label)
              if (is(loci[[i]], 'GRanges'))
                  {
                      vix = which(gr.in(kag$segstats, loci[[i]]))
                      cat('Number of vertices:', length(vix), '\n')
                      eix = which((h$e.ij[,1] %in% vix | h$e.ij[,2] %in% vix) & h$e>0)
                      Bc = as.matrix(h$B)[vix, eix]
                      K = tryCatch(convex.basis(Bc, interval = cb.interval, chunksize = cb.chunksize, verbose = T, maxchunks = cb.maxchunks), error = function(e) as.character(e))
                      if (is.character(K))
                          return(list(README = K))
                      prior = rep(1, ncol(K))
                  }
              else ## assume we are re-heating a previous all.paths = TRUE output (and presumably adding a prior)
                  {
                      K = loci[[i]]$K
                      h = loci[[i]]$h
                      eix = loci[[i]]$eix
                      Bc = loci[[i]]$Bc
                      vix = loci[[i]]$vix
                      prior = rep(1, ncol(K))

                      if (is.matrix(loci[[i]]$mprior))
                      {
                        if (verbose)
                          cat(paste('Adding a matrix prior!!!!!!\n'))

                        ## initialize matrix
                        mprior = array(0, dim = c(nrow(loci[[i]]$mprior), ncol(K)))
                        loci
                        mprior[, values(loci[[i]]$allpaths.og)$kix] = loci[[i]]$mprior
                        mprior[, values(loci[[i]]$allpaths.og)$kix2] = loci[[i]]$mprior
                        colnames(mprior) = 1:ncol(mprior)
                        colnames(mprior)[values(loci[[i]]$allpaths.og)$kix] = as.character(1:ncol(mprior))
                        colnames(mprior)[values(loci[[i]]$allpaths.og)$kix2] = as.character(-(1:ncol(mprior)))
                      }

                      if (!is.null(loci[[i]]$prior))
                      {
                        if (verbose)
                          cat(paste('Adding a prior!!!!!!\n'))
                        prior[c(values(loci[[i]]$allpaths.og)$kix,values(loci[[i]]$allpaths.og)$kix2)]  = loci[[i]]$prior
                      }
                      loci[[i]] = loci[[i]]$win
                  }

              is.cyc = Matrix::colSums(K[h$etype[eix] == 'slack', ])==0 & Matrix::colSums((h$B[, eix, drop = F] %*% K)!=0)==0

              karyo.sol = karyoMIP(K, h$e[eix], h$eclass[eix], nsolutions = nsolutions, tilim = tilim, cpenalty = 1/prior, mprior = mprior)
              kag.sol = karyo.sol[[1]]
          p = karyoMIP.to.path(kag.sol, K, h$e.ij[eix, ], sol$segstats, mc.cores = pmin(4, mc.cores))
          values(p$grl)$cn = p$cn
          values(p$grl)$is.cyc = p$is.cyc
          td.rg$stack.gap = 5e6

          if (!is.null(kag$junctions))
            {
              values(kag$junctions)$lwd = sol$adj[kag$ab.edges[,1:2, 1]]
              values(kag$junctions)$lty = 1
              values(kag$junctions)$label = ifelse(sol$adj[kag$ab.edges[,1:2, 1]]>0, sol$adj[kag$ab.edges[,1:2, 1]], '')
              values(kag$junctions)$col = ifelse(sol$adj[kag$ab.edges[,1:2, 1]]>0, alpha('red', 0.3), alpha('white', 0))
            }
          win = streduce(sol$segstats[vix], 1e4)

          y1 = max(sol$segstats$cn[gr.in(sol$segstats, win)], na.rm = T)*1.1
          pdf(outfile.pdf, height = 30, width = 24)
          grs = gr.simplify(grl.unlist(p$grl), 'grl.ix', split = T)
          values(grs) = values(p$grl)
          names(grs) = names(p$grl)

          if (!is.null(sol$gtrack))
            {
              td.seg = sol$gtrack
              td.seg$y1 = y1
              td = c(td.seg, td.rg)
            }
          else
              {
                  td.seg = gTrack(sol$segstats, y.field = 'cn', angle = 0, col ='black', height = 6, labels.suppress = T, y1 = y1)

                                        #          td = c(gTrack(grs, draw.paths = T, path.cex.arrow = 0, border = NA, angle = 0, ywid = 0.5, path.stack.x.gap = 1e6, height = 20, labels.suppress.gr = T),

                  gt.walk = gTrack(grs, draw.paths = T, border = NA, angle = 0, ywid = 0.5, height = 20, labels.suppress.gr = T)
                  gt.walk$path.cex.arrow = 0
                  gt.walk$path.stack.x.gap = 1e6
                  td = c(
                      gt.walk,
                      td.seg,
                      td.rg)
                  gTrack::plot(td,
                       windows = win, links = kag$junctions)
                  dev.off()
              }

          df = data.frame(label = label, cn = p$cn, walk = sapply(grs, function(x) paste(gr.string(x, mb = F), collapse = ',')), widths = sapply(grs, function(x) paste(width(x), collapse = ',')), width = sapply(grs, function(x) sum(width(x))), numpieces = sapply(grs, length), type = 'walk')
          df = rbind(data.frame(label = label, cn = NA, walk = paste(gr.string(win, mb = F), collapse = ','), widths = paste(width(win), collapse = ','), width = sum(width(win)), type = 'window', numpieces = length(win)), df)
          write.tab(df, outfile.txt)
          out = list(
            win = win, grl = p$grl, grls = grs, td = td, sol = karyo.sol,
            K = K, Bc = Bc, eix = eix, vix = vix, h = h,
            README = 'win=windows, grl = raw granges list corresponding to paths, grls = simplified granges list corresponding to paths, td = gTrack object plotting walks, sol = solution object from karyoMIP of local walks, K = incidence matrix input to karyomip, Bc = input to convex.basis, eix = eix input to karyomip, vix = vix input corresponding to rows of Bc, h = h input to karyomip')

          if (all.paths)
            {
              outfile.allpaths.pdf = sprintf('%s/%s.allpaths.pdf', outdir, label)

              if (verbose)
                cat('Generating all walks\n')

              ## repurpose karyoMIP.to.path to generate all paths using "fake solution" i.e. all 1 weights,  to karyoMIP as input
              pallp = karyoMIP.to.path(list(kcn = kag.sol$kcn*0 + 1, kclass = kag.sol$kclass), K, h$e.ij[eix, ], sol$segstats, mc.cores = pmin(4, mc.cores), verbose = verbose)
              allp = pallp$grl

              allps = gr.simplify(grl.unlist(allp), 'grl.ix', split = T)
              allps[values(allp)$is.cycle] = do.call('GRangesList', lapply(which(values(allp)$is.cycle), function(x) c(allps[[x]], allps[[x]])))
              allps.og = allps; ## save for later
              values(allps.og)$kix = pallp$kix
              values(allps.og)$kix2 = pallp$kix2

              ## text encoding of junctions
              if (!is.null(junction.ix))
                junc.pair = paste(sol$ab.edges[junction.ix[i], 1, ], sol$ab.edges[junction.ix[i], 2, ], sep = ',')

              if (trim | prune) ## junction.ix should be not null here (i.e. they were provided as input or loci = NULL)
                {
                  allps.u = grl.unlist(allps)
                  allps.u$ix.s = gr.match(gr.start(allps.u, ignore.strand = F), starts, ignore.strand = F)
                  allps.u$ix.e = gr.match(gr.end(allps.u, ignore.strand = F), ends, ignore.strand = F)
                  allps = split(allps.u, allps.u$grl.ix)
                  allps.ixs = split(allps.u$ix.s, allps.u$grl.ix) ## start indices of walk intervals in sol$segstats
                  allps.ixe = split(allps.u$ix.e, allps.u$grl.ix) ## end indices of walks intervals in sol$segstats
                  allps.w = split(width(allps.u), allps.u$grl.ix)
                  allps.endc = split(levapply(width(allps.u), by = list(allps.u$grl.ix), FUN = cumsum), allps.u$grl.ix)

                  if (trim) ## only include windows around the junction of interest
                    {
                      ## allps.ix.pairs tells us what junction indices are present in a walk collection
                      allps.ix.pairs = mapply(function(x,y) if (length(x)<=1) NULL else which(paste(x[-length(x)], y[-1], sep = ',') %in% junc.pair), allps.ixe, allps.ixs, SIMPLIFY = F)
                      ## first, which windows contain the junction

                      wix = which(sapply(allps.ix.pairs, length)>0)
                      allps = allps[wix]

                      if (length(allps)>0)
                        {
                          allps.ixs = allps.ixs[wix] ## start interval id of kth interval in ith walk
                          allps.ixe = allps.ixe[wix] ## end interval id of kth interval in ith walk
                          allps.endc = allps.endc[wix] ## end walk coordinate of kth interval in ith walk
                          allps.w = allps.w[wix]
                          allps.ix.pairs = allps.ix.pairs[wix]

                          ## start window for trimming
                          values(allps)$allps.junc.first =
                            pmax(0, mapply(function(x, y) y[x[1]], allps.ix.pairs, allps.endc)) ## walk position of first junction
                          values(allps)$allps.junc.last =
                            pmax(0, mapply(function(x, y) y[x[length(x)]], allps.ix.pairs, allps.endc)) ## walk position of last junction

                          ## check for any quasi-palindromic walks that contain both orientations of a junction
                          ## split each of these into two so we can maintain the width limit
                          pal.wix = which(values(allps)$allps.win.firstix != values(allps)$allps.win.lastix)
                          if (length(pal.wix)>0)C
                            {
                              allps.dup = allps[pal.wix]
                              values(allps.dup)$allps.junc.first = values(allps)$allps.junc.last
                              allps = c(allps, allps.dup)
                              allps.endc = c(allps.endc, allps.endc[pal.wix])
                              allps.w = c(allps.w, allps.w[pal.wix])
                            }

                          values(allps)$allps.win.first =
                            pmax(0, values(allps)$allps.junc.first - trim.w) ## walk coordinate of new window start
                          values(allps)$allps.win.last =
                            pmin(sapply(allps.endc, function(x) x[length(x)]), values(allps)$allps.junc.first + trim.w) ## walk coordinate of new window end
                          values(allps)$allps.win.firstix = ## first walk interval to trim to
                            mapply(function(x, y) setdiff(c(which(x>y)[1], 1), NA)[1], allps.endc, values(allps)$allps.win.first)
                          values(allps)$allps.win.lastix = ## last walk interval to trim to
                            mapply(function(x, y) setdiff(c(which(x>y)[1], length(x)), NA)[1], allps.endc, values(allps)$allps.win.last)
                          values(allps)$allps.win.first.keep =
                            mapply(function(p,e,i) e[i] - p, values(allps)$allps.win.first, allps.endc, values(allps)$allps.win.firstix)
                          values(allps)$allps.win.last.keep =
                            mapply(function(p,e,i,w) w[i] - (e[i] - p), values(allps)$allps.win.last, allps.endc, values(allps)$allps.win.lastix, allps.w)
                          ## apply trimming
                          ## we are trimming walks so that they are within trim.w bases of junction
                          allps.u = grl.unlist(allps)
                          iix = mapply(function(x,y) y %in% values(allps)$allps.win.firstix[x]:values(allps)$allps.win.lastix[x], allps.u$grl.ix, allps.u$grl.iix)
                          allps.u = allps.u[iix]
                          allps.u$keep.end = mapply(function(x, y)
                            ifelse(y == values(allps)$allps.win.firstix[x], values(allps)$allps.win.first.keep[x], NA), allps.u$grl.ix, allps.u$grl.iix)
                          allps.u$keep.start = mapply(function(x, y)
                            ifelse(y == values(allps)$allps.win.lastix[x], values(allps)$allps.win.last.keep[x], NA), allps.u$grl.ix, allps.u$grl.iix)

                          if (any(tmp.ix <- !is.na(allps.u$keep.start))) ## we keep the end of the first segment
                            allps.u[tmp.ix] = gr.start(allps.u[tmp.ix], allps.u$keep.start[tmp.ix], ignore.strand = F)

                          if (any(tmp.ix <- !is.na(allps.u$keep.end))) ## we keep the beginning of the last segment
                            allps.u[tmp.ix] = gr.end(allps.u[tmp.ix], allps.u$keep.end[tmp.ix], ignore.strand = F)

                          ## if there are multiple walks with the same aberrant junction set, then pick the longest of these

                          ## first need to find the aberrant walks in each set
                          ij = paste(allps.u$ix.e[-length(allps.u)], allps.u$ix.s[-1], sep = ',') ## indices of all walk adjacent interval pairs
                          names(ij) = 1:length(ij)
                          ij = ij[diff(allps.u$grl.ix)==0] ## only pick intra-walk interval pairs
                          ij.ix = names(all.junc.pair)[match(ij, all.junc.pair)]
                          ## then compute the width of each walk

                          allps = split(allps.u, allps.u$grl.ix)
                          ij.ix.l = split(ij.ix, allps.u$grl.ix[as.numeric(names(ij))])[names(allps)]
                          values(allps)$ab.junc = lapply(ij.ix.l, paste, collapse = ',')
                          values(allps)$wid = vaggregate(width(allps.u), by = list(allps.u$grl.ix), FUN = sum)[names(allps)]
                          ix.w = order(-values(allps)$wid)
                          allps = allps[ix.w[which(!duplicated(values(allps)$ab.junc[ix.w]))]] ## only keep the longest non-duplicate walks
                        }
                    }

                  ## now dedup and trim contigs to locus (mainly useful if loci was provided as argument)
                  if (length(allps)>0)
                    {
                      win = reduce(gr.stripstrand(loci[[i]]))
                      allps.u = grl.unlist(allps)

                      ## trim to locus
                      ix = gr.match(allps.u, win)
                      allps.u = allps.u[!is.na(ix)]
                      ix = ix[!is.na(ix)]
                      start(allps.u) = pmax(start(allps.u), start(win)[ix])
                      end(allps.u) = pmin(end(allps.u), end(win)[ix])

                      allps.u$ix.s = gr.match(gr.start(allps.u, ignore.strand = F), starts, ignore.strand = F)
                      allps.u$ix.e = gr.match(gr.end(allps.u, ignore.strand = F), ends, ignore.strand = F)

                      ## remove dups
                      allps.ixs = split(allps.u$ix.s, allps.u$grl.ix) ## start indices of intervals
                      allps.ixe = split(allps.u$ix.e, allps.u$grl.ix) ## end indices of intervals

                      allps.u = allps.u[allps.u$grl.ix %in% which(!duplicated(paste(sapply(allps.ixs, paste, collapse = ','), sapply(allps.ixe, paste, collapse = ','))))]
                      allps = split(allps.u, allps.u$grl.ix)
                    }


                  if (prune & length(allps)>0)
                    ## this is to prune pseudo-aberrant walks that basically consist of short insertions of non-reference
                    ## sequences in a big reference chunk
                    {
                      ## for each walk create graph of intervals by determining whether pair ij is BOTH near on the walk (<= d1)
                      ## and near on the refernce (<= d2)
                      allps.u = grl.unlist(allps)

                      ## what are the ij pairs we want to test from this collapsed list
                      ij = merge(cbind(i = 1:length(allps.u), ix = allps.u$grl.ix), cbind(j = 1:length(allps.u), ix = allps.u$grl.ix))[, c('i', 'j')]

                      tmp = levapply(width(allps.u), by = list(allps.u$grl.ix), FUN = cumsum)
                      allps.u.ir = IRanges(tmp - width(allps.u) + 1, tmp)

                      ## distance on the walk
                      D1 = sparseMatrix(ij[, 'i'],ij[, 'j'],
                        x = suppressWarnings(
                          distance(IRanges(start = end(allps.u.ir[ij[,'i']]), width = 1),
                                   IRanges(start(allps.u.ir[ij[,'j']]), width = 1))) + 1e-5, dims = rep(length(allps.u.ir), 2))

                      ## distance on the reference
                      D2 = sparseMatrix(ij[, 'i'],ij[, 'j'],
                        x = suppressWarnings(
                          distance(gr.end(allps.u[ij[,'i']], ignore.strand = F),
                                   gr.start(allps.u[ij[,'j']], ignore.strand = F))) + 1e-5, dims = rep(length(allps.u.ir), 2))

                      D1 = pmin(as.matrix(D1), as.matrix(t(D1)))
                      D2 = pmin(as.matrix(D2), as.matrix(t(D2)))

                      tmp = D1>0 & D1<prune.d1 & D2>0 & D2<prune.d2
                      tmp[which(is.na(tmp))] = FALSE
                      G = graph.adjacency(tmp)
                      cl = clusters(G, 'weak')$membership ## clusters based on this adjacency relationship
                      cls = split(1:length(cl), cl)
                      lens = sapply(allps, length)

                      ## check if there any clusters that contain both the first and last member  of a walk
                      cls.fl = cls[mapply(function(x) all(c(1,lens[allps.u$grl.ix[x[1]]]) %in% allps.u$grl.iix[x]), cls)]

                      if (length(cls.fl)>0)
                        {
                          toprune = allps.u$grl.ix[sapply(cls.fl, function(x) x[1])]
                          if (length(toprune)>0)
                            cat('Pruning', length(toprune), 'walks\n')
                          allps = allps[-toprune]
                        }
                    }
                }

              if (length(allps)>0)
                win = streduce(unlist(allps), 0)
#                win = streduce(unlist(allps), sum(width(unlist(allps)))*0)

              values(allps) = NULL
              out$allpaths = allps
              out$allpaths.og = allps.og ## untouched all.paths if we want to reheat eg after computing 10X support
              gt.walk = gTrack(out$allpaths, draw.paths = T,border = NA, angle = 0, ywid = 0.5, height = 20, labels.suppress.gr = T)
              gt.walk$path.cex.arrow = 0
              gt.walk$path.stack.x.gap = 1e6
              out$gtrack.allpaths = c(
                  gt.walk,
                  td.seg,
                  td.rg)
              pdf(outfile.allpaths.pdf, height = 30, width = 24)
              gTrack::plot(out$gtrack.allpaths,
                      windows = win, links = kag$junctions)
              dev.off()
              out$README = paste(out$README, 'allpaths= all paths through windows (not just optimal ones), td.allpaths = gTrack object of plot of all paths')
            }

          ## if junction.ix was specified then label which positions in the walks represent the rearrangement junction
          if (!is.null(junction.ix) & length(out$allpaths)>0)
            {
              allps = out$allpaths
              allps.u = grl.unlist(allps)
              allps.u$ix.s = gr.match(gr.start(allps.u, ignore.strand = F), starts, ignore.strand = F)
              allps.u$ix.e = gr.match(gr.end(allps.u, ignore.strand = F), ends, ignore.strand = F)
              allps.ixs = split(allps.u$ix.s, allps.u$grl.ix) ## start indices of walk intervals in sol$segstats
              allps.ixe = split(allps.u$ix.e, allps.u$grl.ix) ## end indices of walks intervals in sol$segstats
              allps.ix.pairs = sapply(mapply(function(x,y) if (length(x)<=1) NULL else which(paste(x[-length(x)], y[-1], sep = ',') %in% junc.pair), allps.ixe, allps.ixs, SIMPLIFY = F), paste, collapse = ',')
              values(allps)$junction.id = names(junction.ix)[i]
              values(allps)$junction.ix = allps.ix.pairs
              out$allpaths = allps
            }

          if (length(out$allpaths)>0)
            {
              values(out$allpaths)$string = grl.string(out$allpaths)
              values(out$allpaths)$wid = sapply(out$allpaths, function(x) sum(width(x)))
              values(out$allpaths)$wids = sapply(out$allpaths, function(x) paste(width(x), collapse = ','))
              write.tab(as.data.frame(values(out$allpaths)), outfile.allpaths.txt)
            }

          saveRDS(out, outfile.rds)
          return(out)
        }, mc.cores = mc.cores)
    }

  ## awkward workaround to limit the number of processors Cplex will gobble up
  if (customparams)
    {
      system(paste('rm', param.file))
      Sys.setenv(ILOG_CPLEX_PARAMETER_FILE='')
      if (verbose)
        {
          jmessage('Finished')
        }
    }

  return(out)
}



###########################
#' @name proximity
#' @export
#' @rdname internal
#' proximity
#'
#' Takes a set of n "query" elements (GRanges object, e.g. genes) and determines their proximity to m "subject" elements
#' (GRanges object, e.g. regulatory elements) subject to set of rearrangement adjacencies (GRangesList with width 1 range pairs)
#'
#' This analysis makes the (pretty liberal) assumption that all pairs of adjacencies that can be linked on a karyograph path are in
#' cis (i.e. share a chromosome) in the tumor genome.
#'
#' @param query GRanges of "intervals of interest" eg regulatory elements
#' @param subject GRanges of "intervals of interest" eg genes
#' @param ra GRangesList of junctions (each a length 2 GRanges, similar to input to karyograph)
#' @param jab existing JaBbA object (overrides ra input)
#' @param verbose logical flag
#' @param mc.cores how many cores (default 1)
#' @param max.dist maximum genomic distance to store and compute (1MB by default) should the maximum distance at which biological interactions may occur
#' @return
#' list of n x m sparse distance matrices:
#' $ra = subject-query distance in the rearranged genome for all loci < max.dist in tumor genome
#' $wt = subject-query distance in the reference genome for all loci < max.dist in tumor genome
#' $rel = subject-query distance in ra relative to wild type for above loci
#' NOTE: values x_ij in these matrices should be interpreted with a 1e-9 offset to yield the actual value y_ij
#' i.e. y_ij = x_ij-1e-9, x_ij>0, y_ij = NA otherwise (allows for sparse encoding of giant matrices)
############################################
proximity = function(query, subject, ra = GRangesList(), jab = NULL, verbose = F, mc.cores = 1,
  max.dist = 1e6 ## max distance to store / compute in the output matrix.cores
  )
  {
    if (!is.null(jab))
    {
        nnab = which(!ifelse(is.na(jab$ab.edges[,3,1]), TRUE, jab$ab.edges[,3,1]==0))
        ix = nnab[which(jab$adj[jab$ab.edges[nnab,1:2,1]]>0)]
        if (length(ix)>0)
          {
            ra1 = gr.flipstrand(gr.end(jab$segstats[jab$ab.edges[ix,1,1]], 1, ignore.strand = F))
            ra2 = gr.start(jab$segstats[jab$ab.edges[ix,2,1]], 1, ignore.strand = F)
            ra1 = GenomicRanges::shift(ra1, ifelse(as.logical(strand(ra1)=='+'), -1, 0))
            ra2 = GenomicRanges::shift(ra2, ifelse(as.logical(strand(ra2)=='+'), -1, 0))
            ra = grl.pivot(GRangesList(ra1,ra2))
          }
      }

    if (length(ra)==0)
      return(list())

    if (length(query)==0 | length(subject)==0)
        return(list())

    if (is.null(names(query)))
        names(query) = 1:length(query)

    if (is.null(names(subject)))
        names(subject) = 1:length(subject)

    query.nm = names(query);
    subject.nm = names(subject);

    query = query[, c()]
    subject = subject[, c()]

    query$id = 1:length(query)
    subject$id = 1:length(subject)

    qix.filt = gr.in(query, unlist(ra)+max.dist) ## to save time, filter only query ranges that are "close" to RA's
    query = query[qix.filt]

    six.filt = gr.in(subject, unlist(ra)+max.dist) ## to save time, filter only query ranges that are "close" to RA's
    subject = subject[six.filt]

    if (length(query)==0 | length(subject)==0)
        return(list())

    query$type = 'query'
    subject$type = 'subject'

    gr = gr.fix(c(query, subject))

    kg = JaBbA:::karyograph(ra, gr)

    ## node.start and node.end delinate the nodes corresponding to the interval start and end
    ## on both positive and negative tiles of the karyograph
    gr$node.start = gr$node.end = gr$node.start.n = gr$node.end.n = NA;

    ## start and end indices of nodes
    tip = which(as.character(strand(kg$tile))=='+')
    tin = which(as.character(strand(kg$tile))=='-')
    gr$node.start = tip[gr.match(gr.start(gr,2), gr.start(kg$tile[tip]))]
    gr$node.end = tip[gr.match(GenomicRanges::shift(gr.end(gr,2),1), gr.end(kg$tile[tip]))]
    gr$node.start.n = tin[gr.match(GenomicRanges::shift(gr.end(gr,2),1), gr.end(kg$tile[tin]))]
    gr$node.end.n = tin[gr.match(gr.start(gr,2), gr.start(kg$tile[tin]))]

    if (any(nix <<- is.na(gr$node.start)))
        gr$node.start[nix] = tip[nearest(gr.start(gr[nix],2), gr.start(kg$tile[tip]))]

    if (any(nix <<- is.na(gr$node.end)))
        gr$node.end[nix] = tip[nearest(GenomicRanges::shift(gr.end(gr[nix],2),1), gr.end(kg$tile[tip]))]


    if (any(nix <<- is.na(gr$node.end.n)))
      gr$node.end.n[nix] = tin[nearest(gr.start(gr[nix],2), gr.start(kg$tile[tin]))]

    if (any(nix <<- is.na(gr$node.start.n)))
        gr$node.start.n[nix] = tin[nearest(GenomicRanges::shift(gr.end(gr[nix],2),1), gr.end(kg$tile[tin]))]


#    gr$node.start = gr.match(gr.start(gr-1,2), gr.start(kg$tile))
#    gr$node.end = suppressWarnings(gr.match(gr.end(gr+1,2), gr.end(kg$tile)))

    ## so now we build distance matrices from query ends to subject starts
    ## and subject ends to query starts

    ## so for each query end we will find the shortest path to all subject starts
    ## and for each query start we will find the shortest.path from all subject ends
    ix.query = which(gr$type == 'query')
    ix.subj = which(gr$type == 'subject')

    node.start = gr$node.start
    node.end = gr$node.end
    node.start.n = gr$node.start.n
    node.end.n = gr$node.end.n

    w = width(kg$tile)

    E(kg$G)$weight = width(kg$tile)[E(kg$G)$to]

    ## ix.query and ix.subj give the indices of query / subject in gr
    ## node.start, node.end map gr to graph node ids
    ##
    ## these matrices are in dimensions of query and subject, and will hold the pairwise distances between
    ##
    D.rel = D.ra = D.ref = D.which = Matrix(data = 0, nrow = length(ix.query), ncol = length(ix.subj))

    ## "reference" graph (missing aberrant edges)
    G.ref = subgraph.edges(kg$G, which(E(kg$G)$type == 'reference'), delete.vertices = F)

    EPS = 1e-9

    tmp = mclapply(ix.query, function(i)
      {
        if (verbose)
          cat('starting interval', i, 'of', length(ix.query), '\n')

        ## D1 = shortest query to subject path, D2 = shortest subject to query path, then take shortest of D1 and D2
        ## for each path, the edge weights correspond to the interval width of the target node, and to compute the path
        ## length we remove the final node since we are measuring the distance from the end of the first vertex in the path
        ## to the beginning of the final vertex

        u.node.start = unique(node.start[ix.subj]) ## gets around annoying igraph::shortest.path issue (no dups allowed)
        u.node.end = unique(node.end[ix.subj])

        uix.start = match(node.start[ix.subj], u.node.start)
        uix.end = match(node.end[ix.subj], u.node.end)

        tmp.D1 = (shortest.paths(kg$G, node.end[i], u.node.start, weights = E(kg$G)$weight, mode = 'out') - w[u.node.start])[uix.start]
        tmp.D2 = (shortest.paths(kg$G, node.start[i], u.node.end, weights = E(kg$G)$weight, mode = 'in') - w[node.start[i]])[uix.end]
        tmp.D3 = (shortest.paths(kg$G, node.end.n[i], u.node.start, weights = E(kg$G)$weight, mode = 'out') - w[u.node.start])[uix.start]
        tmp.D4 = (shortest.paths(kg$G, node.start.n[i], u.node.end, weights = E(kg$G)$weight, mode = 'in') - w[node.start.n[i]])[uix.end]
        tmp.D = pmin(tmp.D1, tmp.D2, tmp.D3, tmp.D4)
        ix = Matrix::which(tmp.D<max.dist)
        D.ra[i, ix] = tmp.D[ix]+EPS
        D.which[i, ix] = apply(cbind(tmp.D1[ix], tmp.D2[ix], tmp.D3[ix], tmp.D4[ix]), 1, which.min)

        u.node.start = unique(node.start[ix.subj][ix]) ## gets around annoying igraph::shortest.path issue (no dups allowed)
        u.node.end = unique(node.end[ix.subj][ix])

        uix.start = match(node.start[ix.subj][ix], u.node.start)
        uix.end = match(node.end[ix.subj][ix], u.node.end)

        tmp.D1 = (shortest.paths(G.ref, node.end[i], u.node.start, weights = E(G.ref)$weight, mode = 'out') - w[u.node.start])[uix.start]
        tmp.D2 = (shortest.paths(G.ref, node.start[i], u.node.end, weights = E(G.ref)$weight, mode = 'in') - w[node.start[i]])[uix.end]
        tmp.D3 = (shortest.paths(G.ref, node.end.n[i], u.node.start, weights = E(G.ref)$weight, mode = 'out') - w[u.node.start])[uix.start]
        tmp.D4 = (shortest.paths(G.ref, node.start.n[i], u.node.end, weights = E(G.ref)$weight, mode = 'in') - w[node.start.n[i]])[uix.end]
        tmp.D = pmin(tmp.D1, tmp.D2, tmp.D3, tmp.D4)
        D.ref[i, ix] = tmp.D+EPS

        ## if subject and query intersect (on the reference) then we count both RA and Ref distance as 0
        ## (easier to do a simple range query here)
        ix.zero = gr.in(subject[ix], query[i])
        if (any(ix.zero))
          {
            D.ra[i, ix[ix.zero]] = 0
            D.ref[i, ix[ix.zero]] = 0
          }

        D.rel[i, ix] = ((D.ra[i, ix]-EPS) / (D.ref[i, ix]-EPS)) + EPS

        if (verbose)
           cat('finishing interval', i, 'of', length(ix.query), ':', paste(round(D.rel[i, ix],2), collapse = ', '), '\n')

        return(list(D.rel = D.rel, D.ref = D.ref, D.ra = D.ra, D.which = D.which))
       }, mc.cores = mc.cores)

    for (i in 1:length(tmp))
    {
        if (class(tmp[[i]]) != 'list')
            warning(sprintf('Query %s failed', ix.query[i]))
        else
        {
            D.rel = D.rel + tmp[[i]]$D.rel
            D.ra = D.ra + tmp[[i]]$D.ra
            D.ref = D.ref + tmp[[i]]$D.ref
            D.which = D.which + tmp[[i]]$D.which
        }
    }

    ## "full" size matrix
    rel = ra = ref = ra.which =
      Matrix(data = 0, nrow = length(qix.filt), ncol = length(six.filt), dimnames = list(dedup(query.nm), dedup(names(subject.nm))))
    rel[qix.filt, six.filt] = D.rel
    ra[qix.filt, six.filt] = D.ra
    ref[qix.filt, six.filt] = D.ref
    ra.which[qix.filt, six.filt] = D.which

    ## summary is data frame that has one row for each query x subject pair, relative distance, ra distance, and absolute distance
    tmp = which(rel!=0, arr.ind = T)
    colnames(tmp) = c('i', 'j');
    sum = as.data.frame(tmp)

    if (!is.null(query.nm))
      sum$query.nm = query.nm[sum$i]

    if (!is.null(subject.nm))
      sum$subject.nm = subject.nm[sum$j]

    sum$rel = rel[tmp]
    sum$ra = ra[tmp]
    sum$wt = ref[tmp]

    sum = sum[order(sum$rel), ]
    sum = sum[sum$rel<1, ] ## exclude those with rel == 1

    ## reconstruct paths
    vix.query = matrix(NA, nrow = length(qix.filt), ncol = 4, dimnames = list(NULL, c('start', 'end', 'start.n', 'end.n')))
    vix.subject = matrix(NA, nrow = length(six.filt), ncol = 4, dimnames = list(NULL, c('start', 'end', 'start.n', 'end.n')))
    vix.query[qix.filt, ] = cbind(values(gr)[ix.query, c('node.start')], values(gr)[ix.query, c('node.start')], values(gr)[ix.query, c('node.start.n')], values(gr)[ix.query, c('node.end.n')])
    vix.subject[six.filt] = cbind(values(gr)[ix.subj, c('node.start')], values(gr)[ix.subj, c('node.start')], values(gr)[ix.subj, c('node.start.n')], values(gr)[ix.subj, c('node.end.n')])


    if (nrow(sum)==0)
        return(NULL)

    sum.paths = mcmapply(function(x, y, i)
    {
        if (verbose)
        jmessage('path ', i, ' of ', nrow(sum), '\n')
        if ((ra.which[x, y]) == 1)
            get.shortest.paths(kg$G, vix.query[x, 'end'], vix.subject[y, 'start'], weights = E(kg$G)$weight, mode = 'out')$vpath[[1]]
        else if ((ra.which[x, y]) == 2)
            rev(get.shortest.paths(kg$G, vix.query[x, 'start'], vix.subject[y, 'end'], weights = E(kg$G)$weight, mode = 'in')$vpath[[1]])
        else if ((ra.which[x, y]) == 3)
            get.shortest.paths(kg$G, vix.query[x, 'end.n'], vix.subject[y, 'start'], weights = E(kg$G)$weight, mode = 'out')$vpath[[1]]
        else if ((ra.which[x, y]) == 4)
            rev(get.shortest.paths(kg$G, vix.query[x, 'start.n'], vix.subject[y, 'end'], weights = E(kg$G)$weight, mode = 'in')$vpath[[1]])
    }, sum$i, sum$j, 1:nrow(sum), SIMPLIFY = F, mc.cores = mc.cores)

#    sum$paths = lapply(sum.paths, function(x) x[-c(1, length(x))])
    sum$paths = sum.paths
    sum$ab.edges = lapply(sum.paths, function(p) setdiff(E(kg$G, path = p)$bp.id, NA))

    return(list(sum = sum, rel = rel, ra = ra, wt = ref, G = kg$G, G.ref = G.ref, tile = kg$tile, vix.query = vix.query, vix.subject = vix.subject))
  }



####################################################################
#' @name jabba.melt
#' @export
#' @rdname internal
#' jabba.melt
#'
#' @details
#' melt JaBbA graph into "events" that decompose the total ploidy into amplicons (or deleticons, if anti = TRUE)
#' Each amplicons / deleticon is flanked by either (1) junctions (2) loose ends or (3) chromosome ends / telomeres
#'
#'
#' @param jab JaBbA object "undigested"
#' @param kag karyograph (original karyograph input to JaBbA), if NULL then will "redigest" JaBbA object
#' @param verbose logical flag
#' @param keep.all keep.all (default TRUE) whether to keep 0 copy junctions or collapse segments across these as well
####################################################################
jabba.melt = function(jab, anti = FALSE, verbose = FALSE, mc.cores = 1, max.del = 10)
    {
        abs = rbind(jab$ab.edges[,,1], jab$ab.edges[,,2])[, 1:2]
        abs = abs[!is.na(abs[,1]), ]
        adj = data.table(i = abs[,1], j = abs[,2])
        adj[, cn := jab$adj[cbind(i, j)]]
        setkeyv(adj, c('i', 'j'))
        junc.right = adj[, sum(cn), keyby = i]
        junc.left = adj[, sum(cn), keyby = j]

        gr = gr2dt(jab$segstats)[, seg.id := 1:length(seqnames)][loose == FALSE & strand == '+' & !is.na(cn), ]
        gr[, cluster := {
            tmp = rle(cn)
            rep(paste(seqnames[1], 1:length(tmp$values), sep = '.'), tmp$length)
        }, by = seqnames]

        lix = jab$segstats$loose
        gr[, loose.right := Matrix::rowSums(as.matrix(jab$adj[seg.id, lix, drop = FALSE]))]
        gr[, loose.left := Matrix::colSums(as.matrix(jab$adj[lix, seg.id, drop = FALSE]))]

        gr[, loose.right := Matrix::rowSums(as.matrix(jab$adj[seg.id, lix, drop = FALSE]))]
        gr[, loose.left := Matrix::colSums(as.matrix(jab$adj[lix, seg.id, drop = FALSE]))]

        if (nrow(junc.left)>0)
            {
                gr[, ab.left := ifelse(is.na(junc.left[list(seg.id), V1]), 0, junc.left[list(seg.id), V1])]
                gr[, ab.left.ix := gr[ , adj[, ][cn>0, i[1], keyby = j][list(seg.id), V1]]]
            }
        else
            gr[, ab.left := 0]

        if (nrow(junc.right)>0)
            {
                gr[, ab.right := ifelse(is.na(junc.right[list(seg.id), V1]), 0, junc.right[list(seg.id), V1])]
                gr[, ab.right.ix := gr[ , adj[, ][cn>0, j[1], keyby = i][list(seg.id), V1]]]
            }

        else
            gr[, ab.right := 0]

        gr[, id := 1:nrow(gr)]
        setkey(gr, id)

        .flip = function(gr)
            {
                gr = as.data.table(as.data.frame(gr))
                if (nrow(gr) <= 1)
                    return(gr)
                max.cn = gr[, pmin(max(cn), max.del)] ## cap max cn for "anti analysis" since we may not care about dels in high copy contexts, and want to avoid edge jabba cases where a region got infinite copy number

                ## flip copy number upside down
                gr[, cn := max.cn-cn]

                ## assign right.ab[i] and right.loose[i] to left.ab[i+1]

                for (i in 1:(nrow(gr)-1))
                    {
                        tmp.loose = gr[i+1, loose.left]
                        tmp.ab = gr[i+1, ab.left]
                        gr[i+1, loose.left := gr$loose.right[i]]
                        gr[i+1, ab.left := gr$ab.right[i]]
                        gr[i, loose.right := tmp.loose]
                        gr[i, ab.right := tmp.ab]
                    }
                return(gr)
            }

        if (anti)
            gr = rbindlist(lapply(split(gr, gr$seqnames), .flip))

        .fun = function(gr)
            {
                ## this assumes the input is from one seqneme
                gr = as.data.table(as.data.frame(gr))
                cl = split(gr$id, gr$cluster)
                ix = order(-gr[, cn])

                if (nrow(gr)==0)
                    return(NULL)

                gr[, done := FALSE]
                gr[1, loose.left := cn] ## give telomeres temporary loose ends
                gr[nrow(gr), loose.right := cn]

                out = data.table(seqnames = rep(as.character(NA), nrow(gr)), start = as.numeric(NA), end = as.numeric(NA), cn = as.numeric(NA),
                    cn.og = as.numeric(NA),
                    left.ab = as.numeric(NA),
                    right.ab = as.numeric(NA),
                    left.loose = as.numeric(NA),
                    right.loose = as.numeric(NA),
                    left.ix = as.numeric(NA),
                    right.ix = as.numeric(NA))
                out.i = 1

                for (i in gr$id[ix])
                    {
                        if (verbose)
                            cat('.')
                        setkey(gr, id)
                        if (!gr[list(i), done])
                            {
#                                if (i==21)
 #                                   browser()
                                this.cl = cl[[gr[list(i), cluster]]]
                                out[out.i, ":="(left.ix = this.cl[1],
                                                right.ix = this.cl[length(this.cl)])]
                                left.cl = out[out.i, gr[list(c(left.ix-1, left.ix)), setdiff(cluster, NA)[1]]]
                                right.cl = out[out.i, gr[list(c(right.ix+1, right.ix)), setdiff(cluster, NA)[1]]]
                                ## neighbor CN is max of left and right
                                left.cn = out[out.i, gr[list(left.ix-1), cn]]
                                right.cn = out[out.i, gr[list(right.ix+1), cn]]
                                neighbor.cn = max(c(0, left.cn,  right.cn), na.rm = TRUE)

                                out.cn = gr[list(i), cn - neighbor.cn]

                                out[out.i, ":="(cn = out.cn, seqnames = gr[list(i), seqnames], start = gr[list(left.ix), start], end = gr[list(right.ix), end], cn.og = gr[list(i), cn])]
                                tmp = out[out.i, pmin(gr[list(left.ix), ab.left], out.cn)]
                                out[out.i, ':='(left.ab = tmp, left.loose = out.cn - tmp)]

                                tmp = out[out.i, pmin(gr[list(right.ix), ab.right], out.cn)]
                                out[out.i, ':='(right.ab = tmp, right.loose = out.cn - tmp)]

                                        # update gr
                                        # subtracting the copies we have assigned to this outgoing event

                                gr[list(out[out.i, left.ix]), ab.left := ab.left - out[out.i, left.ab]]
                                gr[list(out[out.i, left.ix]), loose.left := loose.left - out[out.i, left.loose]]
                                gr[list(out[out.i, right.ix]), ab.right := ab.right - out[out.i, right.ab]]
                                gr[list(out[out.i, right.ix]), loose.right := loose.right - out[out.i, right.loose]]
                                gr[list(this.cl), ":="(cn = cn - out.cn)]

                                ## update clusters - merge left and right depending on which cn are equals
                                ## the new copy number for this interval should be equal to the left or right cn
                                this.cli = gr[list(i), cluster]
                                left.cn = ifelse(is.na(left.cn), -Inf, left.cn)
                                right.cn = ifelse(is.na(right.cn), +Inf, right.cn)
                                if (left.cn == right.cn)
                                    {
                                        gr[list(this.cl), ":="(cluster = left.cl)]
                                        gr[list(cl[[right.cl]]), ":="(cluster = left.cl)]
                                        cl[[left.cl]] = do.call('c', cl[unique(c(left.cl, this.cli, right.cl))])
                                    }
                                else if (left.cn == gr[list(i), cn])
                                    {
                                        gr[list(this.cl), ":="(cluster = left.cl)]
                                        cl[[left.cl]] = do.call('c', cl[unique(c(left.cl, this.cli))])
                                    }
                                else
                                    {
                                        gr[list(this.cl), ":=" (cluster = right.cl)]
                                        cl[[right.cl]] = do.call('c', cl[unique(c(this.cli, right.cl))])
                                    }
                                gr[list(this.cl), done := TRUE] ## don't go back to any interval in this cluster
                                ## advance output

                                out.i = out.i +1
                            }
                    }

                out[, left.tel := 0]
                out[, right.tel := 0]
                out[left.ix == gr[, id[1]], ":="(left.tel = left.loose, left.loose = 0)]
                out[right.ix == gr[, id[nrow(gr)]], ":="(right.tel = right.loose, right.loose = 0)]
                return(out)
            }

        tmp = mclapply(split(gr, gr$seqnames), .fun, mc.cores = mc.cores)
        out = rbindlist(tmp)[!is.na(cn), ][cn>0, ]
        out = seg2gr(out, seqlengths = seqlengths(jab$segstats))
        out$cn.max = gr.val(out, jab$segstats, 'cn', FUN = max, weighted = FALSE)$cn
        out$cn.min = gr.val(out, jab$segstats, 'cn', FUN = min, weighted = FALSE)$cn
        return(out)
    }

####################################################
#' @name jabba.hood
#' @export
#' @rdname internal
#' jabba.hood
#'
#' Given JaBbA  object
#' and seed window "win", outputs a reduced set of windows within neighborhoof of n coordinate (ork nodes)
#' within the seed region(s) on the graph (only includes edges with weight !=0)
#'
#' @param jab JaBbA object
#' @param win GRanges of window of interest
#' @param d = distance in coordinates on graph
#' @param k Neighborhood on graph around window of interest to query
#' @param pad pad level at which to collapse nearly reference adjacent intervals
#' @return a reduced set of windows within neighborhood k
#' of seed on the graph (only includes edges with weight !=0)
#########x############################################
jabba.hood = function(jab, win, d = 0, k = NULL, pad = 0, ignore.strand = TRUE, bagel = FALSE, verbose = FALSE)
{
    if (ignore.strand)
        win = gr.stripstrand(win)

    if (is.null(k)) ## use distance
        {

            ss = tryCatch(c(jab$segstats[jab$segstats$loose == FALSE, c()], win[, c()]), error = function(e) NULL)

            if (is.null(ss))
                ss = grbind(c(jab$segstats[jab$segstats$loose == FALSE, c()], win[, c()]))

            if (ignore.strand)
                ss = gr.stripstrand(ss)

            ss = disjoin(ss)
            win = gr.findoverlaps(ss, win, ignore.strand = ignore.strand)

            seg.s = suppressWarnings(gr.start(ss, ignore.strand = TRUE))
            seg.e = suppressWarnings(gr.end(ss, ignore.strand = TRUE))
            D.s = suppressWarnings(jabba.dist(jab, win, seg.s, verbose = verbose))
            D.e = suppressWarnings(jabba.dist(jab, win, seg.e, verbose = verbose))

            min.s = apply(D.s, 2, min, na.rm = TRUE)
            min.e = apply(D.e, 2, min, na.rm = TRUE)
            s.close = min.s<=d
            e.close = min.e<=d

            ## now for all "left close" starts we add whatever distance to that point + pad
            gr.start(ss)[s.close]

            out = GRanges()
            if (any(s.close))
                out = c(out, GenomicRanges::flank(seg.s[s.close], -(d-min.s[s.close])))

            if (any(e.close))
                out = c(out, GenomicRanges::shift(flank(seg.e[e.close], d-min.e[e.close]),1))

            if (!bagel)
                out = streduce(c(win[, c()], out[, c()]))

            return(streduce(out, pad))
        }
    else ## use graph connections
        {
            G = tryCatch(graph.adjacency(jab$adj!=0), error = function(e) NULL)

            ix = which(jab$segstats %^% win)
            if (is.null(G)) ## sometimes igraph doesn't like Matrix
                G = graph.edgelist(which(jab$adj!=0, arr.ind = TRUE))
            vix = unique(unlist(neighborhood(G, ix, order = k)))
            return(streduce(jab$segstats[vix], pad))
        }
}


######################################################
#' @name jabba.dist
#' @export
#' @rdname internal
#' jabba.dist
#'
#' Computes distance between pairs of intervals on JaBbA graph
#'
#' Given "jabba" object and input granges gr1 and gr2 of (signed) intervals
#'
#'
#' @param jab JaBbA object
#' @param gr1 interval set 1 GRanges
#' @param gr2 interval set 2 GRanges
#' @param matrix flag whteher to output a matrix
#' @param max.dist numeric (default = Inf), if non-infinity then output will be a sparse matrix with all entries that are greater than max.dist set to zero
#' @return a length(gr1) x length(gr2) matrix whose entries ij store the distance between
#' the 3' end of gr1[i] and 5' end of gr2[j]
#######################################################
jabba.dist = function(jab, gr1, gr2,
                      matrix = T, ## if false then will return a data frame with fields $i $j $dist specifying distance between ij pairs
                      directed= FALSE, ## flag specifying whether we are computing a "directed distance" across only paths FROM gr1 TO gr2 on graph (ie gr2-->gr1 paths do not count
                      max.dist = Inf, ## if max.dist is not Inf then a sparse matrix will be returned that has 0 at all locations greater than max.dist
                      include.internal = TRUE, ## includes internal connections eg if a junction lies inside a feature then that feature is "close" to another feature
                      verbose = FALSE,
                      EPS = 1e-9  ## the value used for "real 0" if a sparse matrix is returned
  )
{
    if (verbose)
        now = Sys.time()

    intersect.ix = gr.findoverlaps(gr1, gr2, ignore.strand = FALSE)

    ngr1 = length(gr1)
    ngr2 = length(gr2)

    if (is.null(jab$segstats))
      tiles = jab$tile
    else
      tiles = jab$segstats;

    if (is.null(jab$G))
      G = graph.adjacency(jab$adj!=0)
    else
      G = jab$G

    ## keep track of original ids when we collapse
    gr1$id = 1:length(gr1)
    gr2$id = 1:length(gr2)

    ## check for double stranded intervals
    ## add corresponding nodes if present
    if (any(ix <- as.logical( strand(gr1)=='*')) )
    {
        strand(gr1)[ix] = '+'
        gr1 = c(gr1, gr.flipstrand(gr1[ix]))
    }

    if (any(ix <- as.logical( strand(gr2)=='*')))
    {
        strand(gr2)[ix] = '+'
        gr2 = c(gr2, gr.flipstrand(gr2[ix]))
    }

    ## expand nodes by jabba model to get internal connectivity
    if (include.internal)
    {
        gr1 = gr1[, 'id'] %**% jab$segstats
        gr2 = gr2[, 'id'] %**% jab$segstats
    }

    if (verbose)
        {
            jmessage('Finished making gr objects')
            print(Sys.time() -now)
        }

    tmp = igraph::get.edges(G, E(G))
    E(G)$from = tmp[,1]
    E(G)$to = tmp[,2]
    E(G)$weight = width(tiles)[E(G)$to]

    gr1.e = gr.end(gr1, ignore.strand = FALSE)
    gr2.s = gr.start(gr2, ignore.strand = FALSE)


    if (!directed)
        {
            gr1.s = gr.start(gr1, ignore.strand = FALSE)
            gr2.e = gr.end(gr2, ignore.strand = FALSE)
        }

    gr1.e$ix = gr.match(gr1.e, tiles, ignore.strand = F) ## graph node corresponding to end of gr1.ew
    gr2.s$ix= gr.match(gr2.s, tiles, ignore.strand = F) ## graph node corresponding to beginning of gr2

    if (!directed)
        {
            gr1.s$ix = gr.match(gr1.s, tiles, ignore.strand = F) ## graph node corresponding to end of gr1.ew
            gr2.e$ix= gr.match(gr2.e, tiles, ignore.strand = F) ## graph node corresponding to beginning of gr2
        }

    ## 3' offset from 3' end of query intervals to ends of jabba segs  to add / subtract to distance when query is in middle of a node
    off1 = ifelse(as.logical(strand(gr1.e)=='+'), end(tiles)[gr1.e$ix]-end(gr1.e), start(gr1.e) - start(tiles)[gr1.e$ix])
    off2 = ifelse(as.logical(strand(gr2.s)=='+'), end(tiles)[gr2.s$ix]-end(gr2.s), start(gr2.s) - start(tiles)[gr2.s$ix])

    ## reverse offset now calculate 3' offset from 5' of intervals
    if (!directed)
        {
            off1r = ifelse(as.logical(strand(gr1.s)=='+'), end(tiles)[gr1.s$ix]-start(gr1.s), end(gr1.s) - start(tiles)[gr1.s$ix])
            off2r = ifelse(as.logical(strand(gr2.e)=='+'), end(tiles)[gr2.e$ix]-start(gr2.e), end(gr2.e) - start(tiles)[gr2.e$ix])
        }

    ## compute unique indices for forward and reverse analyses
    uix1 = unique(gr1.e$ix)
    uix2 = unique(gr2.s$ix)

    if (!directed)
        {
            uix1r = unique(gr1.s$ix)
            uix2r = unique(gr2.e$ix)
        }

    ## and map back to original indices
    uix1map = match(gr1.e$ix, uix1)
    uix2map = match(gr2.s$ix, uix2)

    if (!directed)
        {
            uix1mapr = match(gr1.s$ix, uix1r)
            uix2mapr = match(gr2.e$ix, uix2r)
        }

    self.l = which(Matrix::diag(jab$adj)>0)

    if (verbose)
        {
            jmessage('Finished mapping gr1 and gr2 objects to jabba graph')
            print(Sys.time() -now)
        }

    if (is.infinite(max.dist)) ## in this case we do not bother making sparse matrix and can compute distances very quickly with one call to shortest.paths
    {
        ## need to take into account forward and reverse scenarios of "distance" here
        ## ie upstream and downstream connections between query and target
        ## edges are annotated with width of target

        ## so for "downstream distance"  we are getting matrix of shortest paths between from uix1 and uix2 node pair
        ## and then correcting those distances by (1) adding the 3' offset of uix1 from its node
        ## and (2) subtracting the 3' offset of uix2
        Df = sweep(
            sweep(
                shortest.paths(G, uix1, uix2, weights = E(G)$weight, mode = 'out')[uix1map, uix2map, drop = F],
                1, off1, '+'), ## add uix1 3' offset to all distances
            2, off2, '-') ## subtract uix2 3' offset to all distances


        if (!directed)
            {
                ## now looking upstream - ie essentially flipping edges on our graph - the edge weights
                ## now represent "source" node widths (ie of the flipped edges)
                                        # need to correct these distances by (1) subtracting 3' offset of uix1 from its node
                ## and (2) adding the 3' offset of uix2
                ## and using the reverse indices
                Dr = sweep(
                    sweep(
                        t(shortest.paths(G, uix2r, uix1r, weights = E(G)$weight, mode = 'out'))[uix1mapr, uix2mapr, drop = F],
                        1, off1r, '-'), ## substract  uix1 offset to all distances but subtract weight of <first> node
                    2, off2r , '+') ## add uix2 offset to all distances

                Df2 = sweep(
                    sweep(
                        shortest.paths(G, uix1r, uix2, weights = E(G)$weight, mode = 'out')[uix1mapr, uix2map, drop = F],
                        1, off1r, '+'), ## add uix1 3' offset to all distances
                    2, off2, '-') ## subtract uix2 3' offset to all distances

                Dr2 = sweep(
                    sweep(
                        t(shortest.paths(G, uix2r, uix1, weights = E(G)$weight, mode = 'out'))[uix1map, uix2mapr, drop = F],
                        1, off1, '-'), ## substract  uix1 offset to all distances but subtract weight of <first> node
                    2, off2r , '+') ## add uix2 offset to all distances
                D = pmin(abs(Df), abs(Dr), abs(Df2), abs(Dr2))
            }
        else
            D = Df

        # then we do the same thing but flipping uix1r vs uix


        if (verbose)
            {
                jmessage('Finished computing distances')
                print(Sys.time() -now)
            }


        ## take care of edge cases where ranges land on the same node, since igraph will just give them "0" distance
        ## ij contains pairs of gr1 and gr2 indices that map to the same node
        ij = as.matrix(merge(cbind(i = 1:length(gr1.e), nid = gr1.e$ix), cbind(j = 1:length(gr2.s), nid = gr2.s$ix)))

        ## among ij pairs that land on the same (strand of the same) node
        ##
        ## several possibilities:
        ## (1) if gr1.e[i] < gr2.s[j] then keep original distance (i.e. was correctly calculated)
        ## (2) if gr1.e[i] > gr2.s[j] then either
        ##   (a) check if there is a self loop and adjust accordingly (i.e. add back width of current tile)
        ##   (b) PITA case, compute shortest distance from i's child(ren) to j

        if (nrow(ij)>0)
          {
            ## rix are present
              rix = as.logical((
                  (as.logical( strand(gr1)[ij[,'i']] == '+' ) &
                   as.logical( strand(gr2)[ij[,'j']] == '+' ) &
                   end(gr1)[ij[,'i']] <= start(gr2[ij[,'j']])) |
                  ( as.logical( strand(gr1)[ij[,'i']] == '-' ) &
                    as.logical( strand(gr2)[ij[,'j']] == '-' ) &
                    start(gr1)[ij[,'i']] >= end(gr2)[ij[,'j']])))

              ij = ij[!rix, , drop = F] ## NTD with rix == TRUE these since they are calculated correctly

            if (nrow(ij)>0) ## any remaining will either be self loops or complicated loops
              {
                selfix = (ij[, 'nid'] %in% self.l)

                if (any(selfix)) ## correct distance for direct self loops (add back width of current node)
                  D[ij[selfix, c('i', 'j'), drop = F]]  = D[ij[selfix, c('i', 'j'), drop = F]] + width(tiles)[ij[selfix, 'nid']]

                ij = ij[!selfix, , drop = F]

                if (nrow(ij)>0) ## remaining are pain in the ass indirect self loops
                  {
                    ch = G[[ij[, 'nid']]] ## list of i nodes children for all remaining ij pairs
                    chu = munlist(ch) ## unlisted children, third column are the child id's, first column is the position of nrix

                    if (ncol(chu)>1)
                        {

                    ## now find paths from children to corresponding j
                            epaths = suppressWarnings(get.shortest.paths(G, chu[, 3], ij[chu[,'ix'], 'nid'], weights = E(G)$weight, mode = 'out', output = 'epath')$epath)
                            epathw = sapply(epaths, function(x,w) if (length(x)==0) Inf else sum(w[x]), E(G)$weight) ## calculate the path weights
                            epathw = epathw + width(tiles)[chu[, 3]] + off1[ij[chu[, 'ix'], 'i']] + off2[ij[chu[,'ix'], 'j']] - width(tiles)[ij[chu[, 'ix'], 'nid']]

                            ## aggregate (i.e. in case there are multiple children per node) by taking min width
                            D[ij[, c('i', 'j'), drop = F]] = vaggregate(epathw, by = list(chu[, 'ix']), min)[as.character(1:nrow(ij))]
                        }
                    }
              }
          }

        if (verbose)
            {
                jmessage('Finished correcting distances')
                print(Sys.time() -now)
            }
      }

    ## need to collapse matrix ie if there were "*" strand inputs and if we are counting internal
    ## connections inside our queries ..
    ## collapsing +/- rows and columns by max value based on their id mapping to their original "*" interval


    ## melt distance matrix into ij
    ij = as.matrix(expand.grid(1:nrow(D), 1:ncol(D)))
    dt = data.table(i = ij[,1], j = ij[,2], value = D[ij])[, id1 := gr1$id[i]][, id2 := gr2$id[j]]

    tmp = dcast.data.table(dt, id1 ~ id2, fun.aggregate = function(x) min(as.numeric(x)))
    setkey(tmp, id1)
    Dtmp = as.matrix(tmp[list(1:ngr1), -1, with = FALSE])
    D = matrix(NA, nrow = ngr1, ncol = ngr2, dimnames = list(NULL,
                                                             1:as.character(ngr2)))
    D[1:nrow(Dtmp), colnames(Dtmp)] = Dtmp


    ## finally zero out any intervals that actually intersect
    ## (edge case not captured when we just examine ends)
    if (length(intersect.ix)>0)
        D[cbind(intersect.ix$query.id, intersect.ix$subject.id)] = 0

    if (verbose)
        {
            jmessage('Finished aggregating distances to original object')
            print(Sys.time() -now)
        }
    return(D)
  }




#' jgraph
#'
#' takes in a jabba object and threshold for clusters and "quasi-reciprocal"
#' junctions
#'
#' @name jgraph
#' @export
#' @rdname internal
jgraph = function(jab, thresh_cl = 1e6, all = FALSE, thresh_r = 1e3, clusters = FALSE)
{
    ## identify sinks and sources
    jix = which(values(jab$junctions)$cn>0)
    so = gr.end(jab$segstats[as.integer(jab$ab.edges[jix,1,1:2])], ignore.strand = FALSE)
    si = gr.start(jab$segstats[as.integer(jab$ab.edges[jix,2,1:2])], ignore.strand = FALSE)
    jid = rep(jix, 2)
    recip = rep(c(FALSE, TRUE), each = length(jix))


    ## compute directed and undirected distances from sinks to sources
    gundir = jabba.dist(jab, si, so, directed = FALSE)
    gdir = jabba.dist(jab, si, so, directed = TRUE)


    gedges = gdir<=thresh_cl
    redges = gundir<=thresh_r & gdir>thresh_cl ## reciprocal edges are indirect that are non-indirect connections


    ## pick only closest g edge leaving a node
    ## this will be asymmetric
    ## since a_j may be closest to b_i
    ## but aa_i may not be closest to bb_j
    gfilt = t(apply(gdir, 1, function(x) sign((1:length(x) %in% which.min(x)))))

    ## key mathematical property:
    ## (identical to stranded adjacency matrix in jabba)
    ## recip, recip is identical to t(!recip, !recip)
    ## BUT recip, !recip and !recip, recip are distinct and symmetric

    ## just try these checks all should be empty
    ## which(gedges[!recip, !recip] != t(gedges[recip, recip]), arr.ind = TRUE)
    ## which(gedges[!recip, recip] != t(gedges[!recip, recip]), arr.ind = TRUE)
    ## which(gedges[recip, !recip] != t(gedges[recip, !recip]), arr.ind = TRUE)


    Gg.dist = pmin(gdir[recip, recip], gdir[!recip, recip], gdir[recip, !recip], gdir[!recip, recip])
    Gr.dist = pmin(gundir[recip, recip], gundir[!recip, recip], gundir[recip, !recip], gundir[!recip, recip])


    ## Gg and Gr should be symmetric
    Gg = sign(gedges[recip, recip]) + sign(gedges[!recip, recip]) + sign(gedges[recip, !recip]) + sign(gedges[!recip, !recip])
    Gr = sign(redges[recip, recip]) + sign(redges[recip, !recip]) + sign(redges[!recip, recip]) + sign(redges[!recip, !recip])
    Gr[cbind(1:nrow(Gr), 1:nrow(Gr))] = 0 ## no self reciprocal edges, as tempting as that might be

    ## Ggf won't be symmetric but we will make it and also make sure we don't remove self loops
    Ggf = sign(gfilt[recip, recip]) + sign(gfilt[!recip, recip]) + sign(gfilt[recip, !recip]) + sign(gfilt[!recip, !recip])
    Ggf = Ggf+t(Ggf) + Matrix::diag(rep(1, nrow(Ggf)))

    M = Matrix(0,
               nrow = length(jab$junctions),
               ncol = length(jab$junctions))

    out = list(Gg = M,
               Gr = M,
               Ggd = M,
               Grd = M)

    out$Gg[jix, jix] = sign(Gg)*sign(Ggf)
    out$Gr[jix, jix] = sign(Gr)

    ## some fake reciprocals leak through somehow - get rid!! (TOFIX)
                                        #    out$Gr[which(out$Gg!=0, arr.ind = TRUE)] = 0

    if (clusters)
    {
        tmp.out = tryCatch(split(1:nrow(out$Gr),
                             igraph::clusters(igraph::graph.adjacency(out$Gg + out$Gr))$membership), error = function(e) NULL)
        if (is.null(tmp.out))
            tmp.out = split(1:nrow(out$Gr),
                                 igraph::clusters(igraph::graph.adjacency(as(out$Gg + out$Gr, 'matrix')))$membership)
        out = tmp.out[rev(order(sapply(tmp.out, length)))]
        names(out) = 1:length(out)
    }

    return(out)
}





#' @name jab2json
#' @export
#' @rdname internal
#' @title jab2json
#'
#' @description
#'
#' Dumps JaBbA graph into json
#'
#'
#'
#' @param jab input jab object
#' @param file output json file
#' @author Marcin Imielinski
jab2json = function(jab, file, maxcn = 100, maxweight = 100)
{

    #' ++ = RL
    #' +- = RR
    #' -+ = LL
    qw = function(x) paste0('"', x, '"')

    ymin = 0;
    ymax = maxcn;

    nodes = jab$segstats %Q% (strand == "+")
    id = rep(1:length(nodes), 2)
    id.type = ifelse(nodes$loose, 'loose_end', 'interval')
    str = ifelse(as.character(strand(jab$segstats))=='+', 1, -1)

    node.dt = data.table(
        iid = 1:length(nodes),
        chromosome = qw(as.character(seqnames(nodes))),
        startPoint = as.character(start(nodes)),
        strand = "*",
        endPoint = as.character(end(nodes)),
        title = as.character(1:length(nodes)),
        type = ifelse(nodes$loose, "loose_end", "interval"),
        y = pmin(maxcn, nodes$cn))

    aadj = jab$adj*0
    rix = which(rowSums(is.na(jab$ab.edges[, 1:2, '+']))==0)
    aadj[rbind(jab$ab.edges[rix, 1:2, '+'], jab$ab.edges[rix, 1:2, '+'])] = 1
    ed = which(jab$adj!=0, arr.ind = TRUE)

    if (nrow(ed)>0)
        {
            ed.dt = data.table(
                so = id[ed[,1]],
                so.str = str[ed[,1]],
                si = id[ed[,2]],
                weight = jab$adj[ed],
                title = "",
                type = ifelse(aadj[ed], 'ALT', 'REF'),
                si.str = str[ed[,2]])[, sig := ifelse(so<si,
                                                      paste0(so * so.str, '_', -si*si.str),
                                                      paste0(-si * si.str, '_', so*so.str)
                                                      )][!duplicated(sig), ][, cid := 1:length(weight), ][,
                                                                                                          ":="(so = so*so.str, si = -si*si.str)]
            connections.json = ed.dt[, paste0(
                c("connections: [", paste(
                                        "\t{",
                                        "cid: ", cid,
                                        ", source: ", so,
                                        ", sink:", si,
                                        ", title: ", qw(title),
                                        ", type: ", qw(type),
                                        ", weight: ", pmin(maxweight, weight),
                                        "}",
                                        sep = "",
                                        collapse = ',\n'),
                  "]"),
                collapse = '\n')
                ]
        }

    intervals.json = node.dt[, paste0(
        c("intervals: [", paste(
                              "\t{",
                              "iid: ", iid,
                              ", chromosome: ", chromosome,
                              ", startPoint: ", startPoint,
                              ", endPoint: ", endPoint,
                              ", y: ", y,
                              ", title: ", qw(title),
                              ", type: ", qw(type),
                              ", strand: ", qw(strand),
                              "}",
                              sep = "",
                              collapse = ',\n'),
          "]"),
        collapse = '\n')
        ]

    meta.json =
        paste('meta: {\n\t',
              paste(
                  c(paste('"ymin:"', ymin),
                  paste('"ymax:"', ymax)),
                  collapse = ',\n\t'),
              '\n}')

    out = paste(c("var json = {",
                  paste(
                      c(meta.json,
                      intervals.json,
                      connections.json),
                      collapse = ',\n'
                  ),"}"),
                  sep = "")

    writeLines(out, file)
}



###############################################################
#' @name karyoMIP
#' @export
#' @rdname internal
#' karyoMIP
#'
#' MIP to locally compute walks in an existing JaBbA reconstruction, note: usually many optimal solutions to a given run.
#' Used by jabba.walk.
#'
#' TODO: Make user friendly, still pretty raw
#'
#' takes |E| x k matrix of k extreme paths (i.e. contigs) across e edges of the karyograph
#' and length |E| vector of edge copy numbers (eclass), length |E| vector of edge equivalence classes (both outputs of jbaMIP.process)
#' and computes most likely karyotypes that fit the edge copy number profile subject to some prior likelihood
#' over the k extreme paths
#'
#' @param K  |E| x k binary matrix of k "extreme" contigs across |E| edges
#' @param e  edge copy numbers across |E| edges
#' @param eclass  edge equivalence classes, used to constrain strand flipped contigs to appear in solutions together, each class can have at most 2 members
#' @param prior  prior log likelihood of a given contig being in the karyotype
#' @param cpenalty karyotype complexity penalty - log likelihood penalty given to having a novel contig in the karyotype, should be calibrated to prior, i.e. higher than the contig-contig variance in the prior, otherwise complex karyotypes may be favored
#' @param tilim time limit to optimizatoin
#' @param nsolutions how many equivalent solutions to report
#' @return
#' Rcplex solution list object with additional field $kcn for path copy number, $kclass for k class id, $mval for mval
###############################################################
karyoMIP = function(K, # |E| x k binary matrix of k "extreme" contigs across |E| edges
  e, # edge copy numbers across |E| edges
  eclass = 1:length(e), # edge equivalence classes, used to constrain strand flipped contigs to appear in solutions together,
                        # each class can have at most 2 members
  kclass = NULL,
  prior = rep(0, ncol(K)), # prior log likelihood of a given contig being in the karyotype
  mprior = NULL, # matrix prior which should be a binary matrix of m x k, eg mapping contigs to their read / barcode support
                 # will result in the addition of a quadratic objective in addition to the complexity penalty
  cpenalty = 1, # karyotype complexity penalty - log likelihood penalty given to having a novel contig in the karyotype,
                # should be calibrated to prior, i.e. higher than the contig-contig variance in the prior,
                # otherwise complex karyotypes may be favored
  tilim = 100, epgap = 1, nsolutions = 50, objsense = 'max', ...)
  {
    M = 1e7;
    K = as(K, 'sparseMatrix')

    if (length(prior)!=ncol(K))
      stop('prior must be of the same length as number of columns in K')

    # variable indices
    v.ix = 1:ncol(K)
    M.ix = max(v.ix) + (1:ncol(K))
    n = max(M.ix);

    # add big M constraints
    Zero = sparseMatrix(1, 1, x = 0, dims = c(n, n)) # upper bound is "infinity" if indicator is positive
    Amub = Zero[1:length(M.ix), ]
    Amub[cbind(1:length(M.ix), v.ix)] = 1
    Amub[cbind(1:length(M.ix), M.ix)] = -M

    Amlb = Zero[1:length(M.ix), ] # lower bound a little > 0 if indicator is positive
    Amlb[cbind(1:length(M.ix), v.ix)] = 1
    Amlb[cbind(1:length(M.ix), M.ix)] = -0.1

    if (is.null(kclass))
      kclass = .e2class(K, eclass)

    kclass.counts = table(kclass)
    if (any(kclass.counts>1)) ## any equiv i.e. strand flipped contig pairs? then make sure they appear in solutions togethrer
      {
        bikclass = which(kclass.counts>1)
        Ac = Zero[1:length(bikclass), ]
        pairs = matrix(unlist(split(1:length(kclass), kclass)[as.character(bikclass)]), ncol = 2, byrow = T)
        Ac[cbind(1:nrow(pairs), pairs[,1])] = 1
        Ac[cbind(1:nrow(pairs), pairs[,2])] = -1
      }
    else
      Ac = Zero[1,,drop = FALSE]

    # combine constraints
    A = rBind(cBind(K, Zero[rep(1, nrow(K)), M.ix]), Amub, Amlb, Ac);
    b = c(e, rep(0, nrow(Amlb)*2), rep(0, nrow(Ac)));
    sense = c(rep('E', nrow(K)), rep('L', nrow(Amlb)), rep('G', nrow(Amlb)), rep('E', nrow(Ac)))
    vtype = c(rep('I', length(v.ix)), rep('B', length(M.ix)))
    cvec = c(rep(0, length(v.ix)), prior-cpenalty*rep(1, length(M.ix)))

    if (is.null(mprior))
       sol = Rcplex(cvec = cvec, Amat = A, bvec = b, sense = sense, Qmat = NULL, lb = 0, ub = Inf, n = nsolutions, objsense = objsense, vtype = vtype, control = c(list(...), list(tilim = tilim, epgap = epgap)))
    else
    {
      if (!is.matrix(mprior))
        stop('mprior must be matrix')

      if (ncol(mprior) != ncol(K))
        stop('mprior must be matrix with as many columns as there are walks')

      m = nrow(mprior)
      jmessage('Adding mprior to karyoMIP')

      ## column cat matrices with blank Zero matrix on left, with
      ## constraints acting on binary variables, then identity matrix on right to capture
      ## the new "barcode residual" variables and their associated indicator variables

      ## goal is to maximize barcode usage
      ## the "residual" of each row of binary constraints
      ## by adding (linear objective) that is a weighted function of the
      ## number

      ## barcode constraints + 2*m additional variabes
      Ap = cbind(Zero[rep(1, nrow(mprior)), rep(1, length(M.ix))], sign(mprior), -diag(rep(1, nrow(mprior))), 0*diag(rep(1, nrow(mprior))))
      prix = n + 1:m ## indices of prior residuals
      iprix = n + m + 1:m ## indices of indicators of prior residuals
      pb = rep(0, nrow(mprior))
      psense = rep('E', nrow(mprior))

      Mplb = Mpub = 0*Ap ## upper and lower bounds for indictors
      Mpub[cbind(1:length(prix), prix)] = 1
      Mpub[cbind(1:length(prix), iprix)] = -M
      Mplb[cbind(1:length(prix), prix)] = 1
      Mplb[cbind(1:length(prix), iprix)] = -0.1
      pmb = rep(0, 2*nrow(Mpub))
      pmsense = c(rep('L', nrow(Mpub)), rep('G', nrow(Mplb)))

      ## define additional variables
      pvtype = c(rep('C', nrow(mprior)), rep('B', nrow(mprior)))

      ## objective function weighs the rows of mprior (barcodes) according to their max weight
      ## so user can weigh importance of individual barcodes
      ## or tune the overall importance of barcodes vs parsimony
      pcvec = c(rep(0, m), apply(mprior, 1, max))

      A = rBind(cBind(A, sparseMatrix(1, 1, x = 0, dims = c(nrow(A), 2*m))), Ap, Mpub, Mplb)
      b = c(b, pb, pmb)
      sense = c(sense, psense, pmsense)
      vtype = c(vtype, pvtype)
      cvec = c(cvec, pcvec)

      jmessage('Solving optimization with additional ', m, ' matrix prior terms')

      sol = Rcplex(cvec = cvec, Amat = A, bvec = b, sense = sense, Qmat = NULL, lb = 0, ub = Inf, n = nsolutions, objsense = objsense, vtype = vtype, control = c(list(...), list(tilim = tilim, epgap = epgap)))
    }

    if (!is.null(sol$xopt))
      sol = list(sol)

    sol = lapply(sol, function(x)
      {
        x$kcn = round(x$xopt[v.ix])
        x$kclass = kclass
        x$mval= round(x$xopt[M.ix])
        return(x)
      })

    return(sol)
  }


##############################################################
#' @name karyoMIP.to.path
#' @export
#' @rdname internal
#' karyoMIP.to.path
#'
#' for a karyoMIP solution and associated K matrix of n x e elementary paths  (input to karyoMIP), and v x e edge signed incidence matrix
#'
#'
#' @param sol solution to karyoMIP
#' @param K matrix of elementary paths (input to karyoMIP)
#' @param e nrow(K) x 2 edge matrix representing vertex pairs (i.e. edges to which K is referring to)
#' @param gr optional GRanges whose names are indexed by rownames of B
#' @param mc.cores integer number of cores
#' @param verbose flag
#' @return
#' A list with following items:
#' $path length k list of paths, cycles (each item i is vector denoting sequence of vertices in G )
#' $is.cycle length k logical vector whose component i denotes whether path i is cyclic
#' $cn  length k integer vector whose component i denotes copy number of contig i
#' $path.grl if path.grl == T
##############################################################
karyoMIP.to.path = function(sol, ## karyoMIP solutions, i.e. list with $kcn, $kclass (edges vectors)
  K, ## K matrix input to karyomip (edges x paths)
  e, ## nrow(K) x 2 edge matrix representing vertex pairs (i.e. edges to which K is referring to)
  gr = NULL, ## optional GRanges who names are indexed by <<rownames>> of B
  mc.cores = 1,
  verbose = T
  )
{
  contigs = which(sol$kcn!=0)
  c1 =  contigs[!duplicated(sol$kclass[contigs])]
  c2 = setdiff(contigs, c1)
  c2 = c2[match(sol$kclass[c1], sol$kclass[c2])]
  contigs = c1
  contigs2 = c2

  nm.gr = names(gr)
  names(gr) = NULL

  if (is.null(nm.gr))
    nm.gr  = 1:length(gr)

  if (any(duplicated(nm.gr)))
    nm.gr = 1:length(gr)

  if (!is.character(e))
    e = matrix(as.character(e), ncol = 2)

  out = list();

  i1 = which(!is.na(e[,1]))
  i2 = which(!is.na(e[,2]))
  B = sparseMatrix(as.numeric(c(e[i1,1], e[i2,2])),  c(i1,i2), x = c(rep(-1, length(i1)), rep(1, length(i1))))
  rownames(B) = 1:nrow(B)

  ## tells us whether the given contig is a cycle .. cycles represent any path lacking net flow in a
  ## non-slack vertex

  is.slack = rowSums(is.na(e))!=0

  out$is.cyc = Matrix::colSums(K[is.slack, contigs, drop = F])==0 & Matrix::colSums((B %*% K[, contigs, drop = F])!=0)==0
  out$cn = sol$kcn[contigs]
  out$kix = contigs;
  out$kix2 = contigs2;

  K = K[, contigs, drop = F]
  out$paths = mclapply(1:length(contigs),
    function(i)
    {
      if (verbose)
        cat('contig', i, 'of', length(contigs), '\n')

      k = K[, i]
      v.all = setdiff(as.vector(e[k!=0,]), NA)
##      v.all = rownames(B)[which(rowSums(abs(B) %*% k)>0)]  ## vertices associated with edges in path / cycle  k

      if (length(v.all)==1) ## this is a slack to slack path involving 1 node
        return(v.all)

      ## make subgraph corresponding to edges in this path / cycle
##       B.tmp = B[, which(!is.slack)[k[!is.slack]!=0], drop = F] ##
##       so = rownames(B.tmp)[apply(B.tmp, 2, function(x) which(x<0))]
##       si = rownames(B.tmp)[apply(B.tmp, 2, function(x) which(x>0))]
##       sG = graph(rbind(so, si))
##       sG = graph(rbind(so, si))

      tmp.e = e[k!=0, ,drop = F]
      tmp.e = tmp.e[rowSums(is.na(tmp.e))==0,,drop = F]
      sG = graph(t(tmp.e))

      if (out$is.cyc[i])
        {
          p.fwd = names(get.shortest.paths(sG, v.all[1], v.all[pmin(length(v.all), 2)])$vpath[[1]])
          p.bwd = names(get.shortest.paths(sG, v.all[pmin(length(v.all), 2)], v.all[1])$vpath[[1]])
          return(unique(unlist(c(p.fwd, p.bwd))))
        }
      else
        {
          io = as.numeric(B[, !is.slack, drop = F] %*% k[!is.slack])
          v.in = rownames(B)[io<0][1]
          v.out = rownames(B)[io>0][1]
          return(names(get.shortest.paths(sG, v.in, v.out)$vpath[[1]]))
        }
    }, mc.cores = mc.cores)

  if (!is.null(gr))
      {
      if (is.null(nm.gr))
        nm.gr = names(B)
      names(gr) = NULL
      out$grl = do.call('GRangesList', lapply(out$paths, function(x) gr[match(x, nm.gr), c()]))  ## match non-slack vertices
      names(out$grl) = paste('Contig ', out$kix, ' (CN = ', out$cn, ')', sep = '')
      values(out$grl)$is.cycle = out$is.cyc
    }

  return(out)
}


#' @name karyotrack
#' @export
#' @rdname internal
#' @title karyo track
#' @details
#'
#' Takes karyograph and outputs gTrack +/- highlighting of one or more paths defined as GRanges or GRangesList (for multiple paths)
#' Edges will only be highlighted when the exact interval pair corresponding to the edge is included in
#' the graph
#'
#' @param kag  output of karyograph
#' @param paths GRanges or GRangesList
#' @return gTrack of karyograph with particular nodes / edges colored with specified colors
############################################
karyotrack = function(kag, paths = NULL, col = 'red', pad = 0)
    {
        if (length(paths)==0)
            paths = NULL

        edge.ix = which(kag$adj!=0, arr.ind = T) ## collect aberrant edges

        ## convert to "simplified form"
        edges = data.frame(from = edge.ix[,1], to = edge.ix[,2])

        estr = paste(edge.ix[,1], edge.ix[,2])
        abestr = paste(kag$ab.edges[,1,1:2], kag$ab.edges[,2,1:2])

        if (nrow(edges)>0)
            {
                edges$type = 'reference'
                if (any(ix <- estr %in% abestr))
                    edges$type[ix] = 'aberrant'

                edges$col = ifelse(edges$type == 'aberrant', alpha('purple', 0.4), alpha('gray10', 0.4))
                edges$h = 1
                edges$lwd = ifelse(edges$type == 'aberrant', 2, 1)
                edges$lty = 1
                edges$cex.arrow = 0
                edges$v = 1
                edges$not.flat = edges$type == 'aberrant'
                edges$v[edges$type == 'aberrant'] = 2
                edges$h[edges$type == 'aberrant'] = 2
                edges$dangle.w = 0.5
            }

        pos.ix = which( as.logical( strand(kag$tile)=='+') )
        kag$tile$tile.id = match(gr.stripstrand(kag$tile), gr.stripstrand(kag$tile[pos.ix]))
        ss = kag$tile
        ss$col = alpha('gray', 0.2)
        ss$border = alpha('black', 0.5)
        ss$ywid = 0.8
        ss$bin = suppressWarnings(disjointBins(ss+1+pad, ignore.strand = FALSE))

        out = gTrack()
        if (!is.null(paths))
            {
                if (is(paths, 'GRanges'))
                    paths = split(paths, 1)

                if (is.null(values(paths)$col))
                    values(paths)$col = col

                paths.u = grl.unlist(paths)[, c('col', 'grl.ix')] %**% ss[, c()]

                ss$col[paths.u$subject.id] = paths.u$col

                edges = as.data.table(edges)

                edges[, end.from := ifelse(as.logical(strand(ss)[from]=='+'), end(ss)[from], start(ss)[from])]
                edges[, start.to := ifelse(as.logical(strand(ss)[to]=='+'), start(ss)[to], end(ss)[to])]

                paths.u = gr2dt(paths.u)[, ix := 1:length(seqnames)]
                paths.up = paths.u[ , list(ix.from = (ix)[-(length(ix))], ix.to = (ix)[-1], color = col[1]), by = grl.ix]
                paths.up[, from := paths.u$subject.id[ix.from]][, to := paths.u$subject.id[ix.to]]
                paths.up[, strand.from := paths.u$strand[ix.from]][, strand.to := paths.u$strand[ix.to]]
                paths.up[ , end.from := ifelse(strand.from == '+', paths.u$end[ix.from], paths.u$start[ix.from])]
                paths.up[ , start.to := ifelse(strand.to == '+', paths.u$start[ix.to], paths.u$end[ix.to])]
                setkeyv(edges, c('from', 'to', 'end.from', 'start.to'))
                setkeyv(paths.up, c('from', 'to', 'end.from', 'start.to'))

                edges[, col := alpha(col, 0.2)];
                edges[, path.match := 0]
                                        #      edges[paths.up, col := color] ## update colors if there is a match on from to, end and start
                edges[paths.up, lwd := lwd*3] ## update colors if there is a match on from to, end and start
                edges[paths.up, col := alpha(col, 0.8)] ## update colors if there is a match on from to, end and start
                edges[paths.up, path.match := 1] ## update colors if there is a match on from to, end and start
                edges = edges[order(path.match), ]
                out = gTrack(paths, draw.paths = TRUE, name = 'Paths')
            }

        out = c(gTrack(ss, y.field = 'bin', edges = edges, name ='Karyograph', angle = 0, labels.suppress = TRUE, yaxis = FALSE), out)

        return(out)
    }




#' @name junction.paths
#' @rdname internal
#' junction.paths
#'
#' Applies "pigeonhole principle" to enumerate all junction paths
#' in a karyograph that can be proven to have copy number greater than 0
#'
#' Takes as input adjacency matrix specifying junction copy numbers
#' and numeric vector specifying node copy numbers.
#'
#' @param cn length n vector of integer copy numbers
#' @param adj nxn matrix of junction copy numbers
#' @return  list#' 
#' with fields
#' $paths  list of n paths, each path i consisting of an n_i x 2 matrix specifying sequences of n_i junctions (each an ij node pair)
#' $mcn minimum copy number associated with path i
#' @export
#################################################
junction.paths = function(cn, adj)
  {
    ## preallocate, preallocate, preallocate
    ed = Matrix::which(adj!=0)
    NMAX = length(cn)*3 ## should be larger than the number of anticipated paths
    EMAX = 1000
    BOOSTER.ROW = 1e4
    BOOSTER.COL = 500
    paths = array(NA, dim = c(NMAX, EMAX))
    firstnode = lastnode = lastix = mcn = rep(NA, NMAX)
    numpaths = 0
    numrows = nrow(paths) ## yes this is a very slow query for giant arrays
    numcols = ncol(paths)
    torem = rep(F, NMAX)

    .sub2ind = function(dim, r, c) (c-1)*dim[1] + r

    if (nrow(adj) != ncol(adj) | nrow(adj)!=length(cn))
      stop('Adjacency matrix must be n x n and cn must be a length n vector')

    for (i in which(!is.na(cn)))
      {
        outgoing.nodes = Matrix:which(adj[i, ]>0)
        incoming.nodes = Matrix:which(adj[, i]>0)
        outgoing.cn = adj[i, outgoing.nodes]
        incoming.cn = adj[incoming.nodes, i]
        outgoing.edges = .sub2ind(dim(adj), i, outgoing.nodes)
        incoming.edges = .sub2ind(dim(adj), incoming.nodes, i)

        ## augment existing paths and adjust their minimum cn
        if (length(outgoing.edges)>0)
          {
            if (numpaths>0)
              {
                ix =  which(lastnode[1:numpaths] == i)
                if( length(ix)>0 )
                  for (j in ix)
                    {
                      new.lastix = rep(lastix[j]+1, length(outgoing.edges))
                      new.paths = paths[rep(j, length(outgoing.edges)), 1:new.lastix[1], drop = F]
                      new.paths[, new.lastix] = outgoing.edges
                      new.mcn =  mcn[j]-(cn[i]-outgoing.cn)
                      new.paths = new.paths[new.mcn>0,, drop = F]
                      new.lastix = new.lastix[new.mcn>0]
                      new.firstnode = rep(firstnode[j], sum(new.mcn>0))
                      new.lastnode = outgoing.nodes[new.mcn>0]
                      new.mcn = new.mcn[new.mcn>0]

                      if (any(is.na(new.lastix)))
                        browser()

                      if (length(new.mcn)>0)
                        {
                          if (any(mcn[j]==new.mcn)) ## only throw away an old path if the mcn of the new path did not decay
                            {
                              torem[j] = T
                              firstnode[j] = lastnode[j] = lastix[j] = mcn[j] = NA
                            }

                          new.ix = numpaths + (1:nrow(new.paths))
                          paths[new.ix, 1:new.lastix[1]] = new.paths[, 1:new.lastix[1]]
                          mcn[new.ix] = new.mcn
                          lastix[new.ix] = new.lastix
                          lastnode[new.ix] = new.lastnode
                          firstnode[new.ix] = new.firstnode
                          numpaths = numpaths + nrow(new.paths)
                        }
                    }
              }

            ## add brand new paths
            new.ix = numpaths + (1:length(outgoing.edges))
            paths[new.ix, 1] = outgoing.edges
            mcn[new.ix] = outgoing.cn
            firstnode[new.ix] = i
            lastnode[new.ix] = outgoing.nodes
            lastix[new.ix] = 1
            numpaths = numpaths + length(outgoing.edges)
          }

        ## augment existing paths and adjust their minimum cn
        if (length(incoming.edges)>0)
          {
            if (numpaths>0)
              {
                ix =  which(firstnode[1:numpaths] == i)
                if( length(ix)>0 )
                  for (j in ix)
                    {
                      now = Sys.time()
                      new.paths = cbind(0, paths[rep(j, length(incoming.edges)), 1:lastix[j], drop = F])
                      new.lastix = rep(lastix[j]+1, length(incoming.edges))
                      new.paths[, 1] = incoming.edges
                      new.mcn =  mcn[j]-(cn[i]-incoming.cn)
                      new.paths = new.paths[new.mcn>0,, drop = F]
                      new.lastix = new.lastix[new.mcn>0]
                      new.firstnode = incoming.nodes[new.mcn>0]
                      new.lastnode = rep(lastnode[j], sum(new.mcn>0))
                      new.mcn = new.mcn[new.mcn>0]

                      if (any(is.na(new.lastix)))
                        browser()

                      if (length(new.mcn)>0)
                        {
                          if (any(mcn[j]==new.mcn)) ## only throw away an old path if the mcn of a new path did not decay
                            {
                              torem[j] = T
                              firstnode[j] = lastnode[j] = lastix[j] = mcn[j] = NA
                            }

                          new.ix = numpaths + (1:nrow(new.paths))
                          paths[new.ix, 1:new.lastix[1]] = new.paths[, 1:new.lastix[1]]
                          mcn[new.ix] = new.mcn
                          lastix[new.ix] = new.lastix
                          lastnode[new.ix] = new.lastnode
                          firstnode[new.ix] = new.firstnode
                          numpaths = numpaths + nrow(new.paths)
                        }
                    }
              }

            ## add brand new paths
            new.ix = numpaths + (1:length(incoming.edges))
            paths[new.ix, 1] = incoming.edges
            mcn[new.ix] = incoming.cn
            firstnode[new.ix] = incoming.nodes
            lastnode[new.ix] = i
            lastix[new.ix] = 1
            numpaths = numpaths + length(incoming.edges)
          }

        ## augment if necessary
        if (numpaths>(numrows-BOOSTER.ROW))
          {
            cat('Allocating more row space\n')
            paths = rbind(paths, array(NA, dim = c(BOOSTER.ROW, ncol(paths))))
            numrows = numrows + BOOSTER.ROW
            torem = c(torem, rep(NA, BOOSTER.ROW))
            firstnode = c(firstnode, rep(NA, BOOSTER.ROW))
            lastnode = c(lastnode, rep(NA, BOOSTER.ROW))
            lastix = c(lastix, rep(NA, BOOSTER.ROW))
            mcn = c(mcn, rep(NA, BOOSTER.ROW))
          }

        if (max(c(0, lastix[1:numpaths]), na.rm = T)>(numcols-BOOSTER.COL))
          {
            cat('Allocating more column space\n')
            paths = cbind(paths, array(NA, dim = c(nrow(paths), BOOSTER.COL)))
            numcols = numcols + BOOSTER.COL
          }

        if ((i %% 500)==0)
          {
            cat(i, numpaths, '\n')
#            print(table(rowSums(!is.na(paths[1:numpaths,]))))
            cat('all last nodes\n')
            print(sort(table(seqnames(this.asol$asegstats[setdiff(lastnode, NA)]))))

            cat('all traversed nodes\n')
            print(sort(table(seqnames(this.asol$asegstats[1:i]))))
            rc = ind2sub(dim(adj), setdiff(as.vector(paths[1:numpaths, 1:max(lastix, na.rm = T)]), NA))

            cat('all path nodes\n')
            print(sort(table(seqnames(this.asol$asegstats[unique(as.numeric(rc))]))))

            cat('adj\n')
            print(sort(table(adj[rc])))

            keep.ix = which(!torem[1:numpaths])
            tmp.lastix = lastix[keep.ix]
            tmp.paths = paths[keep.ix, 1:max(tmp.lastix), drop = F]
            tmp.firstnode = firstnode[keep.ix]
            tmp.lastnode = lastnode[keep.ix]
            tmp.mcn = mcn[keep.ix]
            tmp.numpaths = length(keep.ix)

            paths[1:numpaths, 1:max(tmp.lastix)] = NA
            lastix[1:numpaths] = mcn[1:numpaths] = firstnode[1:numpaths] = lastnode[1:numpaths] = NA
            torem[1:numpaths] = F

            paths[1:tmp.numpaths, 1:max(tmp.lastix)] = tmp.paths[, 1:max(tmp.lastix), drop = F]
            firstnode[1:tmp.numpaths] = tmp.firstnode
            lastnode[1:tmp.numpaths] = tmp.lastnode
            mcn[1:tmp.numpaths] = tmp.mcn
            lastix[1:tmp.numpaths] = tmp.lastix
            numpaths = tmp.numpaths
 #           print(lapply(order(-rowSums(!is.na(paths[1:numpaths, ])))[1:2], function(x) paths[x, !is.na(paths[x, ])]))

            saveRDS(list(paths = paths[1:numpaths,], mcn = mcn[1:numpaths]), 'paths.rds')
            cat(i, numpaths, '\n')
          }
      }

    keep = which(rowSums(!is.na(paths)) !=0 & !torem)

    paths = paths[keep, 1:max(lastix, na.rm = T)]
    paths = lapply(1:nrow(paths), function(x) paths[x, !is.na(paths[x, ])])
    mcn = mcn[keep]

    return(list(paths = paths, mcn = mcn))
  }




#################################################
#' @name loose.ends
#' @rdname internal
#' loose.ends
#'
#' takes jbaMIP output and outputs a vector of ranges
#' on the right or left end of the intervals
#' that have type 1-4 labels where
#'
#' type1 = cn drop, no junction on this side, slack
#' type2 = cn drop, no junction on this side slack
#' type3 = no cn diff, used junction on other side, slack
#' type4 = no cn diff, unused junction on other side, no slack
#'
#'
#' @param sol JaBbA object
#' @param kag karyograph object
#' @return vector of ranges
#' on the right or left end of the intervals
#' that have type 1-4 labels where
#################################################
loose.ends = function(sol, kag)
  {
    if(!any(sol$segstats$eslack.in>0 | sol$segstats$eslack.out>0, na.rm = T))
      return(GRanges(seqinfo = seqinfo(kag$segstats)))

    nnab = !ifelse(is.na(sol$ab.edges[,3,1]), TRUE, sol$ab.edges[,3,1]==0)

    if (any(nnab))
        adj.ab = sparseMatrix(as.numeric(kag$ab.edges[nnab,1,]), as.numeric(kag$ab.edges[nnab,2,]),
            x = sol$adj[cbind(as.numeric(kag$ab.edges[nnab,1,]), as.numeric(kag$ab.edges[nnab,2,]))], dims = dim(sol$adj))
    else
        adj.ab = sol$adj*0

    ss = sol$segstats
    ss$num = 1:length(ss)
    n = length(ss)
    ss$left.ab = ss$cn.diff = ss$right.ab = -1;
    neg.ix = which( as.logical( strand(ss)=='-') )
    adj.ab[neg.ix, ] = adj.ab[rev(neg.ix), ]
    adj.ab[ ,neg.ix] = adj.ab[, rev(neg.ix)]
    ix.right = 1:n %in% as.numeric(kag$ab.edges[,1,])
    ix.left = 1:n %in% as.numeric(kag$ab.edges[,2,])
    ix.right[neg.ix] = ix.right[rev(neg.ix)]
    ix.left[neg.ix] = ix.left[rev(neg.ix)]

    ss[ as.logical( strand(ss)=='-' )] = rev(ss[ as.logical( strand(ss)=='-' ) ])
    tmp.right = Matrix::rowSums(adj.ab)
    tmp.left = Matrix::colSums(adj.ab)
    mask = c(as.numeric(diff(as.numeric(as.factor(seqnames(ss))))==0 & diff(as.numeric(as.factor(strand(ss))))==0), 0)
    ss$left.ab[ix.left] = tmp.left[ix.left]
    ss$right.ab[ix.right] = tmp.right[ix.right]
    ss$cn.diff = c(diff(ss$cn), 0)* mask

    ## now classify loose ends

    ## (1)
    ix.next = 1+(1:n)

    type1 = (rowSums(cbind(ss$eslack.out>0, ss$eslack.in[ix.next]>0))>0 &
             rowSums(cbind(ss$right.ab>0, ss$left.ab[ix.next]>0))==0 &
             ss$cn.diff != 0) * mask

    type2 = (rowSums(cbind(ss$eslack.out>0, ss$eslack.in[ix.next]>0))>0 &
             rowSums(cbind(ss$right.ab>0, ss$left.ab[ix.next]>0))>0 &
             ss$cn.diff != 0) * mask

    type3 = (rowSums(cbind(ss$eslack.out>0, ss$eslack.in[ix.next]>0))>0 &
             rowSums(cbind(ss$right.ab>0, ss$left.ab[ix.next]>0))==1 &
             ss$cn.diff == 0) * mask

    type4 = (rowSums(cbind(ss$right.ab==0, ss$left.ab[ix.next]==0))>0 &
             rowSums(cbind(ss$right.ab>0, ss$left.ab[ix.next]>0))==0 &
             ss$cn.diff == 0) * mask

    ss.p = ss[ as.logical( strand(ss)=='+' ) ]
    win.size = 1

    ss.p$num = 1:length(ss.p)

    slacks.tmp1 = gr.end(ss.p[which(ss.p$eslack.out>0)], win.size, force = T)
    if (length(slacks.tmp1)>0)
      {
        slacks.tmp1$type = '?'
        slacks.tmp1$type[which(ss.p$eslack.out>0) %in% which(type1!=0)] = 'type1'
        slacks.tmp1$type[which(ss.p$eslack.out>0) %in% which(type2!=0)] = 'type2'
        slacks.tmp1$type[which(ss.p$eslack.out>0) %in% which(type3!=0)] = 'type3'
        slacks.tmp1$cn = slacks.tmp1$eslack.out
        slacks.tmp1$sink = TRUE
      }

    slacks.tmp2 = gr.flipstrand(gr.start(ss.p[which(ss.p$eslack.in>0)], win.size, force = T))
    if (length(slacks.tmp2)>0)
      {
        slacks.tmp2$type = '?'
        slacks.tmp2$type[which(ss.p$eslack.in>0) %in% (which(type1!=0)+1)] = 'type1'
        slacks.tmp2$type[which(ss.p$eslack.in>0) %in% (which(type2!=0)+1)] = 'type2'
        slacks.tmp2$type[which(ss.p$eslack.in>0) %in% (which(type3!=0)+1)] = 'type3'
        slacks.tmp2$cn = slacks.tmp2$eslack.in
        slacks.tmp2$sink = FALSE ## i.e. source
      }

    ## with type 4 there is either a right facing ab w copy 0 or a left facing ab w copy 0 (from the "next" interval)
    t4.ix = intersect(which(type4!=0), 1:length(ss.p))
    t4.ix = ifelse(ss.p[t4.ix]$right.ab==0, -(t4.ix+1), t4.ix) ##

    slacks.t4 = c(gr.end(ss.p[t4.ix[t4.ix>0]], win.size, force = T), gr.flipstrand(gr.start(ss.p[-t4.ix[t4.ix<0]], win.size, force = T)))

    if (length(slacks.t4)>0)
      {
        slacks.t4$type = 'type4'
        slacks.t4$cn = 0
        slacks.t4$sink = t4.ix>0 ## t4.ix>0 means an incoming (i.e. left facing) ab edge w copy 0, hence the slack is a "sink" node leaving the other side of the breakpoint
      }

    loose.ends = grbind(slacks.tmp1, slacks.tmp2, slacks.t4)[, c('cn', 'num', 'type', 'sink')]

    return(loose.ends)
  }


#################################################
#' @name jbaMIP.allelic
#' @rdname internal
#' jbaMIP.allelic
#'
#' Takes adj and segstats from output from jbaMIP and
#' a granges of het.sites with $ref.count and $alt.count
#'
#' assumes segstats has fields $cn populated
#' and adj has copy states
#'
#' @param adj adjacency matrix populated with total copy counts on junctions
#' @param segstats granges tiling genome populated with total copy counts on interval, mu_high, sd_high, mu_low, sd_low variables on alleles
#' @param purity purity from solution
#' @param gamma gamma param from jbaMIP
#' @param slack.prior 1/slack.prior = penalty for each slack i.e. loose end copy in solution
#' @export
#' @return
#'
#' list with fields
#' $segstats = input segstats annotated with fitted cn.high, cn.low columns
#' $asegstats = output "allelic" segstats, with $cn, $parent.node, $eslack.in, $eslack.out, $phased fields filled in
#' $adj = output length(segstats) x length(segstats) x 2  "allelic" adjacency matrix with inferred allelic copy numbers on edges
#' $aadj = flattened output 2*length(segstats) x 2* length(segstats)  "allelic" adjacency matrix with inferred allelic copy numbers on edges
#'
#################################################
jbaMIP.allelic = function(
  adj, ## adjacency matrix populated with total copy counts on junctions
  segstats,  ## granges tiling genome populated with total copy counts on interval, mu_high, sd_high, mu_low, sd_low variables on alleles
  purity, ## purity from solution
  gamma, ## gamma param from jbaMIP
  partition = T,
  slack.prior = 0.001
  )
  {
    ploidy = sum(segstats$cn)/sum(as.numeric(width(segstats)))

    mu = c(segstats$mu_high, segstats$mu_low)
    ix = !is.na(mu)
    total = sum((mu * 2 * width(segstats))[ix])
    sw = sum(as.numeric(2*width(segstats))[ix])

    gamma = 2*(1-purity)/purity
    beta = (2*(1-purity)*sw + purity*ploidy*sw) / (purity * total)

    mu_high = segstats$mu_high*beta + gamma
    sd_high = segstats$sd_high*beta
    mu_low = segstats$mu_low*beta + gamma
    sd_low = segstats$sd_low*beta

    ## find the reference junctions
    ord.ix = order(segstats)
    rev.ix = as.logical(strand(segstats[ord.ix]) == '-')
    ord.ix = c(ord.ix[!rev.ix], rev(ord.ix[rev.ix]))

    ref.jun = cbind(ord.ix[-length(ord.ix)], ord.ix[-1])
    ref.jun = ref.jun[adj[ref.jun]>0, ]
    ab.adj = adj
    ab.adj[ref.jun] = 0
    ab.jun = Matrix::which(ab.adj!=0, arr.ind = T) ## aberrant junctions

    ## this will map vertices to their (positive) duplicate
    ## doing this will help contain some of the dimensionality
    dup.vmap = 1:length(segstats)
    pos.ix = which(as.logical( strand(segstats)=='+') )# "neg vertices" duplicates of pos vertices
    neg.ix = which( as.logical( strand(segstats)=='-') )
    dup.ix = suppressWarnings(neg.ix[gr.match(segstats[pos.ix], segstats[neg.ix])])
    dup.vmap[dup.ix] = pos.ix ## map neg vertices to their positive parent

    ## find dup ref and ab junctions
    ## these are neg-neg junctions (dup of pos-pos)
    ## (all neg-pos and pos-neg junctions are unique)
    dup.ref.emap = 1:nrow(ref.jun)
    ref.pos.ix = ref.jun[,1] %in% pos.ix & ref.jun[,2] %in% pos.ix
    ref.neg.ix = ref.jun[,1] %in% neg.ix & ref.jun[,2] %in% neg.ix
    ref.dup.ix = mmatch(ref.jun[ref.pos.ix, ], cbind(dup.vmap[ref.jun[ref.neg.ix,2]], dup.vmap[ref.jun[ref.neg.ix,1]])) ## ij ~ n(j)n(i)
    dup.ref.emap[ref.dup.ix] = ref.pos.ix

    dup.ab.emap = 1:nrow(ab.jun)
    ab.pos.ix = ab.jun[,1] %in% pos.ix & ab.jun[,2] %in% pos.ix
    ab.neg.ix = ab.jun[,1] %in% neg.ix & ab.jun[,2] %in% neg.ix
    ab.dup.ix = mmatch(ab.jun[ab.pos.ix, ], cbind(dup.vmap[ab.jun[ab.neg.ix,2]], dup.vmap[ab.jun[ab.neg.ix,1]])) ## ij ~ n(j)n(i)
    dup.ab.emap[ab.dup.ix] = ab.pos.ix

    ## we will need the following variables
    ## a1, a2 = length(segstats) allelic copy states for low (1) and high (2) state
    ## is1, is2 = incoming allele slack for low vs high state
    ## os1, os2 = outgoing allele slack for low vs high state
    ## r11, r12, r21, r22 = allelic copies on reference junctions for {low, high} x {low, high} combos
    ## n11, n12, n21, n22 = allelic copies on non-reference junctions for {low, high} x {low, high} combos
    ## i_r11, i_r12, i_r21, i_r22 = binary indicator variables representing positivity of the allelic copy state for ref junctions
    ## i_n11, i_n12, i_n21, i_n22 = binary indicator variables representing positivity of the allelic copy state for nonref junctions
    ## ns = (linearly penalized) non reference allelic junction slack (length non reference edges)
    ## eps1, eps2 = (quadratic penalized) epsilon residual between observed allelic segment means and the integer fit
    varmeta = data.frame( ## meta data of all variables we will be using
      var = c(
        rep(c('a1', 'a2'), each = length(segstats)),
        rep(c('is1', 'is2'), each = length(segstats)),
        rep(c('os1', 'os2'), each = length(segstats)),
        rep(c('r11', 'r12', 'r21', 'r22'), each = nrow(ref.jun)),
        rep(c('n11', 'n12', 'n21', 'n22'), each = nrow(ab.jun)),
        rep(c('i_r11', 'i_r12', 'i_r21', 'i_r22'), each = nrow(ref.jun)),
        rep(c('i_n11', 'i_n12', 'i_n21', 'i_n22'), each = nrow(ab.jun)),
        rep(c('i_is1', 'i_is2', 'i_os1', 'i_os2'), each = length(pos.ix)),
        rep('ns', nrow(ab.jun)),
        rep(c('eps1', 'eps2'), each = length(pos.ix))
        ),
      parent = c(
        rep(1:length(segstats), 2),
        rep(1:length(segstats), 2),
        rep(1:length(segstats), 2),
        rep(1:nrow(ref.jun), 4),
        rep(1:nrow(ab.jun), 4),
        rep(1:nrow(ref.jun), 4),
        rep(1:nrow(ab.jun), 4),
        rep(pos.ix, 4),
        1:nrow(ab.jun),
        rep(pos.ix, 2)
        ),
      stringsAsFactors = F)
    varmeta$vtype = 'I'
    varmeta$vtype[grepl('i_', varmeta$var)] = 'B'
    varmeta$vtype[grepl('eps', varmeta$var)] = 'C'
    varmeta$lb = 0
    varmeta$lb[varmeta$vtype == 'C'] = -Inf
    varmeta$ub = Inf
    varmeta$id = 1:nrow(varmeta)
    rownames(varmeta) = paste(varmeta$var, varmeta$parent)

    var = split(1:nrow(varmeta), varmeta$var) ## handy structure to keep track of variables
    Zero = sparseMatrix(1, 1, x = 0, dims = c(nrow(varmeta), nrow(varmeta)))

    ## we will have the following sets of constraints:
    ##
    ## copy1, copy2 = constraints linking counts to copy numbers via alpha, beta, and residual
    ## acopysum = constraints constraining allelic reference copy sums to equal total sums
    ## rcopysum, ncopysum = constraints constraining allelic junction sums to equal total sums
    ## iscopysum = constraints constraining incoming slack sums to equal total sums
    ## oscopysum = constraints constraining outoing slack sums to equal total sums
    ## rphase** = reference edge phase constraints that limit only one pair of r11, r12, r21, r22 to be positive
    ## nphase = aberrant edge phase constraints that limit only one allele pair to be non-negative
    ## oedge1, oedge2 = outgoing allelic edge conservation on each allele
    ## iedge1, iedge2 = incoming allelic edge conservation on each allele
    ## adup1, adup2 = duplicate constraints coupling positive interval copy to reverse complement interval copy
    ## rdup11, rdup12, rdup21, rdup22 = duplicate constraints coupling positive reference edge to reverse complement reference edge copy
    ## adup11, adup12, adup21, adup22 = duplicate constraints coupling positive aberrant edge copy to reverse complement aberrant edge copy
    ## isphase, osphase = phasing for incoming and outgoing slack
    ## i_r**, i_n**, i_is*, i_os* = "big M" constraints instantiating indicator variables
    consmeta = data.frame(
      cons = c(
        rep(c('copy1', 'copy2'), each = length(pos.ix)),
        rep(c('adup1', 'adup2'), each = length(pos.ix)),
        rep(c('acopysum'), each = length(pos.ix)),
        rep('rcopysum', each = length(ref.pos.ix)),
        rep('ncopysum', each = length(ab.pos.ix)),
        rep(c('rdup11', 'rdup12', 'rdup21', 'rdup22'), each = length(ref.pos.ix)),
        rep(c('ndup11', 'ndup12', 'ndup21', 'ndup22'), each = length(ab.pos.ix)),
        rep(c('oedge1', 'oedge2', 'iedge1', 'iedge2'), each = length(pos.ix)),
        rep(c('rphase*1', 'rphase*2', 'rphase1*', 'rphase2*'), each = length(ref.pos.ix)),
        rep(c('nphase'), each = length(ab.pos.ix)),
        rep(c('isphase', 'osphase'), each = length(pos.ix)),
        rep(c('Mlb_r11', 'Mlb_r12', 'Mlb_r21', 'Mlb_r22'), each = length(ref.pos.ix)),
        rep(c('Mub_r11', 'Mub_r12', 'Mub_r21', 'Mub_r22'), each = length(ref.pos.ix)),
        rep(c('Mlb_n11', 'Mlb_n12', 'Mlb_n21', 'Mlb_n22'), each = length(ab.pos.ix)),
        rep(c('Mub_n11', 'Mub_n12', 'Mub_n21', 'Mub_n22'), each = length(ab.pos.ix)),
        rep(c('Mlb_is1', 'Mlb_is2', 'Mlb_os1', 'Mlb_os2'), each = length(pos.ix)),
        rep(c('Mub_is1', 'Mub_is2', 'Mub_os1', 'Mub_os2'), each = length(pos.ix))
        ),
      num = c(
        rep(1:length(pos.ix), 2),
        rep(1:length(pos.ix), 2),
        rep(1:length(pos.ix), 1),
        rep(1:length(ref.pos.ix), 1),
        rep(1:length(ab.pos.ix), 1),
        rep(1:length(ref.pos.ix), 4),
        rep(1:length(ab.pos.ix), 4),
        rep(1:length(pos.ix), 4),
        rep(1:length(ref.pos.ix), 4),
        rep(1:length(ab.pos.ix), 1),
        rep(1:length(pos.ix), 2),
        rep(1:length(ref.pos.ix), 4),
        rep(1:length(ref.pos.ix), 4),
        rep(1:length(ab.pos.ix), 4),
        rep(1:length(ab.pos.ix), 4),
        rep(1:length(pos.ix), 4),
        rep(1:length(pos.ix), 4)
        ),
      stringsAsFactors = F
      )
    consmeta$sense = 'E'
    consmeta$sense[grepl('(phase)|(Mlb)|(Mub)', consmeta$cons)]= 'L'

    ## populate constraints
    ##
    n = nrow(varmeta)
    m = nrow(consmeta)
    A = sparseMatrix(1, 1, x = 0, dims = c(m, n)) #
    consmeta$b = rep(NA, length(n))

    M = 1e7;

    browser()

    ## copy state + eps constraints
    ## a1 = mu1 + eps
    cix = which(consmeta$cons == 'copy1')
    A[cbind(cix, varmeta[paste('a1', pos.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('eps1', pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = mu_high[varmeta[paste('eps1', pos.ix[consmeta$num[cix]]), 'parent']]

    ## a2 = mu2 + eps
    cix = which(consmeta$cons == 'copy2')
    A[cbind(cix, varmeta[paste('a2', pos.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('eps2', pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = mu_low[varmeta[paste('eps2', pos.ix[consmeta$num[cix]]), 'parent']]

    ## dup vertex constraints
    ## a1[neg.ix] = a1[dup.vmap[neg.ix]]
    cix = which(consmeta$cons == 'adup1')
    A[cbind(cix, varmeta[paste('a1', neg.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('a1', dup.vmap[neg.ix[consmeta$num[cix]]]), 'id'])] = -1
    consmeta$b[cix] = 0

    ## a2[neg.ix] = a2[dup.vmap[neg.ix]]
    cix = which(consmeta$cons == 'adup2')
    A[cbind(cix, varmeta[paste('a2', neg.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('a2', dup.vmap[neg.ix[consmeta$num[cix]]]), 'id'])] = -1
    consmeta$b[cix] = 0

    ## vertex allelic copy sum constraints
    ## a = a1 + a2
    cix = which(consmeta$cons == 'acopysum')
    A[cbind(cix, varmeta[paste('a1', pos.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('a2', pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = segstats$cn[varmeta[paste('a1', pos.ix[consmeta$num[cix]]), 'parent']]

    ## reference junction copy sum constraints
    ## r = r11 + r12 + r21 + r22
    cix = which(consmeta$cons == 'rcopysum')
    A[cbind(cix, varmeta[paste('r11', consmeta$num[cix]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('r12', consmeta$num[cix]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('r21', consmeta$num[cix]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('r22', consmeta$num[cix]), 'id'])] = 1
    consmeta$b[cix] = adj[ref.jun[ref.pos.ix[varmeta[paste('r22', consmeta$num[cix]), 'parent']], ]]

    ## aberrant junction copy sum constraints
    ## n = n11 + n12 + n21 + n22
    cix = which(consmeta$cons == 'ncopysum')
    A[cbind(cix, varmeta[paste('n11', consmeta$num[cix]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('n12', consmeta$num[cix]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('n21', consmeta$num[cix]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('n22', consmeta$num[cix]), 'id'])] = 1
    consmeta$b[cix] = adj[ab.jun[ab.pos.ix[varmeta[paste('n22', consmeta$num[cix]), 'parent']], ]]

    ## reference junction dup constraints
    ##
    cix = which(consmeta$cons == 'rdup11')
    A[cbind(cix, varmeta[paste('r11', ref.neg.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('r11', dup.ref.emap[ref.neg.ix[consmeta$num[cix]]]), 'id'])] = 1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'rdup12')
    A[cbind(cix, varmeta[paste('r12', ref.neg.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('r12', dup.ref.emap[ref.neg.ix[consmeta$num[cix]]]), 'id'])] = 1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'rdup21')
    A[cbind(cix, varmeta[paste('r21', ref.neg.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('r21', dup.ref.emap[ref.neg.ix[consmeta$num[cix]]]), 'id'])] = 1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'rdup22')
    A[cbind(cix, varmeta[paste('r22', ref.neg.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('r22', dup.ref.emap[ref.neg.ix[consmeta$num[cix]]]), 'id'])] = 1
    consmeta$b[cix] = 0

    ## aberrant junction dup constraints
    ##
    cix = which(consmeta$cons == 'rdup11')
    A[cbind(cix, varmeta[paste('r11', ab.neg.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('r11', dup.ab.emap[ab.neg.ix[consmeta$num[cix]]]), 'id'])] = 1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'rdup12')
    A[cbind(cix, varmeta[paste('r12', ab.neg.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('r12', dup.ab.emap[ab.neg.ix[consmeta$num[cix]]]), 'id'])] = 1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'rdup21')
    A[cbind(cix, varmeta[paste('r21', ab.neg.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('r21', dup.ab.emap[ab.neg.ix[consmeta$num[cix]]]), 'id'])] = 1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'rdup22')
    A[cbind(cix, varmeta[paste('r22', ab.neg.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('r22', dup.ab.emap[ab.neg.ix[consmeta$num[cix]]]), 'id'])] = 1
    consmeta$b[cix] = 0

    ## outgoing allelic edge constraints
    ##
    ## a1 = r11 + r12 + sum_k {n11}_k + sum_k {n12}_k + os1
    ##
    cix = which(consmeta$cons == 'oedge1')
    rj.ix = lapply(pos.ix[consmeta$num[cix]], function(x) which(ref.jun[,1] %in% x))
    aj.ix = lapply(pos.ix[consmeta$num[cix]], function(x) which(ab.jun[,1] %in% x))
    A[cbind(cix, varmeta[paste('a1', pos.ix[consmeta$num[cix]]), 'id'])] = -1
    for (i in 1:length(rj.ix))
      if (length(rj.ix[[i]])>0)
        {
          A[cbind(cix[i], varmeta[paste('r11', rj.ix[[i]]), 'id'])] = 1
          A[cbind(cix[i], varmeta[paste('r12', rj.ix[[i]]), 'id'])] = 1
        }

    for (i in length(aj.ix))
      if (length(aj.ix[[i]])>0)
        {
          A[cbind(cix[i], varmeta[paste('n11', aj.ix[[i]]), 'id'])] = 1
          A[cbind(cix[i], varmeta[paste('n12', aj.ix[[i]]), 'id'])] = 1
        }

    ## a2 = r21 + r22 + sum_k {n21}_k + sum_k {n22}_k + os2
    cix = which(consmeta$cons == 'oedge2')
    rj.ix = lapply(pos.ix[consmeta$num[cix]], function(x) which(ref.jun[,1] %in% x))
    aj.ix = lapply(pos.ix[consmeta$num[cix]], function(x) which(ab.jun[,1] %in% x))
    A[cbind(cix, varmeta[paste('a2', pos.ix[consmeta$num[cix]]), 'id'])] = -1
    for (i in 1:length(rj.ix))
      if (length(rj.ix[[i]])>0)
        {
          A[cbind(cix[i], varmeta[paste('r21', rj.ix[[i]]), 'id'])] = 1
          A[cbind(cix[i], varmeta[paste('r22', rj.ix[[i]]), 'id'])] = 1
        }
    for (i in length(aj.ix))
      if (length(aj.ix[[i]])>0)
        {
          A[cbind(cix[i], varmeta[paste('n21', aj.ix[[i]]), 'id'])] = 1
          A[cbind(cix[i], varmeta[paste('n22', aj.ix[[i]]), 'id'])] = 1
        }

    ## incoming allelic edge constraints
    ##
    ## a1 = r11 + r21 + sum_k {n11}_k + sum_k {n21}_k + is1
    ##
    cix = which(consmeta$cons == 'iedge1')
    rj.ix = lapply(pos.ix[consmeta$num[cix]], function(x) which(ref.jun[,2] %in% x))
    aj.ix = lapply(pos.ix[consmeta$num[cix]], function(x) which(ab.jun[,2] %in% x))
    A[cbind(cix, varmeta[paste('a1', pos.ix[consmeta$num[cix]]), 'id'])] = -1
    for (i in 1:length(rj.ix))
      if (length(rj.ix[[i]])>0)
        {
          A[cbind(cix[i], varmeta[paste('r11', rj.ix[[i]]), 'id'])] = 1
          A[cbind(cix[i], varmeta[paste('r21', rj.ix[[i]]), 'id'])] = 1
        }

    for (i in length(aj.ix))
      if (length(aj.ix[[i]])>0)
        {
          A[cbind(cix[i], varmeta[paste('n11', aj.ix[[i]]), 'id'])] = 1
          A[cbind(cix[i], varmeta[paste('n21', aj.ix[[i]]), 'id'])] = 1
        }

    ## a2 = r12 + r22 + sum_k {n12}_k + sum_k {n22}_k + is2
    cix = which(consmeta$cons == 'oedge2')
    rj.ix = lapply(pos.ix[consmeta$num[cix]], function(x) which(ref.jun[,2] %in% x))
    aj.ix = lapply(pos.ix[consmeta$num[cix]], function(x) which(ab.jun[,2] %in% x))
    A[cbind(cix, varmeta[paste('a2', pos.ix[consmeta$num[cix]]), 'id'])] = -1
    for (i in 1:length(rj.ix))
      if (length(rj.ix[[i]])>0)
        {
          A[cbind(cix[i], varmeta[paste('r12', rj.ix[[i]]), 'id'])] = 1
          A[cbind(cix[i], varmeta[paste('r22', rj.ix[[i]]), 'id'])] = 1
        }
    for (i in length(aj.ix))
      if (length(aj.ix[[i]])>0)
        {
          A[cbind(cix[i], varmeta[paste('n12', aj.ix[[i]]), 'id'])] = 1
          A[cbind(cix[i], varmeta[paste('n22', aj.ix[[i]]), 'id'])] = 1
        }

    ## reference phase constraints
    ##

    ## i_r11 + i_r12 <=1
    cix = which(consmeta$cons == 'rphase*1')
    A[cbind(cix, varmeta[paste('i_r11', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('i_r21', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 1

    ## i_r12 + i_r22 <=1
    cix = which(consmeta$cons == 'rphase*2')
    A[cbind(cix, varmeta[paste('i_r12', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('i_r22', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 1

    ## i_r11 + i_r12 <=1
    cix = which(consmeta$cons == 'rphase1*')
    A[cbind(cix, varmeta[paste('i_r11', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('i_r12', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 1

    ## i_r21 + i_r22 <=1
    cix = which(consmeta$cons == 'rphase2*')
    A[cbind(cix, varmeta[paste('i_r21', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('i_r22', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 1

    ## aberrant phase constraints
    ## i_n11 + i_n12 + i_n21 + i_n22 <= 1
    cix = which(consmeta$cons == 'nphase')
    A[cbind(cix, varmeta[paste('i_n11', ab.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('i_n12', ab.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('i_n21', ab.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('i_n22', ab.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 1

    ## incoming slack constraints
    cix = which(consmeta$cons == 'isphase')
    A[cbind(cix, varmeta[paste('i_is1', pos.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('i_is2', pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 1

    ## outgoing slack constraints
    cix = which(consmeta$cons == 'osphase')
    A[cbind(cix, varmeta[paste('i_os1', pos.ix[consmeta$num[cix]]), 'id'])] = 1
    A[cbind(cix, varmeta[paste('i_os2', pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 1

    ## "Big M" lower and upper bound constraints
    ##
    ## -M*i_r11 <= r11 <= M*i_r11
    cix = which(consmeta$cons == 'Mlb_r11')
    A[cbind(cix, varmeta[paste('i_r11', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('r11', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'Mub_r11')
    A[cbind(cix, varmeta[paste('i_r11', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('r11', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 0

    ## -M*i_r12 <= r12 <= M*i_r12
    cix = which(consmeta$cons == 'Mlb_r12')
    A[cbind(cix, varmeta[paste('i_r12', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('r12', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'Mub_r12')
    A[cbind(cix, varmeta[paste('i_r12', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('r12', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1

    ## -M*i_r21 <= r21 <= M*i_r21
    cix = which(consmeta$cons == 'Mlb_r21')
    A[cbind(cix, varmeta[paste('i_r21', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('r21', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'Mub_r21')
    A[cbind(cix, varmeta[paste('i_r21', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('r21', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1

    ## -M*i_r22 <= r22 <= M*i_r22
    cix = which(consmeta$cons == 'Mlb_r22')
    A[cbind(cix, varmeta[paste('i_r22', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('r22', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -1
    sconsmeta$b[cix] = 0

    cix = which(consmeta$cons == 'Mub_r22')
    A[cbind(cix, varmeta[paste('i_r22', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('r22', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1

    ## -M*i_n11 <= n11 <= M*i_n11
    cix = which(consmeta$cons == 'Mlb_n11')
    A[cbind(cix, varmeta[paste('i_n11', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('n11', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'Mub_n11')
    A[cbind(cix, varmeta[paste('i_n11', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('n11', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 0

    ## -M*i_n12 <= n12 <= M*i_n12
    cix = which(consmeta$cons == 'Mlb_n12')
    A[cbind(cix, varmeta[paste('i_n12', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('n12', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'Mub_n12')
    A[cbind(cix, varmeta[paste('i_n12', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('n12', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1

    ## -M*i_n21 <= n21 <= M*i_n21
    cix = which(consmeta$cons == 'Mlb_n21')
    A[cbind(cix, varmeta[paste('i_n21', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('n21', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'Mub_n21')
    A[cbind(cix, varmeta[paste('i_n21', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('n21', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1

    ## -M*i_n22 <= n22 <= M*i_n22
    cix = which(consmeta$cons == 'Mlb_n22')
    A[cbind(cix, varmeta[paste('i_n22', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('n22', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -1
    sconsmeta$b[cix] = 0

    cix = which(consmeta$cons == 'Mub_n22')
    A[cbind(cix, varmeta[paste('i_n22', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('n22', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1

    ## -M*i_is1 <= is1 <= M*i_is1
    cix = which(consmeta$cons == 'Mlb_is1')
    A[cbind(cix, varmeta[paste('i_is1', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('is1', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'Mub_is1')
    A[cbind(cix, varmeta[paste('i_is1', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('is1', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 0

    ## -M*i_is2 <= is2 <= M*i_is2
    cix = which(consmeta$cons == 'Mlb_is2')
    A[cbind(cix, varmeta[paste('i_is2', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('is2', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'Mub_is2')
    A[cbind(cix, varmeta[paste('i_is2', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('is2', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 0

    ## -M*i_os1 <= os1 <= M*i_os1
    cix = which(consmeta$cons == 'Mlb_os1')
    A[cbind(cix, varmeta[paste('i_os1', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('os1', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'Mub_os1')
    A[cbind(cix, varmeta[paste('i_os1', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('os1', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 0

    ## -M*i_os2 <= os2 <= M*i_os2
    cix = which(consmeta$cons == 'Mlb_os2')
    A[cbind(cix, varmeta[paste('i_os2', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('os2', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -1
    consmeta$b[cix] = 0

    cix = which(consmeta$cons == 'Mub_os2')
    A[cbind(cix, varmeta[paste('i_os2', ref.pos.ix[consmeta$num[cix]]), 'id'])] = -M
    A[cbind(cix, varmeta[paste('os2', ref.pos.ix[consmeta$num[cix]]), 'id'])] = 1
    consmeta$b[cix] = 0

  }


##############################################################
#' @name jbaMIP.summarize
#' @rdname internal
#' jbaMIP.summarize
#'
#' summarizes miqp result (i.e. multiple solutions) outputting data frame
#' with summary info
#'
##############################################################
jbaMIP.summarize = function(sol)
{
  if (!is.null(sol$obj))
    sol = list(sol)

  df = data.frame(purity = sapply(sol, function(x) x$purity),
    ploidy = sapply(sol, function(x) x$ploidy),
    obj = sapply(sol, function(x) x$obj),
    max.cn = sapply(sol, function(x) max(round(x$vcn), na.rm = T)),
    max.ecn = sapply(sol, function(x) max(round(x$ecn), na.rm = T)),
    tot.eslack = sapply(sol, function(x) c(sum(x$segstats$eslack.out) + sum(x$segstats$eslack.in), NA)[1]),
    num.eslack = sapply(sol, function(x) c(sum(x$segstats$eslack.out!=0) + sum(x$segstats$eslack.in!=0), NA)[1]),
    max.eslack = sapply(sol, function(x) pmax(max(c(x$segstats$eslack.out, 0)), max(c(x$segstats$eslack.in, 0)))))

  return(df)
}




#' @name ra_tier
#' @rdname internal
#' @title ra_tier
#' @description
#'
#' Classify full set of dRanger rearrangements into "tiers" of confidence
#'
#' (1) Tier 1 BPresult>0 and somatic_score>min.score1
#' (2) Tier 2 BPresult=0 and somatic_score>min.score1
#' (3) Tier 3 min.score2<=somatic_score<=min.score2 & tumreads>min.reads
#'
#'
##################
ra_tier = function(dra, min.score1 = 10, min.score2 = 4, min.treads1 = 10, min.treads2 = 3, max.nreads = Inf)
{
    if (is(dra, 'GRangesList')){
        dra = values(dra)
    }
    dra$BPresult[is.na(dra$BPresult)] = -1
    dra$T_SWreads[is.na(dra$T_SWreads)] = 0
    dra$N_SWreads[is.na(dra$N_SWreads)] = 0
    tier = rep(3, nrow(dra))
    tier[dra$BPresult==1 & dra$somatic_score>=min.score1 & (dra$tumreads >= min.treads1 | dra$T_SWreads > min.treads1 | dra$T_BWAreads > min.treads1) & dra$normreads==0 & dra$N_SWreads < 0] = 1
    tier[dra$BPresult>=0 & tier!=1 & dra$normreads == 0 & ((dra$tumreads + dra$T_SWreads) >= min.treads2 | (dra$somatic_score >= min.score2 | (dra$tumreads >= min.treads2)) )] = 2

    return(tier)
}




#' .e2class
#'
#' edge to contig class conversion
#'
#' given matrix K of k contigs over e edges, each belonging to cardinality 1 or cardinality 2 equivalence classes,
#' assigns id's to equivalent contigs
#'
####################################
.e2class = function(K, eclass)
{
  eclass = factor(as.character(eclass))

  if (length(eclass)!=nrow(K))
    stop('eclass must be of the same length as number of rows in K')

  eclass = factor(as.character(eclass))
  class.count = table(eclass);

  if (any(class.count)>2)
    stop('Edge equivalence classes can have at most 2 members')

  biclasses = names(class.count)[class.count==2];  # classes with two edges

  if (length(biclasses)>0)
  {
                                        # edge class rotation matrix
    R = diag(!(eclass %in% biclasses));  ## edges belonging to classes of cardinality 1 are on the diagonal

    ix = matrix(unlist(split(1:length(eclass), eclass)[biclasses]), ncol = 2, byrow = T); # index pairs corresponding to edges in biclasses
    R[ix[, 1:2]] = 1
    R[ix[, 2:1]] = 1

    Kr = R %*% K
    eix = mmatch(t(Kr), t(K))
    eix[is.na(eix)] = 0
    pairs = t(apply(cbind(1:length(eix), eix), 1, sort))
    pairs = pairs[!duplicated(pairs) & rowSums(pairs==0)==0, , drop = FALSE]

    kclass = rep(NA, ncol(K))
    kclass[pairs[,1]] = 1:nrow(pairs);
    kclass[pairs[,2]] = 1:nrow(pairs);
    kclass[is.na(kclass)] = nrow(pairs) + (1:sum(is.na(kclass)))
  }
  else
    kclass = 1:ncol(K)

  return(kclass)
}



#' @name abs2rel
#' @rdname internal
#' abs2rel
#'
#' rescales CN values from relative to "absolute" (i.e. per cancer cell copy) scale given purity and ploidy
#' By default, output is normalized to 1 (i.e. assumes that the total relative copy number signal mass over the genome is 1)
#'
#' takes in gr with signal field "field"
#' @param gr GRanges input with meta data field corresponding to mean relative copy "mean" in that interval
#' @param purity purity of sample
#' @param ploidy ploidy of sample
#' @param gamma gamma fit of solution (over-rides purity and ploidy)
#' @param beta beta fit of solution (over-rides purity and ploidy)
#' @param field meta character specifying meta data field in "gr" variable from which to extract signal, default "mean"
#' @param field.ncn character specifying meta data field in "gr" variable from which to extract germline integer copy number, default "ncn", if doesn't exist, germline copy number is assumed to be zero
#' @return
#' numeric vector of integer copy numbers
############################
abs2rel = function(gr, purity = NA, ploidy = NA, gamma = NA, beta = NA, field = 'cn', field.ncn = 'ncn', total = 1)
{
  abs = values(gr)[, field]
  w = width(gr)
  sw = sum(as.numeric(w))

  ncn = rep(2, length(mu))
  if (!is.null(field.ncn))
    if (field.ncn %in% names(values(gr)))
      ncn = values(gr)[, field.ncn]

  if (is.na(gamma))
    gamma = 2*(1-purity)/purity

  ploidy_normal = sum(w * ncn, na.rm = T) / sw  ## this will be = 2 if ncn is trivially 2

  if (is.na(beta))
    beta = ((1-purity)*ploidy_normal + purity*ploidy) * sw / (purity * total)
  ##  beta = (2*(1-purity)*sw + purity*ploidy*sw) / (purity * total)

                                        #    return((abs + gamma) / beta)
  return((abs + ncn*gamma/2) / beta)
}



#################################################
#' @name adj2inc
#' @rdname internal
#' adj2inc
#'
#' converts adjacency matrix (of directed graph) into incidence matrix - ie
#' an nodes x edges matrix, for each edge i connecting node j to k, column i will have -1 at position
#' j and +1 at position k
#'
#################################################
adj2inc = function(A)
{
  ij = which(A!=0, arr.ind = T)
  return(sparseMatrix(c(ij[,1], ij[,2]), rep(1:nrow(ij), 2), x = rep(c(-1, 1), each = nrow(ij)), dims = c(nrow(A), nrow(ij))))
}



######################################################
#' @name mmatch
#' @rdname internal
#' mmatch
#'
#' Low level utility function to match rows of matrix A to matrix B
#'
######################################################
mmatch = function(A, B, dir = 1)
{
  SEP = '|';

  if (is.null(dim(A)))
    A = rbind(A)

  if (is.null(dim(B)))
    B = rbind(B)

  if (dim(A)[(dir %% 2)+1] != dim(B)[(dir %% 2)+1])
    stop('Dimensions of A and B matrices mismatch')


  if (inherits(A, 'sparseMatrix') | inherits(B, 'sparseMatrix'))
  {
    if (dir == 2)
    {
      A = t(A)
      B = t(B)
    }
    ixA = Matrix::which(A!=0, arr.ind = T)
    ixB = Matrix::which(B!=0, arr.ind = T)

    if (nrow(ixA)>0)
      ixAl = split(1:nrow(ixA), ixA[,1])
    else
      ixAl = c()

    if (nrow(ixB)>0)
      ixBl = split(1:nrow(ixB), ixB[,1])
    else
      ixBl = c()

    Atxt = rep('', nrow(A))
    Btxt = rep('', nrow(B))

    if (length(ixAl))
    {
      tmp.ix = as.numeric(names(ixAl))
      Atxt[tmp.ix] = sapply(1:length(ixAl), function(x) paste(ixA[ixAl[[x]], 2], ':',
                                                              as.character(A[tmp.ix[x], ixA[ixAl[[x]], 2], drop = FALSE]), collapse = SEP))
    }

    if (length(ixBl)>0)
    {
      tmp.ix = as.numeric(names(ixBl))
      Btxt[tmp.ix] = sapply(1:length(ixBl), function(x) paste(ixB[ixBl[[x]], 2], ':',
                                                              as.character(B[tmp.ix[x], ixB[ixBl[[x]], 2], drop = FALSE]), collapse = SEP))
    }
  }
  else
  {
    Atxt = apply(A, dir, function(x) paste(x, collapse = SEP))
    Btxt = apply(B, dir, function(x) paste(x, collapse = SEP))
  }

  return(match(Atxt, Btxt))
}




####################################################################
#' @name jbaMIP.process
#' @rdname internal
#' jbaMIP.process
#'
#' process jbaMIP solution "sol" given original graph "g" (karyograph() list output)
#' into JaBbA object
#'
#' output is
#'
#' @param sol JaBbA object
#' @param allelic logical flag specifying whether object is allelic
#' @return
#' list with items:
#' $B incidence matrix of augmented graph (including slack vertices) (vertices x edges)
#' rownames of $B are vertex names of $G and colnames of B are named with character version of their $G indices
#' (i.e. column order of B  respects the original edge order in the solution)
#'
#' $e edge constraints for downstream karyoMIP, i.e the copy numbers at the edges
#' $e.ij numedges x 2 vertex pair matrix denoting what are the vertex pairs corresponding to the cols of $B and entries of $e, $eclass, $etype etc
#' $eclass id for each unique edge / anti-edge equivalence class
#' $etype specifies whether edge is slack or nonslack
###################################################################
jbaMIP.process = function(
                          ## output of jbaMIP, sol$segstats needs to have field $tile.id whose unique values appear exactly twice in the object,
                          ## corresponding to + and - strands of the same interval
                          sol,
                          allelic = F
                          )
{
  if (allelic)
    sol = list(segstats = sol$asegstats, adj = sol$aadj)

  if (!all(c('segstats', 'adj') %in% names(sol)))
    stop('sol must be output of jbaMIP()')

  if (is.null(sol$segstats$tile.id))
    stop('sol$segstats must be populated with tile.id')
  else
  {
    if (!all(table(sol$segstats$tile.id)==2))
      stop('sol$segstats$tile.id are malformed, there should be exactly two instances of each tile.id in sol$segstats, one for the positive and one for the negative strand of the same interval')

    tmp = lapply(split(1:length(sol$segstats$tile.id), sol$segstats$tile.id), rev)

    recip.ix = rep(NA, length(sol$segstats))
    recip.ix[order(sol$segstats$tile.id)] = unlist(tmp)
  }

  if (is.null(sol$segstats$eslack.in))
    sol$segstats$eslack.in = sol$segstats$slack.in

  if (is.null(sol$segstats$eslack.out))
    sol$segstats$eslack.out = sol$segstats$slack.out

  ed.ij = Matrix::which(sol$adj!=0, arr.ind = T)

  ## B is vertices x edges (i.e. signed incidence matrix)
  B = sparseMatrix(c(ed.ij[,1], ed.ij[,2]), rep(1:nrow(ed.ij), 2), x = rep(c(-1.00001, 1), each = nrow(ed.ij)), dims = c(nrow(sol$adj), nrow(ed.ij)))

  rownames(B) = 1:nrow(B)

  tmp.ix = Matrix::which(abs(B)>=1)
  B[tmp.ix] = round(B[tmp.ix]) ## "0.00001" hack to take care of eclass matching below, these are length 1 self loop edge cases

  ix.tel.5 = Matrix::which(Matrix::colSums(sol$adj!=0)==0)  ## make fake slacks for telomeres
  sol$segstats$eslack.in[ix.tel.5] = sol$segstats$cn[ix.tel.5]

  ix.tel.3 = Matrix::which(Matrix::rowSums(sol$adj!=0)==0)
  sol$segstats$eslack.out[ix.tel.3] = sol$segstats$cn[ix.tel.3]  ## make fake slacks for telomeres

  ix.eslack.out = Matrix::which(sol$segstats$eslack.out!=0);
  names(ix.eslack.out) = paste('out slack', ix.eslack.out)
  ix.eslack.in = which(sol$segstats$eslack.in!=0);
  names(ix.eslack.in) = paste('in slack', ix.eslack.in)

  names(ix.eslack.in)[ix.eslack.in %in% ix.tel.3] = paste(names(ix.eslack.in)[ix.eslack.in %in% ix.tel.3], 'tel')
  names(ix.eslack.out)[ix.eslack.out %in% ix.tel.5] = paste(names(ix.eslack.out)[ix.eslack.out %in% ix.tel.5], 'tel')

  ## we add "slack edges" and "slack nodes" to incidence matrix
  Zero = sparseMatrix(1, 1, x = 0, dims = c(length(ix.eslack.in) + length(ix.eslack.out), ncol(B)))

  if (nrow(Zero)>0)
    rownames(Zero) = c(paste('slack in', 1:length(ix.eslack.in)), paste('slack out', 1:length(ix.eslack.out)))

  Bs = rBind(B, Zero)
  ed.ij = rbind(ed.ij, cbind(ix.eslack.out, NA), cbind(NA, ix.eslack.in))

  Is = Diagonal(n = nrow(Bs), rep(1, nrow(Bs)))

  Bs = cBind(Bs, -Is[, ix.eslack.out], Is[, ix.eslack.in])
  colnames(Bs) = c(as.character(1:ncol(B)), names(ix.eslack.out), names(ix.eslack.in))

  ## map new "slack nodes" to their reciprocals
  recip.ix = c(recip.ix,
               nrow(B) + length(ix.eslack.out) +  match(recip.ix[ix.eslack.out], ix.eslack.in),
               nrow(B) + match(recip.ix[ix.eslack.in], ix.eslack.out)
               )

  ## match matrix against its reverse complement (i.e. rotation) to find reciprocal edges
  erecip.ix = mmatch(t(Bs), t(-Bs[recip.ix, ])) ## maps edges to their reciprocals

  tmp.na = which(is.na(erecip.ix))
  if (length(tmp.na)>0) ## fix the self loops so that they match
    erecip.ix[tmp.na] = tmp.na[mmatch(t(Bs[1:nrow(Bs), tmp.na]), t(Bs[recip.ix,tmp.na, drop = F]))]

  ## now use this mapping to define edge equivalence classes
  rmat = t(apply(cbind(erecip.ix, erecip.ix[erecip.ix]), 1, sort)) ## length(erecip.ix) x 2 matrix of edge ids and their reciprocal, sorted

  ## eclass will map length(erecip.ix) edges to length(erecip.ix)/2 edge equivalence class ids
  eclass = mmatch(rmat, rmat[!duplicated(rmat), ])

  Bs = round(Bs) ## remove the 0.0001 dummy coefficients (i.e. for self loops)

  ## e will store observed copy states corresponding to edges (i.e. columns of Bs)
  e = c(sol$adj[which(sol$adj!=0)], sol$segstats$eslack.out[ix.eslack.out],  sol$segstats$eslack.in[ix.eslack.in])

  return(list(e = e, e.ij = ed.ij, B = Bs, eclass = eclass, etype = c(ifelse(grepl('slack', colnames(Bs)), 'slack', 'nonslack'))))
}


#' @name spmelt
#' @title sparse matrix melt
#'
#' @description
#' Melts sparse matrix into data.table
#' 
#' @param A 
#' @return data.table of all non
#' @author Marcin Imielinski
spmelt = function(A, baseval = 0) {
  if (!is.null(baseval))
  {
    ij = Matrix::which(A!=baseval, arr.ind = TRUE)
  }
  else ## take everything
  {
    ij = as.matrix(expand.grid(1:nrow(A), 1:ncol(A)))
  }
  dt = data.table(i = ij[,1], j = ij[,2], val = A[ij])
}


#' @name dtt
#' @title dtt
#'
#' @description
#' Easy to type shortcut for setDTthreads
#'
#' @export
#' @param numthreads number of threads to set setDTthreads
#' @author Marcin Imielinski
dtt = function(threads = 1) setDTthreads(threads)


#' @name ggjs
#' @title ggjs
#' @description
#'
#' Deploys a gGnome.js by downloaded the js / node source code from mskilab.com
#' and then dumping json / csv files corresponding to the graph and coverage
#' data.
#'
#' Easiest way is to use pairs= argument which is a data.table with columns $pair, $jabba_rds, $cov_rds, $headline, $description
#' where $headline is a character vector of the headline that will be displayed in the search bar and $description is the information that will be shown when the graph is loaded. 
#' 
#' Takes gGraph objects or lists whose first item (or $graph) item is
#' a gGraph, second optional item (or $cov) is a GRanges of coverage
#' (where coverage data is specified by argument "field")
#'
#' If arguments are named then these will be the name of the resulting object
#' in the deployed gGnome.js app, otherwise they will be given default names
#' (ie graph1, graph 2).
#' 
#' Dumps the app to a directory where to which user can navigate and deploy the app
#' (instructions provided at runtime), default path is public_html/ggjs
#'
#' Note: coverage should be less than 1e6 bins per file, if binsize is not NULL
#' (default = 5000) then coverage data will be aggregated prior to dumping. 
#' 
#' example argument:
#' ggjs(path = gg, graph2 = list(gg2, cov2), list(cov = cov3, graph = gg3))
#' 
#' @param ... gGnome objects of list
#' @param field field of coverage GRanges that will be used to dump out the coverage data, default "ratio"
#' @param binsize binsize to aggregate coverage (default 5000)
#' @param skip.download skip download (useful if dumping files to existing directory)
#' @param web whether to dump web version which will update datafiles.csv and dump files to json and coverage subfolders of path (FALSE)
#' @param pairs data.table with $pair, $jabba_rds, $cov_rds field to dump, can have optional fields $headline and $description
#' @param path  path to dump app to, this will also be the command that will be run to deloy the app  ~/public_html/ggjs
#' @export
ggjs = function(...,
                pairs = NULL,
                field = 'ratio', path = '~/public_html/ggjs',
                binwidth = 5e3,
                mc.cores = 1,
                web = TRUE,
                skip.dl = FALSE,
                clean.up = TRUE, 
                force = FALSE,
                win = NULL,
                ggjs.url = 'http://mskilab.com/gGnome.js/ggjs.tar.gz',
                verbose = TRUE)
{
  if (verbose)
    message('Deploying gGnome.js app to ', path)


  if (!web)
  {
    if (file.info(path)$isdir)
      path = paste0(path, '/ggjs')

    ggjs.dir = paste0(path, '_files')
    system(paste('mkdir -p', ggjs.dir))

    if (!skip.dl)
    {
      if (verbose)
        message('Downloading and decompressing app code from ', ggjs.url)
      
      system(paste0('cd ', ggjs.dir, '; wget -q ', ggjs.url, '; tar xfz ', basename(ggjs.url)))
    }

    writeLines(sprintf('open http://localhost:8080/index.html; cd %s; npm start', basename(ggjs.dir)), path)
    system(paste('chmod +x', path))
    }
  else
  {

    if (is.null(pairs))
      stop('Publishing to web server only supported with pairs table argument')

    if (!dir.exists(path))
      {
        warning(sprintf('Creating fresh directory %s for web deployment', path))
        system(paste('mkdir -p', path))
      }
    ##      stop(sprintf('For publishing to web server, path must be an existing directory with gGnome.js app code deployed: please check path %s', path))
    
    
    ggjs.dir = path
    if (!file.exists(paste(path, 'js/server.js', sep = '/')))
    {
      warning(sprintf('Downloading app code from %s', ggjs.url))
      system(paste0('cd ', ggjs.dir, '; wget -q ', ggjs.url, '; tar xfz ', basename(ggjs.url)))
      ##        stop('For publishing to web server, path must exist and already have gGnome.js app code deployed: please check dir contents or reclone into path from https://github.com/mskilab/gGnome.js')
    }

    if (!(all(c('pair', 'jabba_rds', 'cov_rds') %in% names(pairs))))
    {
      stop('pairs must contain fields $pair, $jabba_rds, and $cov_rds')
    }

    pairs$description = paste0('', pairs$description)
    if (is.null(pairs$headline))
      pairs$headline = pairs$description

    datafiles = pairs[, .(datafile = paste0(pair, '.json'), description = paste0(headline))]

    datafiles.csv = paste(normalizePath(path), 'datafiles.csv', sep = '/')
    if (!file.exists(datafiles.csv))
    {
      warning('datafiles.csv missing from web server directory: creating')
    }
    else
    {
      if (verbose)
        message('reading datafiles.csv from gGnome.js web server deploy path ', datafiles.csv)

      datafiles = unique(rbind(datafiles, fread(datafiles.csv, sep = ',')), by = "datafile")      
    }

    fwrite(datafiles, datafiles.csv)
  }

  system(sprintf('mkdir -p %s/coverage %s/json', ggjs.dir, ggjs.dir))
  
  args = list(...)
  if (length(args)>0)
    {
      if (is.null(names(args)))
        names(args) = 1:length(args)
    }

  if (!is.null(pairs))
  {
    if (!(all(c('pair', 'jabba_rds', 'cov_rds') %in% names(pairs))))
    {
     stop('pairs must contain fields $pair, $jabba_rds, and $cov_rds')
    }

    ## optional argument
    pairs$description = paste0('', pairs$description)      

    setkey(pairs, pair)
    new.args = lapply(pairs$pair, function(p) as.list(pairs[.(p), list(id = pair, graph = jabba_rds, cov = cov_rds, description = description)]))
    names(new.args) = pairs$pair
    args = c(args, new.args)    
  }

  if (length(args)==0)
    stop('Either arguments or pairs data.table must be provided')
      
  if (any(nchar(names(args))==0))
    names(args)[nchar(names(args))==0] = "graph"
  
  if (any(duplicated(names(args))))
    names(args) = dedup(names(args))

  mclapply(names(args), function(nm)
  {   
    if (verbose)
      message('Dumping graph ', nm)

    if (!is.list(args[[nm]]))
      arg = list(args[[nm]])
    else
      arg = args[[nm]]

    if (is.null(names(arg)))
      names(arg) = c('graph', 'cov', 'walks')[1:length(arg)]

    if (is.character(arg$graph))
      {
        if (!is.na(arg$graph) && file.exists(arg$graph))
          arg$graph = readRDS(arg$graph)
        else
        {
          warning(sprintf('Provided graph file for sample %s is NA or does not exist', nm))
          arg$graph = NULL
        }
      }

    if (!is.null(arg$graph) && !inherits(arg$graph, 'gGraph'))
    {
      if (!inherits(arg$graph, 'gGraph'))
        arg$graph = suppressWarnings(gGnome::gG(jab = arg$graph))
      else
        arg$graph = gGnome:::refresh(arg$graph)

      arg$graph$set(description = paste0('', arg$description), name = paste0('', arg$id))
    }

    if (clean.up)
    {
      gg = arg$graph
      chrs = c(1:22, 'X', 'Y', paste0('chr', c(1:22, 'X', 'Y')))
      gg = gg[seqnames %in% chrs, ]
      gg = gG(nodes = gg$nodes$gr %>% gr2dt %>% cc(seqnames := gsub('chr', '', seqnames)) %>% dt2gr,
              edges = gg$edges$dt)
      gg$set(y.field = 'cn')
      arg$graph = gg 
    }
    
    graph.json = paste0(ggjs.dir, '/json/', nm, '.json')
    if (!is.null(arg$graph))
    {
      if (!is.null(win) & inherits(win, 'GRanges'))
        arg$graph = arg$graph[arg$graph$nodes$gr %^% win]
      sl = as.character(unique(seqnames(arg$graph$nodes$gr)))
      arg$graph$json(graph.json, seqlevels = sl)
    }

    cov.csv = paste0(ggjs.dir, '/scatterPlot/', nm, '.csv')
    if (is.character(arg$cov))
    {
      if (!is.na(arg$cov) && file.exists(arg$cov) && (force | !file.exists(cov.csv)))
        {
          arg$cov = readRDS(arg$cov)
          if (!is.null(win) & inherits(win, 'GRanges'))
            arg$cov = arg$cov[arg$cov %^% win]

          if (clean.up)            
            arg$cov = arg$cov %Q% (seqnames %in% c(1:22, 'X', 'Y', paste0('chr', c(1:22, 'X', 'Y')))) %>% gr2dt %>% cc(seqnames := gsub('chr', '', seqnames)) %>% dt2gr
        }
      else
        {
          arg$cov = NULL
        }
    }

    if (!is.null(arg$cov))
    {
      if (!inherits(arg$cov, 'GRanges'))
        stop('Coverage must be GRanges')

      if (!(field %in% names(values(arg$cov))))
        stop(paste0('field ', field, ' not found in provided coverage track'))

      if (verbose)
        message('Dumping coverage ', field, ' to binwidth of ', binwidth)

      ## trimming ranges to seqlengths of gg$gr
      arg$cov = arg$cov %&% si2gr(seqinfo(arg$graph))
      cov.dt = as.data.table(arg$cov[, field])
      setnames(cov.dt, field, 'coverage')

      if (!is.null(binwidth))
      {
        cov.dt = cov.dt[, .(coverage = sum(coverage*width, na.rm = TRUE)/sum(as.numeric(width)*(0*coverage+1), na.rm = TRUE)), by = .(seqnames, start = floor(start/binwidth)*binwidth+1)]
        cov.dt = cov.dt[!is.na(coverage), ][!is.infinite(coverage), ]
        if (verbose)
          message('Aggregating coverage via field ', field, ' to binwidth of ', binwidth, ' (from ', length(arg$cov), ' to ', nrow(cov.dt), ' bins)')
        
      }
      fwrite(cov.dt[, .(x= start, y = coverage, chromosome = seqnames)], cov.csv)
    }
  }, mc.cores = mc.cores, mc.preschedule = FALSE)

  if (verbose)
  {
    if (!web)
      message('App deployed! Navigate via terminal to ', ggjs.dir, ' and type "./', basename(path), '" in MacOS command line to launch.')
    else
      message('App deployed! Navigate via web to URL pointing to directory ', ggjs.dir)
  }
}
       
#' @name gr.eval
#' @title gr.eval
#' @description
#'
#' Evaluate a scalar expression on the metadata of a subject, returning
#' a vector of the same length as query with the values populated for
#' matching value (or fill value for non matching)
#'
#' @param query GRanges of intervals that we are interested in populating
#' @param subject GRanges of interval that has metadata that we want to aggregate / query
#' @param expr expression on subject metadata that we want to evaluate
#' @param fill flll value for empty / non overlapping intervas
#' @return vector of length query with expr evaluated on each overlapping interval
#' @export
gr.eval = function(query, subject, expr, fill = NA, ignore.strand = TRUE)
{
  if (ignore.strand)
    ov = query[, c()] %*% subject[, c()]
  else
    ov = query[, c()] %**% subject[, c()]

  ovdt = as.data.table(ov)
  ovdt = cbind(ovdt[, .(query.id, subject.id)], as.data.table(subject)[ovdt$subject.id, ])

  cmd = sprintf('ovdt[, .(V1 = %s, subject.id = subject.id[1]), keyby = query.id]', deparse(substitute(expr)))
  res = eval(parse(text = cmd))[.(1:length(query)), ]
  res[is.na(subject.id), V1 := fill]
  return(res$V1)  
}



#' @name rebin
#' @rdname internal
#' @title reaggregate WGS bins around a new target value
#' @description 
#' 
#' Given GRanges of bins will aggregate around new bin width. 
#' 
#' @param cov GRanges of binned genome-wide coverage
#' @param binwidth new binwidth
#' @return GRanges of binned genome-wide coverage at new bin
#' @export
#' @author Marcin Imielinski
rebin = function(cov, binwidth, field = names(values(cov))[1], FUN = mean, na.rm = TRUE)
{
  tmp = as.data.table(cov[, c()])
  tmp$value =  values(cov)[[field]]
  outdt = tmp[
    , FUN(value, na.rm = na.rm),
      by = .(seqnames, start = floor(start/binwidth)*binwidth + 1)]
  ## xtYao update:
  ## by = .(seqnames, start = ceiling(start/binwidth)*binwidth)]
  outdt[, end := start + binwidth -1]
  out = dt2gr(outdt)
  names(values(out)) = field
  return(out)
}


#' @name ssegment
#' @title wrapper around DNAcopy to segment a coverage profile
#' @description 
#' 
#' Internal function utilizing DNAcopy to segment a coverage profile
#' @param tcov GRanges of binned genome-wide coverage
#' @return GRanges of genomewise segments of piecewise constant coverage
#' @export
#' @author Marcin Imielinski
ssegment = function(cov, field = NULL, log = TRUE, verbose = TRUE, alpha = 1e-5){
  if (is.null(field))
    field = names(values(cov))[1] 
  cov$y = values(cov)[[field]]
  new.sl = seqlengths(cov)
  ix = which(!is.na(cov$y))
  if (verbose)
    message('sending ', length(ix), ' segments to DNAcopy')
  cov = cov[ix]
  if (log == TRUE)
    {
      cov$y = log(cov$y)
    }
  cna = DNAcopy::CNA(cov$y, as.character(seqnames(cov)), start(cov), data.type = 'logratio')
  gc()
  seg = DNAcopy::segment(DNAcopy::smooth.CNAsmooth.CNA(cna), alpha = alpha, verbose = 0)
  out = seg2gr(seg$out, new.sl) ## remove seqlengths that have not been segmented
  out = gr.fix(out, new.sl, drop = T)
  if (verbose)
    message('\t ..finished segmentation')
  return(out)
}


#' @name jhom
#' @title wrapper around GxG homeology for junctions
#' @description 
#' 
#' 
#' @param event Junction object around which to compute homeology
#' @param pad padding to put around breakpoint window
#' @param thresh maximum string distance threshold for calling a bin homeologous
#' @param pad2 pad to put around an output bin in order to measure homeology
#' @param flip flag whether to measure homeology with bins and their reverse complements
#' @param stats flag whether to return stats or gMatrix (default)
#' @param mc.cores  integer number of cores [1]
#' @param anchor logical flag whether to return homeology results around native (inputted) coordinates or Anchored coordinates [TRUE]
#' @param mat logical whether to return results as a mat [FALSE] or gMatrix
#' @return list of gMatrix objects (if stats = FALSE) and data.table of homeology stats (if stats = TRUE)
#' @export
#' @author Marcin Imielinski
jhom = function(event, pad = 100, thresh = 2, stride = 1, pad2 = 5, flip = FALSE, stats = FALSE, mc.cores = 1, anchor = TRUE, mat = FALSE)
{
  event = data.table(bp1 = gr.string(event$left),
                     bp2 = gr.string(event$right))

  gmfeat = function(gm, thresh = 3, op = "<=")
  {
    library(imager)
    if (is(gm, 'gMatrix'))
    {
      mat = gm$mat
      mat = (mat+t(mat))/2
    }
    else
      mat = gm
    im = mat %>% as.matrix %>% as.cimg
    cmd = sprintf("im %s %s", op, thresh)
    px = eval(parse(text = cmd))  
    if (sum(px)==0)
      return(data.table())
    sp = split_connected(px)
    if (length(sp)==0)
      return(data.table())
    res = rbindlist(lapply(1:length(sp), function(k) as.data.table(which(as.matrix(sp[[k]]), arr.ind = TRUE))[, .(k= k, i = pmin(row,col), j = pmax(row,col))]), fill = TRUE)
    return(res)
  }

  gmstats = function(res)
  {
    if (nrow(res)>0)
    {
      res[, .(N = .N, r = cor(i, j)), by = k][, .(
           numfeat = sum(N>0),
           numfeat2 = sum(N>2),
           numfeat5 = sum(N>5),
           numfeat10 = sum(N>10),
           maxfeat = max(N),
           numlines5 = sum(r>0.5, na.rm = TRUE),
           maxlines5 =  max(c(0,N[r>0.5]), na.rm = TRUE),
           numlines = sum(r>0.9, na.rm = TRUE),
           maxlines =  max(c(0,N[r>0.9]), na.rm = TRUE),
           maxcor = max(r, na.rm = TRUE)
         )]
    }
    else
      data.table(numfeat = 0, maxfeat = 0)
  }

  win = c(GRanges('Left:0')+pad,
          GRanges('Right:0')+pad)
  evbp1 = gr.end(gr.flipstrand(parse.gr(event$bp1)))
  if (flip)
    evbp2 = gr.flipstrand(gr.end(parse.gr(event$bp2)))
  else
    evbp2 = gr.end(parse.gr(event$bp2))

  seq1 = ffTrack::get_seq('~/DB/GATK/human_g1k_v37.fasta.2bit', evbp1+pad)
  seq2 = ffTrack::get_seq('~/DB/GATK/human_g1k_v37.fasta.2bit', evbp2+pad)
  res = mclapply(1:length(seq1), function(i, stats)
  {
    message(i)
    if (!anchor)
      win = c(evbp1[i], evbp2[i])+pad
    win$seq = c(seq1[i], seq2[i])
    hom = homeology(win, stride = stride, pad = pad2)
    if (!stats)
      return(hom)
    ix1 = which(hom$gr %^% win[1])
    ix2 = which(hom$gr %^% win[2])
    return(gmstats(gmfeat(hom$mat[ix1,ix2], thresh = thresh))[, seq := i])
  }, stats = stats, mc.cores = mc.cores)
  if (stats)
    res = cbind(event, rbindlist(res, fill = TRUE))
  return(res)
}


#' @name write_ggjs
#' @title dumps ggjs formatted csv file from GRanges
#' @description 
#' 
#' 
#' @param gr GRanges input
#' @param filename filename to dump csv to
#' @param field field name to dump as (numeric) data ["score"]
#' @export
#' @author Marcin Imielinski
write_ggjs = function(gr, filename, field = 'score')
{
  out = data.table(
    chromosome = gr %>% seqnames %>% as.character,
    start = start(gr),
    end = end(gr),
    y =  values(gr)[[field]])

  fwrite(out, filename, sep = ',')
}




#' rel2abs
#'
#' rescales CN values from relative to "absolute" (i.e. per cancer cell copy) scale given purity and ploidy
#'
#' takes in gr with signal field "field"
#'
#' @param gr GRanges input with meta data field corresponding to mean relative copy "mean" in that interval
#' @param purity purity of sample
#' @param ploidy ploidy of sample
#' @param gamma gamma fit of solution (over-rides purity and ploidy)
#' @param beta beta fit of solution (over-rides purity and ploidy)
#' @param field meta data field in "gr" variable from which to extract signal, default "mean"
#' @param field.ncn meta data field in "gr" variable from which to extract germline integer copy number, default "ncn", if doesn't exist, germline copy number is assumed to be zero
#' @return
#' numeric vector of integer copy numbers
#' @export
rel2abs = function(gr, purity = NA, ploidy = NA, gamma = NA, beta = NA, field = 'ratio', field.ncn = 'ncn')
{
  mu = values(gr)[, field]
  mu[is.infinite(mu)] = NA
  w = as.numeric(width(gr))
  w[is.na(mu)] = NA
  sw = sum(w, na.rm = T)
  mutl = sum(mu * w, na.rm = T)

  ncn = rep(2, length(mu))
  if (!is.null(field.ncn))
    if (field.ncn %in% names(values(gr)))
      ncn = values(gr)[, field.ncn]

  ploidy_normal = sum(w * ncn, na.rm = T) / sw  ## this will be = 2 if ncn is trivially 2

  if (is.na(gamma))
    gamma = 2*(1-purity)/purity

  if (is.na(beta))
    beta = ((1-purity)*ploidy_normal + purity*ploidy) * sw / (purity * mutl)
                                        #      beta = (2*(1-purity)*sw + purity*ploidy*sw) / (purity * mutl)


                                        # return(beta * mu - gamma)
  return(beta * mu - ncn * gamma / 2)
}


#' @name abs2rel
#' @description abs2rel
#'
#' rescales CN values from relative to "absolute" (i.e. per cancer cell copy) scale given purity and ploidy
#' By default, output is normalized to 1 (i.e. assumes that the total relative copy number signal mass over the genome is 1)
#'
#' takes in gr with signal field "field"
#' @param gr GRanges input with meta data field corresponding to mean relative copy "mean" in that interval
#' @param purity purity of sample
#' @param ploidy ploidy of sample
#' @param gamma gamma fit of solution (over-rides purity and ploidy)
#' @param beta beta fit of solution (over-rides purity and ploidy)
#' @param field meta character specifying meta data field in "gr" variable from which to extract signal, default "mean"
#' @param field.ncn character specifying meta data field in "gr" variable from which to extract germline integer copy number, default "ncn", if doesn't exist, germline copy number is assumed to be zero
#' @return
#' numeric vector of integer copy numbers
#'  @export
abs2rel = function(gr, purity = NA, ploidy = NA, gamma = NA, beta = NA, field = 'cn', field.ncn = 'ncn', total = 1)
{
  abs = values(gr)[, field]
  w = width(gr)
  sw = sum(as.numeric(w))
  
  ncn = rep(2, length(gr))
  if (!is.null(field.ncn))
    if (field.ncn %in% names(values(gr)))
      ncn = values(gr)[, field.ncn]

  if (is.na(gamma))
    gamma = 2*(1-purity)/purity

  ploidy_normal = sum(w * ncn, na.rm = T) / sw  ## this will be = 2 if ncn is trivially 2

  if (is.na(beta))
    beta = ((1-purity)*ploidy_normal + purity*ploidy) * sw / (purity * total)

                                        #    return((abs + gamma) / beta)
  return((abs + ncn*gamma/2) / beta)
}

#' @name fqfl
#' @description fqfl
#'
#' takes character vector of paths with suffix [(R1)|(R2)].fastq or ...fastq.gz and creates
#' a fastq file list data.table with first two columns R1 and R2 and last column read group id.
#' 
#' @param paths character vector of fastq or fastq.gz paths
#' @return
#' fastq file list data.table with first two columns R1 and R2 and last column read group id.
#'  @export
fqfl = function(paths)
{
  tmp = data.table(path = paths)
  tmp[, read := ifelse(grepl('R1\\.fastq(\\.gz)?$', basename(path)), 'R1',
                ifelse(grepl('R2\\.fastq(\\.gz)?$', basename(path)), 'R2', NA))]
  tmp[, group := gsub('\\.((R1)|(R2))\\.fastq(\\.gz)?$', '', basename(path))]
  dcast.data.table(tmp, group ~ read, value.var = 'path')[, .(R1, R2, group)]
}


#' @name llplot
#' @description Lolliplot wrapper 
#'
#' Adaptation of code from https://www.bioconductor.org/packages/release/bioc/vignettes/trackViewer/inst/doc/trackViewer.html#embl-ebi_proteins_api
#' to plot lolliplot of mutations from SNPEff or VEP annotated VCF.
#'
#' @param vars GRanges of mutations ingested (e.g. via grok_vcf(long = TRUE)) from SnpEff style annotated vcf / bcf file, by default needs to have columns $REF, $ALT, $gene, $protein_pos, $feature_id, $variant.p
#' @param gene gene / protein to plot (default is inferred from vars$gene)
#' @param legend named vector mapping unique values of $type column to colors
#' @param domain_types vector of domain types to include in plot (if NULL will include all)
#' @param type.field name of vars field that corresponds to variant type, will determine the color mapping below (default $type), if blank will infer automatically from annotation.field
#' @param protein_pos.field name of vars field that specifies to variant protein position (default $protein_pos)
#' @param annotation.field name of vars field that specifies variant annotation (default $annotation)
#' @param label.field name of label field (default $variant.p)
#' @param gene.field name of vars field that specifies to variant gene (default $gene)
#' @param ... other parameters to pass to lolliplot e.g. type, yaxis, 
#' @return plots lolliplot with trackViewer::lolliplot
#' @export
llplot = function(variants,
                  gene = variants$gene[1] %>% as.character,
                  domain_types = c("DNA_BIND", "MOTIF", "DOMAIN", "REGION", "BINDING", "CHAIN", "TOPOLOGY"), 
                  legend = NULL,
                  type.field = 'type',
                  gene.field = 'gene',
                  annotation.field = 'annotation', 
                  protein_pos.field = 'protein_pos',
  #                feature_id.field = 'feature_id',
                  label.field = 'variant.p',
                  wes = FALSE,
                  verbose = TRUE,
                  na.col = 'gray90',
                  yaxis = TRUE, 
                  ..., 
                  apiurl = "www.ebi.ac.uk/proteins/api/", taxid = '9606', orgdb = "org.Hs.eg.db" , maxaccession = 20)
{
  if (verbose)
    message('Processing variant data')

  if ( protein_pos.field %in% names(values(variants)) )
  {
    values(variants)$protein_pos =  values(variants)[[protein_pos.field]]
    if (!is.integer(variants$protein_pos))
    {
      variants$protein_pos = sapply(strsplit(as.character(variants$protein_pos), '\\/'), '[', 1) %>% as.integer
    }
  }
  else
    stop(sprintf('%s field does not exist in provided vars object', protein_pos.field))
  
  if ( gene.field %in% names(values(variants)) )
  {
    values(variants)$gene = values(variants)[[gene.field]]
  }
  else
    stop(sprintf('%s field does not exist in provided vars object', gene.field))

  variants = variants[variants$gene %in% gene]

  ## if (length(variants)==0)
  ##   {
  ##     message('empty variant set provided, no plot produced')
  ##     return()
  ##   }

  if (annotation.field %in% names(values(variants)))
  {
    values(variants)$annotation =  values(variants)[[annotation.field]]
  }
  
  if (label.field %in% names(values(variants)))
  {
    values(variants)$label = values(variants)[[label.field]]
  }
  else
    values(variants)$label = ''
  
  if (type.field %in% names(values(variants)))
  {
    if (verbose)
    {
      message('found type annotation in ', type.field, ' field of vars')
    }
    variants$type = values(variants)[[type.field]]
  }
  else ## infer type from annotation field
  {
    trunc = c('disruptive', 'splice_region', 'frameshift', 'stop')
    .parenify = function(x) paste0('(', paste(x, collapse = ')|('), ')')
    variants$truncating = grepl(.parenify(trunc), variants$annotation)
    variants$indel = nchar(variants$REF)!=nchar(variants$ALT)
    variants$type = paste(ifelse(variants$truncating, 'truncating', 'missense'), ifelse(variants$indel, 'indel', 'SNV'))    
  }

  if (is.null(legend))
  {
    legend = unique(variants$type)
    legend = structure(brewer.master(length(legend), 'Spectral', wes = FALSE), names = legend)
    legend = c(legend, c(other = na.col))
    legend = legend[names(legend) %in% unique(names(legend))]
  }
  
  if (any(is.na(variants$type)))
    variants$type[is.na(variants$type)] = 'other'

  variants$tag = paste(seqnames(variants), variants$protein_pos, variants$label)

  pvariants = as.data.table(variants)[, .(score = .N),  by = .(seqnames, tag, start = protein_pos, end = protein_pos, label, color = legend[type])][!is.na(start) & !is.na(end), ] %>% dt2gr
  names(pvariants) = pvariants$label


  library(httr) # load library to get data from REST API
  ## org database to get the uniprot accession id
  eval(parse(text = paste('library(', orgdb,')')))
  eid = BiocGenerics::mget(gene, get(sub(".db", "SYMBOL2EG", orgdb)))[[1]]
  if (verbose)
    message('Requesting UniProt protein domain annotation for gene id ', eid)
  chr = BiocGenerics::mget(eid, get(sub(".db", "CHR", orgdb)))[[1]]
  accession = unlist(lapply(eid, function(.ele){
    BiocGenerics::mget(.ele, get(sub(".db", "UNIPROT", orgdb)))
  }))

  dtypes = ''
  if (!is.null(domain_types))
    dtypes = paste0("&types=", paste(domain_types, collapse = '%2C'))
  featureURL = paste0(apiurl,
                       "features?offset=0&size=-1&reviewed=true",
                       dtypes, 
                       "&taxid=", taxid,
                       "&accession=", paste(accession, collapse = ",")
                       )
  response = GET(featureURL)
  stop_for_status(response)
  content = content(response)

  if (verbose)
    message('Processing protein domain data ', eid)
  content = content[[1]]
  acc = content$accession
  sequence = content$sequence
  domains = rbindlist(lapply(content$features, '[', c('type', 'category', 'description' ,'begin', 'end')))[, seqnames := chr][, start := begin %>% as.integer][, end := end  %>% as.integer] %>% dt2gr
  domains$fill = 1+seq_along(domains)
  names(domains) = domains$description
  domains$height = 0.04

  trackViewer::lolliplot(pvariants, domains, ranges = GRanges(chr, IRanges(1, nchar(sequence))), ylab = gene, legend = legend, yaxis = yaxis, main = gene, ...)
}


#' @name cc
#' @title cc
#' @description select columns of data.frame or data.frame-like object 
#'
#' @param x column expression or variables 
#' @param ... additional variables
#' @export
cc = function(x, y = c(), ...)
{
    if (is.data.table(x))
      eval(parse(text = paste0('x[, ', deparse(substitute(y)), ', ...]')))
    else
      x[, y]
}


#' @name rr
#' @title rr
#' @description select row of data.frame or data.frame- like object data.table
#'
#' @param x expression or variables
#' @param ... additional variables
#' @export
rr = function(x, y = c())
{

  if (is.data.table(x))
    y = tryCatch(eval(eval(parse(text = substitute(deparse(substitute(y)))), parent.frame()),x, parent.frame(2)), error = function(e) NULL)

  if (!is.null(dim(x)))
    x[y, ]
  else
    x[y]
}

#' @name dd
#' @title dd
#' @description "dollar sign" usage of expression
#'
#' @param x variable
#' @export
dd = function(x, y = c())
{
  eval(parse(text = paste0('x$', y)))
}


#' @name variants
#' @title variants
#' @description Call substitutions and indels from contigs
#'
#' Calls substitutions from contigs by comparing DNAStringSet ref
#' against reference sequence ref via RSeqLib::BWA
#'
#' The caller is IUPAC "ambiguity code aware" meaning that for every instance of ambiguity in the query and reference
#' it will output the cartesian product of all mismatching variants (iupac = TRUE) otherwise
#' it will treat those bases literally. 
#'
#' Note: not recommended to use on reads, only for contigs, i.e. will not scale to millions of reads 
#'
#' @param query DNAStringSet of query
#' @param ref DNASTringSet of ref
#' @param expand.iupac logical flag (TRUE) specifying whether to expand iupac for computing SNV  in query and reference
#' @return GRanges in ref coordinates of SNV and indels
#' @author Marcin Imielinski
#' @export
variants = function(query, ref, expand.iupac = TRUE, verbose = FALSE)
{
  nmq = names(query)
  nmr = names(ref)

  if (verbose)
    message(nmq)

  if (is.null(nmq))
    nmq = 1:length(nmq)

  if (is.null(nmr))
    nmr = 1:length(nmr)

  if (!is(query, 'DNAStringSet'))
    query = DNAStringSet(as.character(query))

  if (!is(ref, 'DNAStringSet'))
    ref = DNAStringSet(as.character(ref))

  if (is.null(names(query)))
    names(query) = nmq
  
  if (is.null(names(ref)))
    names(ref) = nmr


  ## replace any gap characters in query
  query = suppressWarnings(replaceAt(query, vmatchPattern('-', query), ''))
 
  bw = RSeqLib::BWA(seq = ref)
  aln = bw[query]

  ## build chain connecting alignment to reference
  cg = cgChain(aln)

  ## create every single base on query
  bases = gr.tile(seqinfo(links(cg)$x), 1)

  ## look up those bases in query
  bases$ALT = query[bases]

  ## lift those bases to reference and compare to sequence to find SNV
  basesl = suppressWarnings(lift(cg, bases) %&% si2gr(ref))
  basesl$REF = ref[basesl]
  basesl$qname = names(query)
  basedt = basesl[, c("REF", "ALT", "qname")] %>% gr2dt

  if (expand.iupac)
  {
    iupac = data.table(
      IUPAC = c('A', 'T', 'G', 'C', 'R', 'R', 'Y', 'Y', 'S', 'S', 'W', 'W', 'K', 'K', 'M', 'M', 'B', 'B', 'B', 'D', 'D', 'D', 'H', 'H', 'H', 'V', 'V', 'V', 'N', 'N', 'N', 'N'),
      base = c('A', 'T', 'G', 'C', 'A', 'G', 'C', 'T', 'G', 'C', 'A', 'T', 'G', 'T', 'A', 'C', 'C', 'G', 'T', 'A', 'G', 'T', 'A', 'C', 'T', 'A', 'C', 'G', 'A', 'C', 'G', 'T')
    )
      
    expand.ref = copy(iupac)
    setnames(expand.ref, 'base', 'ref')
    expand.alt = copy(iupac)
    setnames(expand.alt, 'base', 'alt')

    basedt = basedt %>% merge(expand.alt, by.x = 'ALT', by.y = 'IUPAC', all.x = TRUE, allow.cartesian = TRUE) %>%
      merge(expand.ref, by.x = 'REF', by.y = 'IUPAC', all.x = TRUE, allow.cartesian = TRUE)
    basedt$REF = basedt$ref
  }
    
  snv = basedt[REF != alt, ]

  if (length(snv))
    snv$type = 'SNV'

  ## find unaligned / unmapped chunks  of query to define indels
  insertionsl = basesl[c()]
  insertions = reduce(bases[-basesl$query.id])

  if (length(insertions))
  {
    ## pad either to the left or to the right so that there is something to lift
    start(insertions) = ifelse(start(insertions)>1, start(insertions)-1, start(insertions))
    end(insertions) = ifelse(start(insertions)== 1 & end(insertions) < seqlengths(insertions)[as.character(seqnames(insertions))],
                             end(insertions)+1, end(insertions))
        
    insertions$type = 'INS'
    insertions$ALT = query[insertions]
    insertions$alt = NA
    insertions$qname = as.character(seqnames(insertions))
    insertionsl = lift(cg, insertions)
    if (length(insertionsl))
      insertionsl$REF = ref[insertionsl]
  }

  deletions = basesl[c()]
  vb = suppressWarnings(grl.unlist(varbase(aln)))
  if (length(vb))
    {
      deletions = vb[, c('qname', 'type')] %Q% (type == 'D')
    }

  if (length(deletions))
  {
    deletions$alt = NA
    deletions$ALT = ''
    deletions$REF = ref[deletions]
    deletions$type = 'DEL'
  }
  
  variants = rbind(snv, gr2dt(insertionsl), gr2dt(deletions), fill = TRUE)

  if (nrow(variants))
  {
    variants = variants %>% dt2gr
    strand(variants) = '+'
    variants = variants[, intersect(c('qname', 'type','REF', 'ALT', 'alt'), names(values(variants)))]
  } else
    variants = bases[c(), c()]

  return(sort(variants))
}


#' @name contig.support
#' @title contig.support
#' @description
#'
#' Takes as input a GRanges of bam alignments (e.g. outputted from bamUtils::read.bam) and a GRanges of rearranged
#' reference aligned contigs (e.g. output of RSeqLib::BWA).
#'
#' It identifies the subset of reads that support each of the contigs and "lifts" those reads
#' through the read --> contig and contig --> reference alignments, returning supporting reads in reference coordinates.
#'
#' The criteria for support include min.bases aligning to at least two chunks of the rearranged contig, and
#' requirement that min.aligned.frac fraction of bases in every supporting read is aligned to that contig.
#'
#' Additional requirements for support include not allowing split alignment of individual reads to the contigs
#' (note: this does not mean we don't detect split reads that support the structural variant, this is captured
#' by the contig -> reference alignment, we are just requiring the reads align (near) perfectly to the contig).
#' and requiring alla alignments from a read pair (oriented to R1 frame of the fragment) to align to the same
#' strand of the contig.
#'
#' Finally, reads are not included in support if they align better to the reference than their native alignment,
#' which is determined by comparing the $AS of their contig alignment with their original alignment score, stored
#' in the provided metadata $AS field.  If reference AS is not provided as metadata, it will is assumed to be zero. 
#'
#' $AS can be optionally recomputed against a DNAStringSet "ref" that represent the reference
#' sequence.  (Note that this "ref" does not have to be the full genome reference, it is just used to compute
#' the alignment scores, and in fact for this to work  efficiently, it's recommended that the provided
#' reference sequence is local to the regions of interest, e.g. a few kb flanking each SV breakpoint,
#' rather than the whole genome.)
#'
#' The outputted reads include additional metadata including number of bases aligning to each chunk of the aligned contig.
#' 
#' 
#' @param reads GRanges in SAM / BAM format e.g. output of read.bam or BWA, with fields $qname, $cigar, $flag $seq all populated in standard fashion, and optionally $AS
#' @param contig GRanges in SAM / BAM format wth fields $qname, $cigar and $seq all [populated
#' @param ref optional DNAStringSet representing a reference sequence to compute alignments against
#' @param chimeric logical flag whether to require reads to support junctions in chimericcontigs (ie discontiguous chunks on the reference), chimeric = FALSE
#' @param strict strict requires that the alignment score of the read to contig alignment needs to be better for at least one read (and also not worse for any of the reads) 
#' @param
#' @return reads re-aligned to the reference through the contigs with additional metadata describing features of the alignment
#' @export
#' @author Marcin Imielinski
contig.support = function(reads, contig, ref = NULL, chimeric = TRUE, strict = TRUE, cg.contig = gChain::cgChain(contig), isize.diff = 1e3, min.bases = 20, min.aligned.frac = 0.95, new = TRUE, 
                          verbose = TRUE)
{
  if (length(reads)==0)
    stop('reads must be non empty GRanges with $qname, $cigar, $seq, and $flag fields')

  if (length(contig)==0)
    stop('contigs must be non empty GRanges with $qname, $cigar and $seq fields')

  if (verbose)
    message('Prepping reads for contig alignment')
  seq = unique(gr2dt(contig), by = c('qname'))[, structure(as.character(seq), names = as.character(qname))]
  bwa.contig = RSeqLib::BWA(seq = seq)
  chunks = gChain::links(cg.contig)$x
  strand(chunks) = '+'
  chunks = disjoin(chunks)
  if (is.null(reads$R1))
    reads$R1 = bamflag(reads$flag)[,'isFirstMateRead']>0
  reads$read.id = 1:length(reads)
  if (is.null(reads$AS))
  {
    warning('AS not provided in reads, may want to consider using tag = "AS" argument to read.bam or provide a ref sequence to provide additional specificity to the contig support')
    reads$AS = 0
  }
  nix = as.logical(strand(reads) == '-' )
  reads$seq[nix] = reverseComplement(DNAStringSet(reads$seq[nix])) ## flip read sequences to original strand
  reads[!reads$R1] = gr.flipstrand(reads[!reads$R1]) ## flip R2 read orientation to R1 strand
  reads$seq[!reads$R1] = reverseComplement(DNAStringSet(reads$seq[!reads$R1])) ## flip R2 read sequences to R1 strand
  reads = reads %Q% (!duplicated(paste(qname, R1)))

  if (!is.null(ref)) ## realign reads against reference DNAStringSet if provided to get alignment scores
  {
    if (verbose)
      message('Realigning reads against reference DNAStringSet')
    bwa.ref = RSeqLib::BWA(seq = ref)
    tmp = bwa.ref[reads$seq] %>% gr2dt
    tmp$ix = as.numeric(as.character(tmp$qname))
    tmp$R1 = reads$R1[tmp$ix]
    tmp$qname = reads$qname[tmp$ix]
    tmp = unique(tmp, by = c('qname', 'R1'))
    setkeyv(tmp, c('qname', 'R1'))
    if (nrow(tmp))
    {
      tmp[, isize := ifelse(any(seqnames != seqnames[1] | any(strand != strand[1])), NA_integer_, diff(range(start, end))), by = qname]
      reads$isize = pmin(tmp[.(reads$qname, reads$R1), isize], reads$isize, Inf, na.rm = TRUE)
      reads$AS = tmp[.(reads$qname, reads$R1), AS]
    }
  }

  if (verbose)
    message('Aligning reads against derivative contigs')
  

  ## aligning reads to contig
  rdt = as.data.table(reads)
  rdt[, ref.aligned := countCigar(cigar)[, 'M']]
  rdt[, ref.aligned.frac := ref.aligned/qwidth[1], by = .(qname, R1)]

  reads$ref.aligned.frac = rdt$ref.aligned.frac

  readsc = bwa.contig[reads$seq] %>% gr2dt
  readsc$cigar = as.character(readsc$cigar)
  readsc$ix = as.integer(as.character(readsc$qname))
  readsc$R1 = reads$R1[readsc$ix]
  readsc$read.id = reads$read.id[readsc$ix]



  ## these are splits on the contig, not reference --> shouldn't be any for good alignment
  readsc[, nsplit := .N, by = .(qname, R1)] 
  readsc[, aligned := countCigar(cigar)[, 'M']]

  ## these are splits on the contig, not reference --> shouldn't be any for good alignment
  readsc[, aligned.frac := aligned/qwidth[1], by = .(qname, R1)]
  readsc$AS.og = reads$AS[readsc$ix]
  readsc$isize = abs(reads$isize[readsc$ix])


  readsc$seqnames.og = seqnames(reads)[readsc$ix] %>% as.character
  readsc$strand.og = strand(reads)[readsc$ix] %>% as.character
  readsc$start.og = start(reads)[readsc$ix]
  readsc$end.og = end(reads)[readsc$ix]
  readsc$ref.isize = gr2dt(readsc)[, ref.isize := ifelse(
                                       all(seqnames.og == seqnames.og[1]) & all(strand.og == strand.og[1]),
                                       as.numeric(diff(range(c(start.og, end.og)))),                                   
                                       Inf), by = qname]$ref.isize %>% abs

  readsc$ref.aligned.frac = reads$ref.aligned.frac[readsc$ix]
  readsc$AS.og[is.na(readsc$AS.og)] = 0
  readsc$qname = reads$qname[readsc$ix]


  ## new scoring method based on cgChain of reads to contigs
  if (new)
    {
      ## cgChain representing read to contig alignments
      readsc$al.id = 1:nrow(readsc)

      if (verbose)
        message('Generating read to contig cgChain')
      alcg = gChain::cgChain(readsc)
      alchunks = cbind(as.data.table(values(alcg)), as.data.table(gChain::links(alcg)$x), as.data.table(gChain::links(alcg)$y)[, .(contig = seqnames, contig.start = start, contig.end = end, contig.strand = strand)])

      ## strands should be aligned to read / fragment + strand, but if not let's flip
      alchunks[strand == '-', ":="(strand = '+', contig.strand = c('+' = '-', '-' = '+')[contig.strand])]

      ## now for each al.id (ie bam record) let's pick the left most gChain / links record on the read / fragment 
      ## ie this is the lowest coordinate on the query
      ## (note that cgChain will split indels into separate ranges hence giving one to many mapping of al.id
      ## to records in links)
      setkeyv(alchunks, c('qname', 'start', 'end'))
      ## alchunks[, is.min := start == min(start), by = al.id]
      ## alchunks = alchunks[is.min == TRUE, ]

      ## so now we want to find alignments that are
      ## (1) concordant with respect to the contig
      ##  i.e. there is a monotonic increase (decrease) of contig.start if the contig.strand is + (-)
      ## (2) most of the read (aligned.frac) is represented
      ## (3) AS scores are better than original
      ## (4) isize better than original (where isize is the contig. span between the first and last alignment) .. related to (1)


      if (verbose)
        message('Scoring read to contig to alignments')
      alchunks[, contig.sign := ifelse(contig.strand == '+', 1, -1)]
      alchunks[, concordant.sign := all(contig.sign == contig.sign[1]), by = qname]

      ## check if we never go from R1 == FALSE to R1 == TRUE
      alchunks[, concordant.R1R2 := all(diff(!R1)>=0), by = qname]

      ## check to see that our contig.start always increasing or decreasing
      alchunks[, concordant.start := all(diff(contig.sign[1]*contig.start)>0), by = qname]

      alchunks[, contig.isize := diff(range(contig.start, contig.end)), by = qname]
      alchunks[, bases := sum(width), by = qname]

      alchunks[, AS.better := sum(width[AS>AS.og]), by = qname]
      alchunks[, AS.worse := sum(width[AS<AS.og]), by = qname]
      alchunks[, AS.equal := sum(width[AS==AS.og]), by = qname]

      keepq = alchunks[concordant.sign & concordant.R1R2 & concordant.start &
                       bases > min.bases & aligned.frac > min.aligned.frac & aligned.frac >= ref.aligned.frac & 
                      (AS.better>0 | contig.isize<ref.isize) & AS.worse == 0, ]

      
      keepq = keepq[, .(qname, contig, contig.isize, contig.strand, bases, contig.sign, AS.better, AS.worse, AS.equal)] %>% unique(by = 'qname')
    }
  else ## old scoring method
  {
    ## if strict (default) remove any alignments that overlap others in the same qname
    if (strict)
    {
      readsc = dt2gr(readsc)
      readsc = readsc %Q% (rev(order(AS)))
      readsc = readsc[!duplicated(gr.match(readsc, readsc, by = 'qname')), ] %>% gr2dt
    }


    if (verbose)
      message('Computing overlap stats')

    ov = dt2gr(readsc) %*% chunks
    strand(ov) = readsc$strand[ov$query.id]
    ov$subject.id = paste0('chunk', ov$subject.id)
    ovagg = dcast.data.table(ov %>% gr2dt, qname ~ subject.id, value.var = 'width', fun.aggregate = sum)
    ovagg$nchunks = rowSums(ovagg[, -1]>min.bases)  ## good means we hit multiple chunks with sufficient bases
    rstats = gr2dt(ov)[, .(
                    contig.id = unique(seqnames)[1],
                    pos = sum(width[strand == '+']),
                    neg = sum(width[strand == '-']),
                    aligned.frac = min(aligned.frac),
                    num.contigs = length(unique(seqnames)), ### fixing later ... multiple contigs as input could distort results
                    paired = any(R1) & any(!R1), 
                    isize.contig = diff(range(c(start, end))),
                    isize.og = isize[1],
                    qsplit = any(nsplit>1), ## any sequences in this qname split on the contig ie a bad alignment on the contig
                    worse = any(AS.og>AS), ## any alignment in this qname worse than vs reference?
                    better = any(AS>AS.og) ## any alignment in this qname better than vs reference?
                  ), by = qname] %>% merge(ovagg, by = 'qname')

    ## apply filters ie nchunks>1 if chimeric, all alignments have to be of one sign
    ## if not paired then AS < AS.og else isize<isize.og
    keepq = rstats[nchunks>chimeric & (pos == 0 | neg  == 0) & aligned.frac > min.aligned.frac & !worse & (better | !strict | (paired & isize.contig < isize.og - isize.diff)) & !qsplit & num.contigs == 1, ]
    if (nrow(keepq)==0)
      return(reads[c()])

    keepq$aligned.frac = NULL
  }

  readsc = merge(readsc, keepq, by = 'qname') %>% dt2gr
  
  if (verbose)
    message('Lifting reads through contig back to reference')

  out = gChain::lift(cg.contig, readsc)

  if (length(out)) ## add reads metadata back to out
  {
    out[!out$R1] = gr.flipstrand(out[!out$R1])
    out$col = ifelse(out$R1, 'blue', 'gray')

    if (verbose)
      message('Adding metadata to reads')
    metacols = setdiff(names(values(reads)), names(values(out)))
    values(out) = cbind(values(out), values(reads)[match(out$read.id, reads$read.id), metacols])
  }

  if (verbose)
    message('Done')
  out
}

#' @name junction.support
#' @title junction.support
#' @description
#'
#' Takes as input a GRanges of bam alignments (e.g. outputted from bamUtils::read.bam) and a GRanges of rearranged
#' reference aligned contigs (e.g. output of RSeqLib::BWA) and a set of Junction objects, and outputs reads supporting
#' these junctions by building a contig around each junction (from the reference) and then running contig.support (see
#' that functions docuemntation for criteria)
#'
#' @param reads GRanges in SAM / BAM format e.g. output of read.bam or BWA, with fields $qname, $cigar, $flag $seq all populated in standard fashion, and optionally $AS
#' @param junctions Junction object
#' @param bwa RSeqLib BWA object and path to fasta file corresponding to the reference
#' @param ref optional DNAStringSet corresponding to reference genome sequence
#' @param pad padding around the junction breakpoint around  which to analyze contig and reference sequences, this should be several standard deviations above the average insert size (2000)
#' @param realign flag whether to realign or just use existing alignments
#' @param bx logical flag whether data is linked reads, must then have BX flag, and the pad will be set to minimum 1e5
#' @param verbose logical flag (TRUE)
#' @param ... additional parameters to contig support
#' @return reads re-aligned to the reference through the contigs with additional metadata describing features of the alignment
#' @export
#' @author Marcin Imielinski
junction.support = function(reads, junctions = NULL, bwa = NULL, ref = NULL, pad = 500, bx = FALSE, pad.ref = pad*20, realign = TRUE, walks = NULL, verbose = TRUE, ...)
{

  if (!inherits(reads, 'GRanges') || is.null(reads$qname) || is.null(reads$cigar) || is.null(reads$seq) || is.null(reads$flag))
    stop('read input must be GRanges with fields $qname, $cigar, $seq, $flag and optionally $AS')

  if (bx)
    pad = max(pad, 1e5)

  if (!is.null(junctions))
    walks = jJ(junctions$grl)$gw(pad = pad)

  if (is.null(walks))
    stop('Either walks or junctions must be provided')

  if (bx)
  {
    if (is.null(reads$BX))
      stop('reads must have BX tag, may need to read.bam with tag option to extract it')

    if (!length(reads))
      return(reads)

    sc = score.walks(walks$grl, reads = reads, verbose = FALSE, raw = TRUE)$sc
    res = as.data.table(melt(as.matrix(sc)))[value>0, .(BX = Var1, walk = Var2)]
    reads = gr2dt(reads) %>% merge(res, by = 'BX') %>% dt2gr
    return(reads)
  }

  if (!realign)
  {
    if (is.null(junctions))
      junctions = walks$edges$junctions

    ## strand flip since 
    ## read orientation convention
    ## is opposite to junction convention
    reads = gr.flipstrand(reads) 
    reads$R1 = bamUtils::bamflag(reads$flag)[,'isFirstMateRead']>0
    r1 = reads %Q% (R1 == TRUE) %>% as.data.table
    r2 = reads %Q% (R1 == FALSE) %>% as.data.table
    ov = merge(r1, r2, by = 'qname')
    sl = seqlengths(reads)
    grl = grl.pivot(
      GRangesList(dt2gr(ov[, .(seqnames = seqnames.x, start = start.x, end =end.x, strand = strand.x)],
                        seqlengths = sl),
                  dt2gr(ov[, .(seqnames = seqnames.y, start = start.y, end = end.y, strand = strand.y)],
                        seqlengths = sl)))
    values(grl)$qname = ov$qname
    ## make junctions out of reads and cross with "real" junctions
    jn = merge(jJ(grl), junctions, cartesian = TRUE, pad = pad)
    if (!length(jn))
      return(reads[c()])
    out = merge(as.data.table(gr.flipstrand(reads)), unique(jn$dt[, .(qname, junction.id = subject.id)]), by = 'qname') %>% dt2gr(seqlengths = sl)
    return(out)
  }
  
  if (inherits(bwa, 'character') && file.exists(bwa))
  {
    if (verbose)
      message('Loading BWA index')
    bwa = BWA(bwa)
  }

  if (!inherits(ref, 'DNAStringSet'))
  {
    if (verbose)
      message('Loading genome reference as DNAStringSet')

    ref = rtracklayer::import(bwa@reference)
  }

  ## only use the fasta header before the first space as the seqnames of ref 
  names(ref) = strsplit(names(ref), '\\s+') %>% sapply('[', 1)

  if (length(setdiff(seqnames(walks$nodes$gr), seqlevels(ref))))
    stop('seqlevels mismatch between junctions / walks and reference, please check reference (e.g. chr issues)')

  if (length(setdiff(seqnames(walks$nodes$gr), seqlevels(bwa))))
    stop('seqlevels mismatch between junctions / walks and BWA reference, please check reference (e.g. chr issues)')

  if (verbose)
    message('Building and mapping derivative contigs')

  contig = bwa[ref[gr.fix(walks$grl, ref, drop = TRUE)]]

  if (verbose)
    message('Building reference contigs flanking junctions')
  contigref = ref[gr.fix(walks$edges$junctions$footprint + pad.ref, ref, drop = TRUE)]


  if (verbose)
    message('Making gChain mapping contigs to reference')
  cg.contig = gChain::cgChain(contig)

  if (verbose)
    message('Running contig support')

  reads = contig.support(reads, contig, ref = contigref, cg.contig = cg.contig, ...)
  reads$junction.id = reads$contig.id
  return(reads)  
}

#' @name memu
#' @title memu
#' @description
#' check memory usage by user on current server
#' @export
memu = function()
{
  res = pipe('~/scripts/mem')  %>% readLines %>% paste(collapse='\n') %>% fread
  setnames(res, c('user', 'GB'))
  return(res[rev(order(GB)), ])
}


#' @name oncotable
#' @title oncotable
#' @description
#'
#' Takes as input (keyed) "tumors" (aka pairs) table which a metadata table with specific
#' columns pointing to paths corresponding to one or more of the following pipeline outputs:
#'
#' $annotated_bcf  Path to annotated.bcf file that is the primary output of SnpEff module from which TMB and basic mutation
#' descriptions are extracted along with their basic categories (these will comprising the core of the oncoplot are computed)
#' 
#' $fusions  Path to fusion.rds file that is the primary output of the Fusions modjle, from which protein coding fusions will
#' be computed for
#' 
#' $jabba_rds  Path to jabba.simple.rds output representing primary output of JaBbA module from which SCNA and overall
#' junction burden are computed
#' 
#' $complex    Path to complex.rds gGnome cached object that is the primary output of Events module, from which simple
#' and complex event burdens are computed
#' 
#' $signature_counts Path to signature_counts.txt that is the primary output of Signatures module from which SNV signature
#' counts are computed
#' 
#' The function then outputs a melted data.table of "interesting" features that can be saved and/or immediately output
#' into oncoprint.  This data.table will at the very least have fields $id $type (event type), $track, and  $source
#' populated in addition to a few other data type specific columns.
#'
#' The $source column is the name of the column of tumors from which that data was extracted, and track is a grouping
#' variable that allows separation of the various data types. 
#'
#' All the paths above can be NA or non existent, in which case a dummy row is inserted into the table so that downstream
#' applications know that data is missing for that sample. 
#'
#' @param tumors keyed data.table i.e. keyed by unique tumor id with specific columns corresponding to  paths to pipeline outputs(see description)
#' @param gencode path to gencode .gtf or .rds with GRanges object, or a GRanges object i.e. resulting from importing the (appropriate) GENCODE .gtf via rtracklayer, note: this input is only used in CNA to gene mapping
#' @param amp.thresh SCNA amplification threshold to call an amp as a function of ploidy (4)
#' @param del.thresh SCNA deletion threshold for (het) del as a function of ploidy (by default cn = 1 will be called del, but this allows additoinal regions in high ploidy tumors to be considered het dels)
#' @param mc.cores number of cores for multithreading
#' @param verbose logical flag 
#' @author Marcin Imielinski
#' @export
oncotable = function(tumors, gencode = NULL, verbose = TRUE, amp.thresh = 4, filter = 'PASS', del.thresh = 0.5, mc.cores = 1)
{
  if (is.null(gencode))
    gencode = skidb::read_gencode()
  else if (is.character(gencode))
  {
    if (grepl('.rds$', gencode))
      gencode = readRDS(gencode)
    else
      gencode = rtracklayer::import(gencode)
  }

  pge = gencode %Q% (type  == 'gene' & gene_type == 'protein_coding')

  .oncotable = function(dat, x = dat[[key(dat)]][1], pge, verbose = TRUE, amp.thresh = 4, del.thresh = 0.5, filter = 'PASS')
  {
    out = data.table()

    ## collect gene fusions
    if (!is.null(dat$fusions) && file.exists(dat[x, fusions]))
    {
      if (verbose)
        message('pulling $fusions for ', x)
      fus = readRDS(dat[x, fusions])$meta
      if (nrow(fus))
      {
        fus = fus[silent == FALSE, ][!duplicated(genes), ]
        fus[, vartype := ifelse(in.frame == TRUE, 'fusion', 'outframe_fusion')] # annotate out of frame fusions
        fus = fus[, .(gene = strsplit(genes, ',') %>% unlist, vartype = rep(vartype, sapply(strsplit(genes, ','), length)))][, id := x][, track := 'variants'][, type := vartype][, source := 'fusions']
        out = rbind(out, fus, fill = TRUE, use.names = TRUE)
      }
    } 
    else ## signal missing result
      out = rbind(out, data.table(id = x, type = NA, source = 'fusions'), fill = TRUE, use.names = TRUE)

    ## collect complex events
    if (!is.null(dat$complex) && file.exists(dat[x, complex]))
    {
      if (verbose)
        message('pulling $complex events for ', x)
      sv = readRDS(dat[x, complex])$meta$events
      if (nrow(sv))
      {
        sv = sv[, .(value = .N), by = type][, id := x][, track := ifelse(type %in% c('del', 'dup', 'invdup', 'tra', 'inv'), 'simple sv', 'complex sv')][, source := 'complex']
        out = rbind(out, sv, fill = TRUE, use.names = TRUE)
      }
    }
    else
      out = rbind(out, data.table(id = x, type = NA, source = 'complex'), fill = TRUE, use.names = TRUE)

    ## collect copy number / jabba
    if (!is.null(dat$jabba_rds) && file.exists(dat[x, jabba_rds]))
    {
      if (verbose)
        message('pulling $jabba_rds to get SCNA and purity / ploidy for ', x)
      jab = readRDS(dat[x, jabba_rds])
      out = rbind(out,
                  data.table(id = x, value = c(jab$purity, jab$ploidy), type = c('purity', 'ploidy'), track = 'pp'),
                  fill = TRUE, use.names = TRUE)

      gg = gG(jab = jab)

      # get the ncn data from jabba
      kag = readRDS(dat[x, gsub("jabba.simple.rds", "karyograph.rds", jabba_rds)])
      ngr = gg$nodes$gr
      if ('ncn' %in% names(mcols(kag$segstats))){
          ngr = ngr %$% kag$segstats[, c('ncn')]
      } else {
          # if there is no ncn in jabba then assume ncn = 2
          ngr$ncn = 2
      }
      ndt = gr2dt(ngr)

      # we will use the normal ploidy to determine hetdels 
      # so instead of a cutoff of del.thresh * ploidy, we use:
      # del.thresh * ploidy * ncn / normal_ploidy
      # where ncn is the local normal copy number
      seq_widths = as.numeric(width(ngr))
      # since we are comparing to CN data which is integer then we will also round the normal ploidy to the nearest integer.
      normal_ploidy = round(sum(seq_widths * ngr$ncn, na.rm = T) / sum(seq_widths, na.rm = T))

      scna = rbind(
        ndt[cn>=amp.thresh*jab$ploidy, ][, type := 'amp'],
        ndt[cn < ncn | cn<del.thresh*jab$ploidy*ncn/normal_ploidy, ][, type := 'hetdel'],
        ndt[cn == 0, ][, type := 'homdel']
      )

      if (nrow(scna))
      {
        scna = dt2gr(scna, seqlengths = seqlengths(gg)) %*% pge[, 'gene_name'] %>% gr2dt

        if (nrow(scna))
        {
          scna[, track := 'variants'][, source := 'jabba_rds'][, vartype := 'scna']
          out = rbind(out,
                      scna[, .(id = x, value = cn, type, track, gene = gene_name)],
                      fill = TRUE, use.names = TRUE)
        }
      }
    }
    else
      out = rbind(out, data.table(id = x, type = NA, source = 'jabba_rds'), fill = TRUE, use.names = TRUE)

    ## collect signatures
    if (!is.null(dat$signature_counts) && file.exists(dat[x, signature_counts]))
    {
      if (verbose)
        message('pulling $signature_counts for ', x)
      sig = fread(dat[x, signature_counts])
      sig = sig[, .(id = x, value = num_events, type = Signature, etiology = Etiology, frac = frac.events, track = 'signature', source = 'signature_counts')]
      out = rbind(out, sig, fill = TRUE, use.names = TRUE)
    }
    else
      out = rbind(out, data.table(id = x, type = NA, source = 'signature_counts'), fill = TRUE, use.names = TRUE)

    ## collect gene mutations
    if (!is.null(dat$annotated_bcf) && file.exists(dat[x, annotated_bcf]))
    {
      if (verbose)
        message('pulling $annotated_bcf for ', x, ' using FILTER=', filter)
      bcf = grok_bcf(dat[x, annotated_bcf], label = x, long = TRUE, filter = filter)
      if (verbose)
        message(length(bcf), ' variants pass filter')
      genome.size = sum(seqlengths(bcf))/1e6
      nmut = data.table(as.character(seqnames(bcf)), start(bcf), end(bcf), bcf$REF, bcf$ALT) %>% unique %>% nrow
      mut.density = data.table(id = x, value = c(nmut, nmut/genome.size), type = c('count', 'density'),  track = 'tmb', source = 'annotated_bcf')
      out = rbind(out, mut.density, fill = TRUE, use.names = TRUE)
      keepeff = c('trunc', 'cnadel', 'cnadup', 'complexsv', 'splice', 'inframe_indel', 'fusion', 'missense', 'promoter', 'regulatory','mir')
      bcf = bcf[bcf$short %in% keepeff]
      if (verbose)
        message(length(bcf), ' variants pass keepeff')
      vars = NULL
      if (length(bcf))
      {
        bcf$variant.g = paste0(seqnames(bcf), ':', start(bcf), '-', end(bcf), ' ', bcf$ALT, '>', bcf$REF)
        vars = gr2dt(bcf)[, .(id = x, gene, vartype, variant.g, variant.p, distance, annotation, type = short, track = 'variants', source = 'annotated_bcf')] %>% unique
      }
      out = rbind(out, vars, fill = TRUE, use.names = TRUE)
    }
    else
      out = rbind(out, data.table(id = x, type = NA, source = 'annotated_bcf'), fill = TRUE, use.names = TRUE)

    if (verbose)
      message('done ', x)

    return(out)
  }

  if (is.null(key(tumors)))
  {
    if (is.null(tumors$id))
      stop('Input tumors table must be keyed or have column $id')
    else
      setkey(tumors, id)
  }

  out = mclapply(tumors[[key(tumors)]], .oncotable,
                 dat = tumors, pge = pge, amp.thresh = amp.thresh, filter = filter, del.thresh = del.thresh, verbose = verbose, mc.cores = mc.cores)
  out = rbindlist(out, fill = TRUE, use.names = TRUE)

  setnames(out, 'id', key(tumors))
  return(out)
}



#' @name oncoprint
#' @title oncoprint
#' @description
#'
#' Simple wrapper around to oncoPrint from complexHeatmap package to allow quick plotting
#' of patients x genes + metadata.  Uses the data gathered by oncotab to generate quick simple
#' plots that include a core matrix of genes x patients containing data on  SCNA, (complex), fusions,
#' SNV, and indels.  Additional tracks plotting log TMB + 1, log SV burden, complex events,
#' SNV signatures can be provided. 
#' 
#'
#' @param tumors  keyed table of tumors (aka pairs table) with field $oncotable which points to a cached .rds file of an oncotable e.g. produced by oncotable function or Oncotable module / task
#' @param oncotab output from oncotable function with field $id
#' @param genes character vector of genes
#' @param columns additional columns of tumors matrix to plot as horizontal tracks below the main track
#' @param split character of name of column in tumors table to split on (NULL)
#' @param sort  logical flag whether to automatically sort rows i.e. genes and columns i.e. tumors in a "stair step" pattern or default to the provided (TRUE)
#' @param noncoding logical flag whether to show non protein coding mutations
#' @param sort.genes logical flag whether to sort rows i.e. genes with respect to their frequency (TRUE)
#' @param sort.tumors logical flag whether to sort columns i.e. patients in a stairstep pattern with respect to the provided gene order (TRUE)
#' @param sv.stack  logical flag whether to stack bar plot simple and complex SV event counts (FALSE)
#' @param signatures logical flag whether to show signatures (if data is provided / available) (TRUE)
#' @param svevents logical flag whether to show events (if data is provided / available) (TRUE)
#' o=@param tmb logical flag whether to show TMB bar plot (TRUE)
#' @param tmb.log  logical flag whether to log TMB + 1 (TRUE)
#' @param pp logical flag whether to show purity / ploidy (if data is provided / available) (TRUE)
#' @param ppdf whether to print to pdf via ppdf
#' @param track.height height of tracks in 'cm'
#' @param split.gap  gap between splits
#' @param signature.main integer indices of main COSMIC signatures to keep
#' @param signature.thresh lower threshold for non main signature fraction in at least one sample to plot
#' @param outframe.fusions show fusions that are out-of-frame (FALSE)
#' @param cex length 1 or 2 vector canvas expansion factor to apply to the oncoprint itself (relative to 10 x 10 cm) (c(1,3))
#' @param return.mat whether to return.mat
#' @param wes logical flag whether to use wesanderson coolors
#' @param mc.cores multicore threads to use for $oncotable loading from tumors table (not relevant if oncotab provided)
#' @param ... other arguments to ppdf
#' @return ComplexHeatmap object (if ppdf = FALSE), and genotype matrix (if return)
#' @author Marcin Imielinski
#' @export 
oncoprint = function(tumors = NULL,
                     oncotab = NULL,
                     genes = c('KRAS', 'EGFR', 'BRAF', 'TP53', 'TERT', 'CCND1', 'MYC', 'PIK3CA', 'PTEN', 'CDKN2A', 'ARID1A', 'SMARCA4'),
                     split = NULL, 
                     sort = TRUE, sort.genes = sort, sort.tumors = sort,
                     columns = NULL,
                     noncoding = FALSE,
                     cna = TRUE, tmb = TRUE, pp = TRUE, signature = TRUE, svevents = TRUE, basic = FALSE, 
                     ppdf = TRUE,
                     return.oncotab = FALSE,
                     return.mat = FALSE,                     
                     wes = TRUE,
                     drop = TRUE,
                     drop.genes = FALSE, 
                     track.height = 1,
                     signature.thresh = 0.2,
                     signature.main = c(1:5,7,9,13),
                     outframe.fusions = FALSE,
                     track.gap = track.height/2,
                     split.gap = 1,
                     colnames.fontsize = 10,
                     rownames.fontsize = 10,
                     track.fontsize = 10,
                     mc.cores = 1,
                     verbose = FALSE,
                     height = 20,
                     width = 20,
                     ...)
{

  if (basic)
    tmb = svevents = signature = FALSE

  if (!length(genes))
    stop('genes must be provided either as a vector or named list of gene identifiers')

  if (is.list(genes))
    genes = dunlist(genes)[, .(genes = V1, group = listid)]
  else
    genes = data.table(genes = genes, group = NA)

  genes = genes[!duplicated(genes), ]

  if (!is.null(tumors))
  {
    if (!is.null(key(tumors)))
      tumors$id = tumors[[key(tumors)]]

    if (is.null(tumors$id))
      stop('tumors be either keyed or have $id field, if you are resorting e.g. manually sorting your input table the key may get lost so then you should set an $id field explicitly')
    
    if (any(duplicated(tumors$id)))
      stop('check key field in tumors table: duplicated ids present. The key should be unique per row, and matched to the $id field of oncotab')
  }

  missing = c()
  if (is.null(oncotab))
  {
    errmsg = 'Either oncotab or tumors argument must be provided, where tumors is a keyed data table (where each row is a tumor) with column $oncotable of file paths pointing to the cached rds Oncotable results for each tumors'
    if (is.null(tumors) || is.null(tumors$oncotable))
      stop(errmsg)

    fe = file.exists(tumors$oncotable)
    missing = union(missing, tumors$id[!fe])

    if (any(!fe))
      warning(paste(sum(!fe), 'of', length(fe), 'tumors with missing oncotab, will remove if drop = TRUE, otherwise mark'))

    if (!nrow(tumors))
      stop('No tumors with $oncotable field pointing to existing path')

    if (verbose)
      message('Scraping $oncotable paths for oncotable .rds files.  To speed up, consider multi-threading with mc.cores and if you will be creating multiple plots.  Also consider running this with return.oncotab = TRUE and use that for subsequent calls via oncotab = argument.')

    oncotab = mclapply(which(fe), function(x) {y = readRDS(tumors$oncotable[x]); if (nrow(y)) y[, id := tumors$id[x]]; return(y)}, mc.cores = mc.cores) %>% rbindlist(fill = TRUE)
    oncotab$id = factor(oncotab$id, tumors$id)    
  }

  if (!is.null(tumors))
    {
      oncotab$id = factor(oncotab$id, tumors$id)
      missing = union(missing, setdiff(tumors$id, oncotab$id))
    }
  else
    oncotab$id = factor(oncotab$id)
  
  oncotab = oncotab[!is.na(id), ]

  if (!nrow(oncotab))
  {
    if (!is.null(tumors))
      stop('empty oncotable provided, check tumors table, there may be an id mismatch or no non empty files')
    else
      stop('empty oncotable provided, please check inputs')
  }

  ## temp FIX: remove hetdel for genes that have a homdel --> need to fix oncotable itself to deal with this
  if (any(ix <- oncotab$type == 'homdel'))
  {
    oncotab[, rem := FALSE]
    oncotab[type %in% c('amp', 'hetdel', 'homdel'), rem := type == 'hetdel' & any(type == 'homdel'), by = .(gene, id)]
    oncotab = oncotab[rem == FALSE, ]
    oncotab$rem = NULL
  }

  vars = oncotab[track == 'variants', ][gene %in% genes$genes, ][type != 'synonymous', ]

  ## keep track of missing samples ie those that had either SNV, jabba, fusions
  ## will get a gray column in the plot
  missing = union(missing, vars[track == 'variants' & is.na(type), id])

  if (!noncoding)
    vars = vars[!(type %in% c('promoter', 'noncoding', 'regulatory')), ]

  if (!cna)
    vars = vars[!(type %in% c('amp', 'hetdel', 'homdel')), ]  

  vars[, gene := factor(gene, genes$genes)]
  vars = vars[!is.na(gene), ]

  ## convert to matrix format for complex heatmap
  if (nrow(vars))
    {
      varc = dcast.data.table(data = vars, gene ~ id, value.var = "type", fill = 'WT', drop = FALSE, fun.aggregate = function(x) paste(x, collapse = ','))
      varm = as.matrix(varc[, -1])
      rownames(varm) = varc$gene
    }
  else
  {
    varm = matrix('WT', nrow = length(levels(vars$gene)), ncol = length(levels(vars$id)),
           dimnames = list(levels(vars$gene), levels(vars$id)))
  }

  ## prune / label missing genotypes (ie either due to missing or incomplete oncotable entries)
  if (length(missing))
    {
      if (!drop)
        varm[, intersect(colnames(varm), missing)] = 'missing'
      else
        varm = varm[, setdiff(colnames(varm), missing)]
    }

  ## then gene binary order
  if (sort.genes)
    {
      ##ix = skitools::border(varm!='') %>% rev
      ix = rev(order(rowSums(varm!='WT' & varm != 'missing', na.rm = TRUE)))
      varm = varm[ix, , drop = FALSE]
    }
  
  ## then sample binary mutation order
  if (sort.tumors)
    {
      jx = rev(skitools::border(t(varm)!='WT' & t(varm) != 'missing'))
      varm = varm[, jx, drop = FALSE]
    }
    
  ## customize appeagrid appearance with mix of rectangles and circles
  ord = c("amp", "hetdel", "homdel", 'trunc', 'splice', 'inframe_indel', 'fusion', 'missense', 'promoter', 'regulatory')
  if (outframe.fusions == TRUE){
      ord = c("amp", "hetdel", "homdel", 'trunc', 'splice', 'inframe_indel', 'outframe_fusion', 'fusion', 'missense', 'promoter', 'regulatory')
  }
  alter_fun = function(x, y, w, h, v) {
    CSIZE = 0.25
    w = convertWidth(w, "cm")*0.7
    h = convertHeight(h, "cm")*0.7
    l = min(unit.c(w, h))
    grid.rect(x, y, w, h, gp = gpar(fill = alpha("grey90", 0.4), col = NA))
    v = v[ord]
    for (i in which(v)) {
      if (names(v)[i] %in% c('amp', "hetdel", "homdel", 'fusion', 'outframe_fusion'))
        grid.rect(x,y,w,h, gp = gpar(fill = varcol[names(v)[i]], col = NA))
      else if (grepl("missing", names(v)[i]))
        grid.rect(x, y, w, h, gp = gpar(fill = varcol[names(v)[i]], col = NA))
      else if (grepl("trunc", names(v)[i]))
        {
          grid.segments(x - w*0.5, y - h*0.5, x + w*0.5, y + h*0.5,
                        gp = gpar(lwd = 2, col = varcol[names(v)[i]]))
          grid.segments(x - w*0.5, y + h*0.5, x + w*0.5, y - h*0.5,
                        gp = gpar(lwd = 2, col = varcol[names(v)[i]]))
        }
      else if (grepl("(missense)|(promoter)|(regulatory)", names(v)[i]))
      {
        grid.circle(x,y,l*CSIZE, gp = gpar(fill = varcol[names(v)[i]], col = NA))
      }
      else {
        if (grepl("indel", names(v)[i]))
          grid.rect(x,y,w*0.9,h*0.4, gp = gpar(fill = varcol[names(v)[i]], col = NA))
      }
    }
  }

  varcol = c(
    WT = alpha('gray', 0),
    fusion = alpha('green', 0.5),
    outframe_fusion = alpha('greenyellow', 0.5),
    hetdel = 'lightblue',
    missing = 'gray',            
    amp = "red",
    drop = FALSE,
    homdel = "darkblue",
    missense = 'gray40',
    inframe_indel = 'darkgreen',
    promoter  = alpha('red', 0.5),
    regulatory  = alpha('red', 0.2),
    trunc = alpha("blue", 0.8),
    mir = alpha('purple', 0.4),
    splice = "purple"
  )
  
  ids = colnames(varm)
  out.mat = varm ## in case we want to return.mat

  ## generate additional plots if requested / available
  bottom_data = top_data = list()
  if (tmb & any(oncotab$track == 'tmb'))
  {
    tmbd = oncotab[track == 'tmb' & type == 'density', structure(value, names = as.character(id))][ids]
    
    top_data$TMB = tmbd
    out.mat = rbind(TMB = tmbd, out.mat)
  }

  if (pp & any(oncotab$track == 'pp'))
  {
    top_data$Purity = oncotab[track == 'pp' & type == 'purity', structure(value, names = as.character(id))][ids]
    top_data$Ploidy = oncotab[track == 'pp' & type == 'ploidy', structure(value, names = as.character(id))][ids]

    out.mat = rbind(Purity = top_data$Purity, Ploidy = top_data$Ploidy, out.mat)
  }

  ## put together top track from all topdata
  ab = anno_oncoprint_barplot(border = FALSE, height = unit(track.height, "cm"))                
  toptracks = HeatmapAnnotation(column_barplot = ab)
  if (length(top_data))
  {
    topcols = brewer.master(names(top_data), wes = wes)
    tmp = lapply(names(top_data),
                 function(x) anno_barplot(top_data[[x]],
                                          border = FALSE,
                                          axis_param = list(gp = gpar(fontsize = track.fontsize)),
                                          height = unit(track.height, 'cm'),
                                          gp = gpar(fill = topcols[x], col = topcols[x])))
    names(tmp) = names(top_data)
    tmp$gap = unit(track.gap, 'cm')
    toptracks = do.call(HeatmapAnnotation, c(tmp, list(column_barplot = ab)))
  }

  packed_legends = list()
  bottomtracks = list()
  if (signature & any(oncotab$track == 'signature'))
  {
    sigd = oncotab[track == 'signature', ][type != 'Residual', ]

    ## keep any signature outside of keep that has at least signature.thresh in at least
    ## one tumor
    signature.keep = paste('Signature', signature.main, sep = '_') %>%
      union(sigd[frac>signature.thresh, type])
    sigd[, type := ifelse(type %in% signature.keep, as.character(gsub('Signature_', '', type)), 'other')]
    sigdc = dcast.data.table(sigd, id ~ type, value.var = 'frac', fun.aggregate = sum)
    sigdm = as.matrix(sigdc[, -1])
    rownames(sigdm) = sigdc$id
    sigdm = sigdm[ids,, drop = FALSE]
    sigdm = sigdm[, suppressWarnings(order(as.numeric(colnames(sigdm)))), drop = FALSE]
    out.mat = rbind(out.mat, t(sigdm))
    if (wes)
      sigcols = brewer.master(colnames(sigdm), 'BottleRocket1', wes = TRUE)
    else
      sigcols = brewer.master(colnames(sigdm), 'Dark2')

    sigcols['other'] = 'gray'
    bottomtracks$COSMIC = anno_barplot(
      sigdm,
      legend = TRUE,
      axis_param = list(gp = gpar(fontsize = track.fontsize)),
      height = unit(3*track.height, 'cm'),
      border = FALSE,
      gp = gpar(fill = sigcols, col = sigcols)
    )
    packed_legends = c(packed_legends,
      list(Legend(labels = names(sigcols), ncol = 2, legend_gp = gpar(fill = sigcols), title = 'COSMIC')))
  }

  if (svevents & any(oncotab$track %in% c('complex sv', 'simple sv')))
  {
    cx = dcast.data.table(oncotab[track == 'complex sv', ][, type := as.character(type)][, id := factor(id, ids)], id ~ type, fill = 0, drop = FALSE, value.var = 'value')
    simple = dcast.data.table(oncotab[track == 'simple sv', ][, type := as.character(type)][, id := factor(id, ids)], id ~ type, fill = 0, drop = FALSE, value.var = 'value')
    out.mat = rbind(out.mat, t(as.matrix(cx[,-1])), t(as.matrix(simple[,-1])))

    uev = names(cx)[-1]
    if (wes)
    {
      cxcols = brewer.master(names(cx)[-1], 'IsleOfDogs1', wes = TRUE)
      simplecols = brewer.master(names(simple)[-1], 'Zissou1', wes = TRUE)
    }
    else
    {
      cxcols = brewer.master(names(cx)[-1], 'Accent', wes = FALSE)
      simplecols = brewer.master(names(simple)[-1], 'Pastel1', wes = FALSE)
    }

    cxtracks = lapply(names(cx)[-1], function(x)
      anno_barplot(
        cx[[x]],
        legend = TRUE,
        axis_param = list(gp = gpar(fontsize = track.fontsize)),
        height = unit(track.height, 'cm'),
        border = FALSE,
        gp = gpar(fill = cxcols[x], col = NA)
      ))
    names(cxtracks) = names(cx)[-1]

    simpletracks = lapply(names(simple)[-1], function(x)
      anno_barplot(
        simple[[x]],
        legend = TRUE,
        axis_param = list(gp = gpar(fontsize = track.fontsize)),
        height = unit(track.height, 'cm'),
        border = FALSE,
        gp = gpar(fill = simplecols[x], col = NA)
        ))
    names(simpletracks) = names(simple)[-1]

    bottomtracks = c(bottomtracks, simpletracks, cxtracks)
  }

  ## process custom columns if any 
  if (!is.null(tumors) && length(intersect(columns, names(tumors))))
  {
    columns = intersect(columns, names(tumors))
    custom = tumors[match(ids, id), columns, with = FALSE]
    out.mat = rbind(out.mat, t(as.matrix(custom)))
    customcols = brewer.master(columns, wes = wes)
    customtracks = lapply(columns, function(x)
    {
      ## discrete data simple plot ie heatmap
      if (is.character(custom[[x]]) | is.factor(custom[[x]]) | is.logical(custom[[x]]))
      {
        if (is.logical(custom[[x]]))
          cols = c("FALSE" = 'gray', "TRUE" = 'red')
        else
          cols = brewer.master(unique(custom[[x]]), wes = wes)
        list(
          anno = anno_simple(
            as.character(custom[[x]]),
            height = unit(track.height/2, 'cm'),
            col = cols),
          legend = Legend(labels = names(cols),
                          ncol = 2, legend_gp = gpar(fill = cols, col = NA),
                          title = x)
        )
      }
      else ## numeric data barplot
        list(anno = 
               anno_barplot(
                 custom[[x]],
                 legend = TRUE,
                 axis_param = list(gp = gpar(fontsize = track.fontsize)),
                 height = unit(track.height, 'cm'),
                 border = FALSE,
                 gp = gpar(fill = customcols[x], col = NA)
               ))
    })

    customanno = lapply(customtracks, function(x) x$anno)
    names(customanno) = columns
    bottomtracks = c(bottomtracks, customanno)

    ix = lengths(customtracks)==2
    if (any(ix))
      packed_legends = c(packed_legends,
                         lapply(customtracks[ix], function(x) x$legend))
  }
  
  if (length(bottomtracks))
  {
    bottomtracks$gap = unit(track.gap, 'cm')
    bottomtracks = do.call(HeatmapAnnotation, bottomtracks)
  }

  if (length(packed_legends))
    packed_legends = do.call(packLegend, packed_legends)


  if (!is.null(split))
  {
    if (is.null(tumors))
      warning('split variable must be provided along with keyed tumors table')

    if (split %in% names(tumors))
      split = tumors[match(ids, id), ][[split]]
    else
    {
      warning('split column not found in provided tumors table')
      split = NULL
    }
  }

  gene_split = NULL
  if (!all(is.na(genes$group)))
    gene_split = genes[match(rownames(varm), genes), group]

  ## to overcome empty plot issue and also plot pct correctly
  show_pct = TRUE
  if (any(varm!='WT'))
    varm[varm == 'WT'] = ''
  else
    show_pct = FALSE ## if plot has no alterations we keep the WT so oncoPrint doesn't freak 

  if (!length(toptracks))
    toptracks = NULL

  if (!length(bottomtracks))
    bottomtracks = NULL

  op = ComplexHeatmap::oncoPrint(varm,
                      get_type = function(x) unlist(strsplit(x, ",")), ##get type = separating each cell in matrix into vector
                      alter_fun = alter_fun,
                      top_annotation = toptracks,
                      bottom_annotation = bottomtracks,
                      row_split = gene_split,
                      show_pct = show_pct, 
                      row_gap = unit(split.gap, 'cm'),
                      column_split = split,
                      column_gap = unit(split.gap, 'cm'),
                      col = varcol,
                      remove_empty_columns = FALSE,
                      remove_empty_rows = drop.genes, 
                      row_order = 1:nrow(varm),
                      column_order = 1:ncol(varm),
                      pct_gp = gpar(fontsize = rownames.fontsize),
                      row_names_gp = gpar(fontsize = rownames.fontsize),
                      column_names_gp = gpar(fontsize = colnames.fontsize),
                      show_column_names = TRUE,
                      show_heatmap_legend = TRUE
                      )


  if (ppdf)
    if (length(packed_legends))
      skitools::ppdf(draw(op, annotation_legend_list = packed_legends), height = height, width = width, ...)
    else
      skitools::ppdf(draw(op), height = height, width = width, ...)

  if (return.oncotab)
    oncotab
  else if (return.mat)
    out.mat
  else
    op
} 

#' @name alignment_metrics
#' @title alignment_metrics
#' @description
#'
#' Grabs alignment metrics from BWAMemFast like Job object (e.g. BWA mem fast)
#' or paths to bam files aligned with some variant of alignment pipelines run with
#' Picard alignment_summary_metrics postprocessing. 
#' and returns the category = 'PAIR' rows of alignment_summary_metrics indexed by entity id
#'
#' @author Marcin Imielinski
#' @export
alignment_metrics = function(job)
{
  if (is(job, 'Job'))
    fn = dirs(job, 'alignment_summary_metrics$')
  else ## assume this is  a path
  {
    if (is.null(names(job)))
      names(job) = 1:length(job)
    job = job[file.exists(job)]
    fn = lapply(dirname(job), dir, full = TRUE, pattern = 'alignment_summary_metrics')
    names(fn) = names(job)
  }
  if (!length(fn))
    return(data.table())
  res = (lapply(names(fn), function(x) if (length(fn[[x]]>0)) (grep('#', readLines(fn[[x]][1]), invert = TRUE, value = TRUE) %>% fread(text = .))[, sample := x]) %>% rbindlist(fill = TRUE))[CATEGORY == 'PAIR', ]
  return(res)
}



#' @name edge2tip
#' @title edge2tip
#'
#' Returns matrix or data.table mapping edge to tips from ape tree
#'
edge2tip = function(tree, matrix = TRUE)
{
  dt = dunlist(phangorn::Descendants(tree, tree$edge[,2], type = 'tips'))[, tip := factor(tree$tip.label[V1], tree$tip.label)][ , value := 1][, .(edge = listid, tip)] %>% setkey(edge)
  if (matrix)
  {
    dt$value = 1
    tmpdt = dcast.data.table(dt, edge ~ tip, value.var = 'value', fill = 0)
    return(as.matrix(tmpdt[,-1]>0)[, rownames(this.mat)])
  }
  return(dt)
}



#' @name circos
#' @title circos
#'
#' Quick utility function for circos plot with read depth, junctions, and segments
#' 
#' @param junctions Junction object with optional metadata field  $col to specify color
#' @param cov GRanges of scatter points with optional fields $col
#' @param segs GRanges of segments with optional fields $col and $border
#' @param win GRanges window to limit plot to
#' @param cytoband GRanges of cytoband
#' @param y.field field in cov that specifies the y axis to draw
#' @param cex.points cex for cov points
#' @param max.ranges max ranges for cov points (1e4)
#' @param ylim ylim on cov (default automatically computed)
#' @param cytoband.path path to UCSC style cytoband path
#' @param y.quantile quantile normalization
#' @param chr.sum whether to chr.sub everything 
#' @author Marcin Imielinski
#' @export

#' @name circos
#' @title circos
#'
#' Quick utility function for circos plot with read depth, junctions, and segments
#' 
#' @param junctions Junction object with optional metadata field  $col to specify color
#' @param cov GRanges of scatter points with optional fields $col
#' @param segs GRanges of segments with optional fields $col and $border
#' @param win GRanges window to limit plot to
#' @param cytoband GRanges of cytoband
#' @param y.field field in cov that specifies the y axis to draw
#' @param cex.points cex for cov points
#' @param max.ranges max ranges for cov points (1e4)
#' @param ylim ylim on cov (default automatically computed)
#' @param cytoband.path path to UCSC style cytoband path
#' @param y.quantile quantile normalization
#' @param chr.sum whether to chr.sub everything 
#' @author Marcin Imielinski
#' @export
circos = function(junctions = jJ(), cov = NULL, segs = NULL, win = NULL, field = 'ratio', cytoband = NULL, y.field = field, ylim = NA, cytoband.path = '~/DB/UCSC/hg19.cytoband.txt', cex.points = 1, ideogram.outer = TRUE, scatter = TRUE, bar = FALSE, line = FALSE, gap.after = 1, labels.cex = 1, y.quantile = 0.9999, chr.sub = TRUE, max.ranges = 1e4, axis.frac = 0.02, palette = 'BrBg', ...)
{

  if (!file.exists(cytoband.path))
    stop('cytoband not file, must be UCSC style tsv')

  if (is.null(cytoband))
    cytoband = circlize::read.cytoband(cytoband.path)$df

  cytoband = as.data.table(cytoband)
  setnames(cytoband, c('seqnames', 'start', 'end', 'band', 'stain'))

  if (chr.sub)
    cytoband[, seqnames := gsub('chr', '', seqnames)]
  
  if (!is.null(win))
  {
    if (is.character(win) | is.integer(win) | is.numeric(win) | is.factor(win))
      win = parse.gr(as.character(win))

    if (inherits(win, 'data.frame'))
      win = dt2gr(win)

    cytoband  = as.data.table(dt2gr(cytoband) %*% win)[, .(seqnames, start, end, band, stain)]
  }

  total.width = cytoband[, sum(as.numeric(end-start))]
  if (!is.na(axis.frac) && axis.frac>0)
  {
     axis.width = ceiling(axis.frac*total.width)
     cytoband = rbind(cytoband, data.table(seqnames = 'axis', start = 0, end = axis.width, band = '', stain = ''), fill = TRUE)
  }

  if (chr.sub)
  {
    ix = ((junctions$left %>% gr.sub('chr', ''))  %^% dt2gr(cytoband)) &
                          ((junctions$right %>% gr.sub('chr', '')) %^% dt2gr(cytoband))
    junctions = junctions[ix]
  }
  else
  {
    ix = junctions$left %^% dt2gr(cytoband) & junctions$right %^% dt2gr(cytoband)
    junctions = junctions[ix]
  }

  cytoband[, seqnames := as.character(seqnames)]
  args  = list(...)
  ## some important pars
  labels.cex = ifelse(is.null(args$labels.cex), 1, args$labels.cex)
  bands.height = ifelse(is.null(args$bands.height), 0.1, args$bands.height)
  cn.height = ifelse(is.null(args$cn.height), 0.3, args$cn.height)
  link.h.ratio = ifelse(is.null(args$link.h.ratio), 0.75, args$link.h.ratio)
  bpdt = junctions$dt
  bp1 = junctions$left %>% gr2dt
  bp2 = junctions$right%>% gr2dt
  circlize::circos.clear()
  circlize::circos.par(start.degree = 90, gap.after = gap.after*1)
  circlize::circos.genomicInitialize(cytoband, sector.names = unique(cytoband$seqnames), plotType = NULL, 
                                          track.height = bands.height,
                                          labels.cex = labels.cex)

  circlize::circos.genomicTrackPlotRegion(cytoband, stack = TRUE,
                                panel.fun = function(region, value, ...) {
                                  xlim = circlize::get.cell.meta.data("xlim")
                                  ylim = circlize::get.cell.meta.data("ylim")
                                  chr = circlize::get.cell.meta.data("sector.index") %>% gsub('chr', '', .)
                                  if (circlize::get.cell.meta.data("sector.index") != 'axis')
                                  {
                                    circlize::circos.text(mean(xlim), 0.9, chr, cex = 1.5, facing = "clockwise", adj = c(0,1),
                                                  niceFacing = TRUE)
                                  }
                                  }, track.height = 0.1, bg.border = NA)

  ## inner ideogram
  if (ideogram.outer)
    {
      circlize::circos.genomicTrackPlotRegion(cytoband, stack = TRUE,
                                    panel.fun = function(region, value, ...) {
                                      xlim = circlize::get.cell.meta.data("xlim")
                                      ylim = circlize::get.cell.meta.data("ylim")
                                      chr = circlize::get.cell.meta.data("sector.index")
                                      if (circlize::get.cell.meta.data("sector.index") != 'axis')
                                      {
                                        at = pretty(xlim, n = 3)
                                        circlize::circos.axis(direction = "outside", labels.facing = "outside", major.at = at, minor.ticks = 10, labels = (at/1e6) %>% as.integer, labels.cex = labels.cex*0.3)
                                        circlize::circos.genomicRect(region, value, col =  circlize::cytoband.col(value[[2]]), border = NA)
                                        circlize::circos.rect(xlim[1], ylim[1], xlim[2], ylim[2], border = "black")
                                      }
                                    }, track.height = 0.05, bg.border = NA)
    }
      
  ## coverage scatter plot
  if (!is.null(cov))
  {
    if (inherits(cov, 'data.frame'))
      cov = dt2gr(cov)

    cov = cov[!is.na(values(cov)[[y.field]])]
    cov = cov[!is.infinite(values(cov)[[y.field]])]

    if (is.na(ylim))
      ylim = c(0, quantile(values(cov)[[y.field]], y.quantile, na.rm = TRUE))
    
    cov$y = values(cov)[[y.field]] %>% as.numeric
    cov$y = cov$y %>% pmin(ylim[2]) %>% pmax(ylim[1])

    if (is.null(cov$col))
      cov$col = 'black'

    cov = cov[sample(length(cov), pmin(length(cov), max.ranges))]
    uchr = unique(cytoband$seqnames)
    cov = cov %&% dt2gr(cytoband)
    covdt = gr2dt(cov)[, seqnames := factor(seqnames, uchr)]
    circlize::circos.genomicTrackPlotRegion(covdt[, .(seqnames, start, end, y, as.character(col), ytop = y)],
                                  ylim = ylim,
                                  track.height = cn.height,
                                  bg.border = ifelse(uchr == 'axis', NA, alpha('black', 0.2)),
                                  panel.fun = function(region, value, ...) {
                                    if (circlize::get.cell.meta.data("sector.index") != 'axis')
                                    {
                                      if (circlize::get.cell.meta.data("sector.index") == uchr[1])
                                        circlize::circos.yaxis(side = 'left')                                    
                                      if (scatter)
                                        circlize::circos.genomicPoints(region, value, numeric.column = 1, col = value[[2]], pch = 16, cex = cex.points, ...)
                                      if (bar)
                                        circlize::circos.genomicRect(region, value[[1]], ytop.column = 1, border = value[[2]], col = value[[2]], pch = 16, cex = cex.points, ...)
                                      if (line)
                                        circlize::circos.genomicLines(region, value[[1]], col = value[[2]], pch = 16, cex = cex.points, ...)
                                    }
                                  })
  }
  circlize::circos.par(cell.padding = c(0, 0, 0, 0))

  if (!is.null(segs))
  {
    if (inherits(segs, 'data.frame'))
      segs = dt2gr(segs)

    if (chr.sub)
      segs = segs %>% gr.sub('chr', '')

    segs = segs[segs %^% dt2gr(cytoband), ]

    segs = as.data.table(segs)
    if (is.null(segs$col))
      segs$col = 'gray'

    if (is.null(segs$border))
      segs$border = segs$col

    if (chr.sub)
      segs[, seqnames := gsub('chr', '', seqnames)]

    circlize::circos.genomicTrackPlotRegion(segs[, .(seqnames, start, end, col, border)], stack = TRUE,
                                  panel.fun = function(region, value, ...) {
                                    circlize::circos.genomicRect(region, value, col = value[[1]], border = value[[2]])
                                    xlim = circlize::get.cell.meta.data("xlim")
                                    ylim = circlize::get.cell.meta.data("ylim")
                                    chr = circlize::get.cell.meta.data("sector.index")
#                                    circlize::circos.rect(xlim[1], ylim[1], xlim[2], ylim[2], border = "black")
                                  }, track.height = 0.05, bg.border = NA)
  }

  circlize::circos.par(cell.padding = c(0, 0, 0, 0))


  ## inner ideogram
  if (!ideogram.outer)
    {
      circlize::circos.genomicTrackPlotRegion(cytoband, stack = TRUE,
                                    panel.fun = function(region, value, ...) {
                                      xlim = circlize::get.cell.meta.data("xlim")
                                      ylim = circlize::get.cell.meta.data("ylim")
                                      chr = circlize::get.cell.meta.data("sector.index")
                                      if (circlize::get.cell.meta.data("sector.index") != 'axis')
                                      {
                                        at = pretty(xlim, n = 3)
                                        circlize::circos.axis(direction = "outside", labels.facing = "outside", major.at = at, minor.ticks = 10, labels = (at/1e6) %>% as.integer, labels.cex = labels.cex*0.3)
                                        circlize::circos.genomicRect(region, value, col = circlize::cytoband.col(value[[2]]), border = NA)
                                        circlize::circos.rect(xlim[1], ylim[1], xlim[2], ylim[2], border = "black")
                                      }
                                    }, track.height = 0.05, bg.border = NA)
    }

  if (nrow(bpdt))
  {

    if (is.null(bpdt$lwd))
      bpdt$lwd = NA_integer_

    bpdt[is.na(lwd), lwd := 1]

    if (is.null(bpdt$col))
      bpdt$col = NA_character_

    bpdt[is.na(col), col := 'red']

    if (is.null(bpdt$lty))
      bpdt$lty = NA_integer_

    bpdt[is.na(lty), lty := 1]

    if (nrow(bpdt))
      bpdt$span  = cut(junctions$span, c(0, 1e6, 3e8, Inf))

    spmap = structure(c(0.05, 0.2, 1), names = levels(bpdt$span))
    ixs = split(1:nrow(bpdt), bpdt$span)
    lapply(names(ixs), function(i)
      circlize::circos.genomicLink(
                  bp1[ixs[[i]], .(seqnames, start, end)],
                  bp2[ixs[[i]], .(seqnames, start, end)],
                  h = spmap[i],
#                  rou = circlize:::get_most_inside_radius()*c(0.1, 0.5, 1)[bpdt$span[ixs[[i]]] %>% as.integer],
                  col = bpdt[ixs[[i]], ]$col,
                  lwd =  bpdt[ixs[[i]], ]$lwd,
                  lty =  bpdt[ixs[[i]], ]$lty,
                  h.ratio = link.h.ratio,
                  border=NA)
      )
  }
  circlize::circos.clear()
}



#' Plot bin and het copy number histogram as well as a contour plot for allele fractions.
#'
#' @param cov bin coverage depth
#' @param hets het_pileup
#' @param pu purity
#' @param pl ploidy
#' @param xmax maximum value for x-axis of histograms
#' @param outputdir output directory for plots
#' @param prefix a prefix to use for output plots
#' @param suffix a suffix to use for output plots
#' @param N_subsample How many entries to randomly sample from hets when generating the contour plot
#' @export
#' @author Alon Shaiber
PPplots = function(cov, hets, pu, pl, somatic_vars=NA, xmax=10, hist_breaks=1e4, outputdir='.', prefix='', suffix = '', N_subsample=1e4){
    if (prefix!= ''){prefix = paste0(prefix, '_')}
    if (suffix!= ''){suffix = paste0('_', suffix)}
    #' histogram of bin copy number
    cov$cn = rel2abs(cov, field = 'foreground', purity = pu, ploidy = pl)
    output = paste0(outputdir, '/', prefix, 'CN_hist', suffix, '.pdf')
    ppdf(
    {
      hist(cov$cn %>% pmin(xmax), hist_breaks, xlab = 'Rescaled copy number', main = 'Copy number histogram')
      abline(v = 0:xmax, lty = 3, col = 'red')
    }, output
    )
    #' histogram of het copy numbers
    hetsc = rbind(
      hets[, .(seqnames, start, end, count = alt.count.t, type = 'alt')],
      hets[, .(seqnames, start, end, count = ref.count.t, type = 'ref')]
    )
    hetsc$ncn = 1
    hetsc$cn = rel2abs(hetsc %>% dt2gr, field = 'count', purity = pu, ploidy = pl/2)
    output =  paste0(outputdir, '/', prefix, 'het_hist', suffix, '.pdf')
    ppdf(
    {
      hist(hetsc$cn %>% pmin(xmax), hist_breaks, xlab = 'Rescaled copy number', main = 'Copy number histogram for alleles')
      abline(v = 0:xmax, lty = 3, col = 'red')
    }, output)
    #' density plots for hets
    hets2 = dcast.data.table(hetsc, seqnames + start + end ~ type, value.var = 'cn')
    hets2[, low := pmin(alt, ref)]
    hets2[, high := pmax(alt, ref)]
    output =  paste0(outputdir, '/', prefix, 'het_density', suffix, '.pdf')

    binwidths = c(MASS::bandwidth.nrd(hets2$low), MASS::bandwidth.nrd(hets2$high))
    if ((binwidths[1] <= 0) | (binwidths[2] <= 0)){
        print('The density of the het allele count is too dense and so a stat_density_2d plot cannot be generated.')
    }
    else{
        p = ggplot(hets2[sample(.N, N_subsample), ], aes(x = low, y = high, fill = ..level..)) +
                stat_density_2d(geom = "polygon") +
                scale_fill_distiller(palette = 4, direction = 1) +
                theme_bw(base_size = 25)
        ppdf(print(p), output)
    }
}


#' @name file.ready
#' @title file.ready
#' @description
#'
#' Checks if a file exists and whether it is empty or not.
#'
#' @details
#' Returns TRUE if an input file path is not NA, exists, and not an empty file.
#' If the path provided is NA then by default FALSE would be returned, unless dont_raise is set to TRUE
#' and then an error would be raised.
#' @param path path to the file
#' @param dont_raise if set to FALSE then an error would be raised if there was no path provided
#' @export
#' @author Alon Shaiber
file.ready = function(path, dont_raise=TRUE){
    not_nas = !is.na(path)
    if (any(not_nas)){
        if (dont_raise == FALSE){
            stop('There was no file provided.')
        }
    }
    if (class(path) != "character"){
        stop(sprintf('File name must be of type "character", but you provided a "%s" object.', class(path)))
    }
    exists = file.exists(path)
    nonzero = file.size(path) > 0
    return(not_nas & exists & nonzero)
}

