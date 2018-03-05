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
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        {
                          if (class(x) == 'Hits')
                            c(S4Vectors::queryLength(x), S4Vectors::subjectLength(x))
                          else
                            as.numeric(dim(x))[1:2]
                        }))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.dim)
    names(out) <- c("Type", "Size", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
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
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
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
      setnames(tab, dedup(gsub(pattern, sub, names(tab), perl = TRUE), suffix = '.'))
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
#' @param label character vector, optional specifying which data points to label (obs vector has to be named, for this to work)
#' @param plotly toggles between creating a pdf (FALSE) or an interactive html widget (TRUE)
#' @param annotations named list of vectors containing information to present as hover text (html widget), must be in same order as obs input 
#' @param gradient named list that contains one vector that color codes points based on value, must bein same order as obs input 
#' @param titleText title for plotly (html) graph only
#' @author Marcin Imielinski, Eran Hodis, Zoran Gajic
#' @export
qq_pval = function(obs, highlight = c(), exp = NULL, lwd = 1, bestfit=T, col = NULL, col.bg='black', pch=18, cex=1, conf.lines=T, max=NULL, max.x = NULL, max.y = NULL, qvalues=NULL, label = NULL, plotly = FALSE, annotations = list(), gradient = list(), titleText = "", subsample = NA, ...)
{
    if(!(plotly)){
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
        if (!is.null(exp))
            exp = -log10(exp[ix1])

        ix2 = !is.infinite(obs)
        if (!is.null(exp))
            ix2 = ix2 &  !is.infinite(exp)

        obs = obs[ix2]
        col = col[ix2]

        highlight = highlight[ix2]
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

        dat = data.table(x = sort(exp), y = obs[ord], colors = colors[ord], pch = pch, cex = cex)
        if (!is.null(names(obs)))
            {
                dat$names = names(obs[ord])
                setkey(dat, names)
            }

        if (nrow(dat)>1e5) ## rough guide to subsmapling the lower p value part of the plot
            subsample = 5e4/nrow(dat)

        if (is.na(subsample[1]))
            dat[, plot(x, y, xlab = expression(Expected -log[10](italic(P))), ylab = expression(Observed -log[10](italic(P))), xlim = c(0, max.x), col = colors, ylim = c(0, max.y), pch=pch, cex=cex, bg=col.bg, ...)]
        else
            {
                subsample = pmin(pmax(0, subsample[1]), 1)
                dat[ifelse(x<=2, ifelse(runif(length(x))<subsample, TRUE, FALSE), TRUE), plot(x, y, xlab = expression(Expected -log[10](italic(P))), ylab = expression(Observed -log[10](italic(P))), xlim = c(0, max.y), col = colors, ylim = c(0, max.y), pch=pch, cex=cex, bg=col.bg, ...)]
            }

        if (!is.null(label))
            {
                if (length(label)>0)
                    if (is.null(key(dat)))
                        warning('Need to provide names to input vector to draw labels')
                    else
                        dat[list(label), text(x, y, labels=label, pos=3)];
            }

        lines(x=c(0, max(max.y, max.x)), y = c(0, max(max.x, max.y)), col = "black", lwd = lwd)

        if (!is.na(subsample))
            dat = dat[sample(nrow(dat), subsample*nrow(dat)), ]

        lambda = lm(y ~ x-1, dat)$coefficients;

        lines(x=c(0, max.x), y = c(0, lambda*max.y), col = "red", lty = 2, lwd = lwd);
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
      if (is.character(x))
              paste("c('", paste(x, sep = "", collapse = "', '"), "')", sep = "")
      else
          paste("c(", paste(x, sep = "", collapse = ", "), ")", sep = "")
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
                  dir = file.dir(jname)
                  jname = file.name(jname)
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

                        if (bootstrap)
                            {
                                ## asking across bootstrap smaples how does the intersection fluctuate
                                ## among segments contributing to the peak

                                if (!is.null(id.field))
                                    {
                                        peak.gr = seg2gr(gr2dt(peak.gr)[, list(seqnames = seqnames[1], start = min(start),
                                            eval(parse(text = paste(field, '= sum(', field, '*(end-start))/sum(end-start)'))),end = max(end)),
                                            by = eval(id.field)])
                                        names(values(peak.gr))[] = field ## not sure why I need to do this line, should be done above
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
ra_breaks = function(rafile, keep.features = T, seqlengths = hg_seqlengths(), chr.convert = T, snowman = FALSE, swap.header = NULL,  breakpointer = FALSE, seqlevels = NULL, force.bnd = FALSE, skip = NA, 
    get.loose = FALSE ## if TRUE will return a list with fields $junctions and $loose.ends
  )
  {
      if (is.character(rafile))
          {
              if (grepl('(.bedpe$)', rafile))                  
              {
                      ra.path = rafile                      
                      cols = c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'name', 'score', 'str1', 'str2')

                      ln = readLines(ra.path)
                      if (is.na(skip))
                      {
                              nh = min(c(Inf, which(!grepl('^((#)|(chrom))', ln))))-1
                              if (is.infinite(nh))
                                  nh = 1
                          }
                      else
                          nh = skip

                                       
                      if ((length(ln)-nh)==0)
                          if (get.loose)
                              return(list(junctions = GRangesList(GRanges(seqlengths = seqlengths))[c()], loose.ends = GRanges(seqlengths = seqlengths)))
                          else                              
                              return(GRangesList(GRanges(seqlengths = seqlengths))[c()])
#                          return(GRangesList())
                      
                          
                      if (nh ==0)
                          rafile = fread(rafile, header = FALSE)
                      else
                          {
                              
                              rafile = tryCatch(fread(ra.path, header = FALSE, skip = nh), error = function(e) NULL)
                              if (is.null(rafile))
                                  rafile = tryCatch(fread(ra.path, header = FALSE, skip = nh, sep = '\t'), error = function(e) NULL)

                              if (is.null(rafile))
                                  rafile = tryCatch(fread(ra.path, header = FALSE, skip = nh, sep = ','), error = function(e) NULL)

                              if (is.null(rafile))
                                  stop('Error reading bedpe')                                  
                          }
                      setnames(rafile, 1:length(cols), cols)
                      rafile[, str1 := ifelse(str1 %in% c('+', '-'), str1, '*')]
                      rafile[, str2 := ifelse(str2 %in% c('+', '-'), str2, '*')]
#                      rafile[, str1 := ifelse(str1=='+', '-', '+')]
                                        #                      rafile[, str2 := ifelse(str2=='+', '-', '+')]
                      
                  }
              else if (grepl('(vcf$)|(vcf.gz$)', rafile))
                  {
                      library(VariantAnnotation)
                      
                      vcf = readVcf(rafile, Seqinfo(seqnames = names(seqlengths), seqlengths = seqlengths))
                      if (!('SVTYPE' %in% names(info(vcf)))) {
                        warning('Vcf not in proper format.  Is this a rearrangement vcf?')
                          return(GRangesList());
                    }
                      
                      ## vgr = rowData(vcf) ## parse BND format                      
                      vgr = read_vcf(rafile, swap.header = swap.header)
                      
                      ## no events
                      if (length(vgr) == 0)
                        return (GRangesList())

                      ## fix mateids if not included
                      if (!"MATEID"%in%colnames(mcols(vgr))) {
                        nm <- vgr$MATEID <- names(vgr)
                        ix <- grepl("1$",nm)
                        vgr$MATEID[ix] = gsub("(.*?)(1)$", "\\12", nm[ix])
                        vgr$MATEID[!ix] = gsub("(.*?)(2)$", "\\11", nm[!ix])
                        vgr$SVTYPE="BND"
                      }
                      
                      if (!any(c("MATEID", "SVTYPE") %in% colnames(mcols(vgr))))
                        stop("MATEID or SVTYPE not included. Required")
                      
                      vgr$mateid = vgr$MATEID
                      if (is.null(vgr$SVTYPE))
                          vgr$svtype = vgr$SVTYPE
                      else
                          vgr$svtype = vgr$SVTYPE

                      if (!is.null(info(vcf)$SCTG))
                          vgr$SCTG = info(vcf)$SCTG

                      if (force.bnd)
                          vgr$svtype = "BND"
                      
                      if (sum(vgr$svtype == 'BND')==0)
                          warning('Vcf not in proper format.  Will treat rearrangements as if in BND format')

                      if (!all(vgr$svtype == 'BND'))
                          warning(sprintf('%s rows of vcf do not have svtype BND, ignoring these', sum(vgr$svtype != 'BND')))

                      bix = which(vgr$svtype == "BND")
                      vgr = vgr[bix]
                      alt <- sapply(vgr$ALT, function(x) x[1])
                      vgr$first = !grepl('^(\\]|\\[)', alt) ## ? is this row the "first breakend" in the ALT string (i.e. does the ALT string not begin with a bracket)
                      vgr$right = grepl('\\[', alt) ## ? are the (sharp ends) of the brackets facing right or left
                      vgr$coord = as.character(paste(seqnames(vgr), ':', start(vgr), sep = ''))
                      vgr$mcoord = as.character(gsub('.*(\\[|\\])(.*\\:.*)(\\[|\\]).*', '\\2', alt))
                      vgr$mcoord = gsub('chr', '', vgr$mcoord)

                      if (all(is.na(vgr$mateid)))
                          if (!is.null(names(vgr)) & !any(duplicated(names(vgr))))
                              {
                                  warning('MATEID tag missing, guessing BND partner by parsing names of vgr')
                                  vgr$mateid = paste(gsub('::\\d$', '', names(vgr)), (sapply(strsplit(names(vgr), '\\:\\:'), function(x) as.numeric(x[length(x)])))%%2 + 1, sep = '::')
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
                          else
                              stop('MATEID tag missing')

                      vgr$mix = as.numeric(match(vgr$mateid, names(vgr)))

                      pix = which(!is.na(vgr$mix))

                      vgr.pair = vgr[pix]

                      if (length(vgr.pair)==0)
                          stop('No mates found despite nonzero number of BND rows in VCF')
                      vgr.pair$mix = match(vgr.pair$mix, pix)
                      vix = which(1:length(vgr.pair)<vgr.pair$mix )
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
                      ra = grl.pivot(GRangesList(vgr.pair1[, c()], vgr.pair2[, c()]))

                      this.inf = values(vgr)[bix[pix[vix]], ]

                      if (is.null(this.inf$POS))
                          this.inf = cbind(data.frame(POS = ''), this.inf)
                      if (is.null(this.inf$CHROM))
                          this.inf = cbind(data.frame(CHROM = ''), this.inf)

                      if (is.null(this.inf$MATL))
                          this.inf = cbind(data.frame(MALT = ''), this.inf)
                      
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
                      
                      if (is.null(values(ra)$TIER))
                          values(ra)$tier = ifelse(values(ra)$FILTER == "PASS", 2, 3) ## baseline tiering of PASS vs non PASS variants
                      else
                          values(ra)$tier = values(ra)$TIER

                      if (!get.loose)
                          return(ra)
                      else
                          {
                              npix = is.na(vgr$mix)
                              vgr.loose = vgr[npix, c()] ## these are possible "loose ends" that we will add to the segmentation

                              ## NOT SURE WHY BROKEN                              
                              tmp =  tryCatch( values(vgr)[bix[npix], ],
                                  error = function(e) NULL)                             
                              if (!is.null(tmp))
                                  values(vgr.loose) = tmp
                              else
                                  values(vgr.loose) = cbind(vcf@fixed[bix[npix], ], info(vcf)[bix[npix], ])
                              
                              return(list(junctions = ra, loose.ends = vgr.loose))
                          }
                  }
              else
                  rafile = read.delim(rafile)
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
  
    if (snowman) ## flip breaks so that they are pointing away from junction
      {
        rafile$str1 = ifelse(rafile$strand1 == '+', '-', '+')
        rafile$str2 = ifelse(rafile$strand2 == '+', '-', '+')
      }
      
    if (!is.null(seqlevels)) ## convert seqlevels from notation in tab delim file to actual
      {
        rafile$chr1 = seqlevels[rafile$chr1]
        rafile$chr2 = seqlevels[rafile$chr2]        
      }

     
    if (is.null(rafile$str1))
      rafile$str1 = rafile$strand1

    if (is.null(rafile$str2))
      rafile$str2 = rafile$strand2
     if (!is.null(rafile$pos1) & !is.null(rafile$pos2))
         {
             if (breakpointer)
                 {
                     rafile$pos1 = rafile$T_BPpos1
                     rafile$pos2 = rafile$T_BPpos2
                 }
             
             if (!is.numeric(rafile$pos1))
                 rafile$pos1 = as.numeric(rafile$pos1)

             if (!is.numeric(rafile$pos2))
                 rafile$pos2 = as.numeric(rafile$pos2)

             ## clean the parenthesis from the string

             rafile$str1 <- gsub('[()]', '', rafile$str1)
             rafile$str2 <- gsub('[()]', '', rafile$str2)

             ## goal is to make the ends point <away> from the junction where - is left and + is right
             if (is.character(rafile$str1) | is.factor(rafile$str1))
                 rafile$str1 = gsub('0', '-', gsub('1', '+', gsub('\\-', '1', gsub('\\+', '0', rafile$str1))))
             
             if (is.character(rafile$str2) | is.factor(rafile$str2))
                 rafile$str2 = gsub('0', '-', gsub('1', '+', gsub('\\-', '1', gsub('\\+', '0', rafile$str2))))

             
             if (is.numeric(rafile$str1))
                 rafile$str1 = ifelse(rafile$str1>0, '+', '-')

             if (is.numeric(rafile$str2))
                 rafile$str2 = ifelse(rafile$str2>0, '+', '-')
             
             rafile$rowid = 1:nrow(rafile)

             bad.ix = is.na(rafile$chr1) | is.na(rafile$chr2) | is.na(rafile$pos1) | is.na(rafile$pos2) | is.na(rafile$str1) | is.na(rafile$str2) | rafile$str1 == '*'| rafile$str2 == '*' | rafile$pos1<0 | rafile$pos2<0
             
             rafile = rafile[which(!bad.ix), ]
             
             if (nrow(rafile)==0)
                 return(GRanges())
             
             seg = rbind(data.frame(chr = rafile$chr1, pos1 = rafile$pos1, pos2 = rafile$pos1, strand = rafile$str1, ra.index = rafile$rowid, ra.which = 1, stringsAsFactors = F),
                 data.frame(chr = rafile$chr2, pos1 = rafile$pos2, pos2 = rafile$pos2, strand = rafile$str2, ra.index = rafile$rowid, ra.which = 2, stringsAsFactors = F))

             if (chr.convert)
                 seg$chr = gsub('chr', '', gsub('25', 'M', gsub('24', 'Y', gsub('23', 'X', seg$chr))))
             
             out = seg2gr(seg, seqlengths = seqlengths)[, c('ra.index', 'ra.which')];
             out = split(out, out$ra.index)
         }
     else if (!is.null(rafile$start1) & !is.null(rafile$start2) & !is.null(rafile$end1) & !is.null(rafile$end2))
         {
             ra1 = gr.flipstrand(GRanges(rafile$chr1, IRanges(rafile$start1, rafile$end1), strand = rafile$str1))
             ra2 = gr.flipstrand(GRanges(rafile$chr2, IRanges(rafile$start2, rafile$end2), strand = rafile$str2))
             out = grl.pivot(GRangesList(ra1, ra2))             
         }
     
     
     if (keep.features)
         values(out) = rafile[, ]

      if (!get.loose)
          return(out)
      else                              
          return(list(junctions = out, loose.ends = GRanges()))

     return(out)
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
#' @return scores of walks or (if raw == tRUE) raw barcode to walk maps
#' @export
#' @author Marcin Imielinski
score.walks = function(wks, bam = NULL, reads = NULL, win = NULL, wins = NULL, rthresh = 4, thresh = 1e5, pad = 1e4, raw = FALSE, allpaths = TRUE, verbose = TRUE)
{
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
  if (!is.na(rthresh))
    {
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
  bzthresh = 1
  bthresh = reads.dt[bx.diffz>bzthresh, min(bx.diff)]
  bmean = reads.dt[, mean(bx.diff, na.rm = TRUE)]

  ## concordant read pairs
  readsc = dt2gr(reads.dt[which(!discordant), ])

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
    if (is.null(reads$bx.diff))
      reads[, bx.diff := c((start-shift(end))[-1], NA), by = .(seqnames, BX)]  
    reads[, rl := label.runs(bx.diff<thresh | is.na(bx.diff)), by = .(seqnames, BX)]
    reads[is.na(rl), rl := -(1:.N)] ## label loners
    reads[, rll := paste(seqnames, BX, rl, sep = '_')]
    reads = dt2gr(reads[, .(seqnames = seqnames[1], start = start[1], end = end[.N], BX = BX[1]), by = rll])
    return(reads)
  }

  readsc = .reads2clouds(readsc)
  #readsc = dt2gr(as.data.table(grl.reduce(split(readsc + bthresh/2, readsc$BX)))[, BX := group_name])

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
  wks.chain = spChain(wks)

  ## now we want to ask what are the read pairs that now become concordant
  ## on each lifted walk??
  ## then score per BX and walk, how many discordant pairs are made concordant post-lift 
  ## and finally note which barcodes have the maximum number of their discordant pairs 
  ## lifted onto the walk
  readsdl = wks.chain * readsd
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
  rovcl = wks.chain * readsc
  rovclb = t(wks.chain) * rovcl
  values(rovclb)$walk = seqnames(rovcl)[rovclb$query.id]
  values(rovclb)$bx.walk = paste(values(rovclb)$BX, values(rovclb)$walk)

  tmp = rep(readsc, length(wks))
  tmp$walk = rep(1:length(wks), each = length(readsc))
  leftover.ix = setdiff(1:length(tmp), gr.findoverlaps(tmp, rovclb, by = c("BX", "walk"))$query.id)
  leftovers = tmp[leftover.ix]

  ## count "leftover" clouds per BX per walk
leftovers

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
  bxwid.max = data.table(BX = names(rsc), wid = grl.eval(rsc, max(width)))
  setkey(bxwid.max, BX)
  
                                        #      bxwid.lift = as.data.table(rovcl)[, .(wid.lifted = sum(width)), keyby = .(BX, seqnames)]
  bxwid.lift[bxwid, wid.left := wid-wid.lifted, on = 'BX']

  ## assume 10 "fragments" per BX genome wide, and nclbx clouds per bx
  nclbx = length(unique(readsc$rll)) / length(unique(readsc$BX))
  gsize = sum(as.numeric(seqlengths(reads)))
  wsize = sum(as.numeric(width(wins)))
  lambda = nclbx*10*wsize/gsize ## lambda for clouds per window

  ## neg.mat = negative overlap matrix from lift
  neg.mat = sparseMatrix(as.integer(factor(bxwid.lift$BX, bxlev)), as.numeric(as.character(bxwid.lift$seqnames)), x = bxwid.lift$wid.left, dims = c(length(bxlev), length(wks)), dimnames = list(bxlev, 1:length(wks)))

  ## reduce the footprint of each BX on each walk + bthresh pad
  #  rovcl.fp = as.data.table(grl.reduce(split(rovcl + bthresh, rovcl$BX)))[, BX := group_name]
  rovcl.fp = as.data.table(.reads2clouds(rovcl, bthresh))[width>bthresh ,]
  setkeyv(rovcl.fp, c("seqnames", "BX"))  
  rovcl.fp[, gaps := start-shift(end),  by = .(seqnames, BX)]
  rovcl.fp[is.na(gaps), gaps := 0]
  rovcl.fp[, pgap := dexp(gaps, 1/bmean, log = TRUE)]

  log.sum.exp =  function(x){
    offset = max(x)
    log(sum(exp(x - offset))) + offset
  }

  ## calculate widths of (largest vs next largest) footprints per walk
  bxstats = rovcl.fp[, .(jpgap = sum(pgap)), keyby = .(BX, seqnames)]
  bxstats[, lse := log.sum.exp(jpgap), by = BX]
  bxstats[, pw := exp(jpgap-lse), by = BX]

  ## bxstats = rovcl.fp[rev(order(width)), ][, .(wid.lifted = width[1]), keyby = .(BX, seqnames)]
  ## bxstats[, wid.og := bxwid.max[.(bxstats$BX), wid]] ## compare reduced wid in lifted to max wid in non lifted
  ## bxstats[, wid.rel := wid.lifted / wid.og]
  ## bxstats[, wid.rel := (wid.rel-mean(wid.rel))/range(wid.rel), by = BX]
  ## bxstats[is.na(wid.rel), wid.rel := 0]  
  ## ## pick only the walk x barcode combos with max lifted footprint 
  ## bxstats[, max.lifted := max(wid.lifted), by = BX]
  ## #bxstats = bxstats[wid.lifted == max.lifted , ]

  if (verbose)
    message("Creating barcode x walk matrices")

  ## rovc.mat = convert bxstats to barcode x walk matrix    
  ##  rovc.mat = sparseMatrix(as.integer(factor(bxstats$BX, bxlev)), as.numeric(as.character(bxstats$seqnames)), x = bxstats$wid.rel, dims = c(length(bxlev), length(wks)), dimnames = list(bxlev, 1:length(wks)))

  rovc.mat = sparseMatrix(as.integer(factor(bxstats$BX, bxlev)), as.numeric(as.character(bxstats$seqnames)), x = bxstats$pw, dims = c(length(bxlev), length(wks)), dimnames = list(bxlev, 1:length(wks)))
  
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
  neg.mat = sweep(neg.mat, 1, apply(neg.mat, 1, min), '-')

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

  sc = provd*provc*(1-pneg)
  sc = sweep(sc, 1, rowSums(sc), '/')

  browser()
  ## NA all rows that are equivalently distributed across all walks
  sc[apply(sc, 1, function(x) all(diff(x)==0)), ] = NA

  if (!is.null(wk.nm))
    colnames(sc) = wk.nm

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

      if (nchar(file.dir(file))==0)
          file = paste0('./', file)      
      
      if (!file.exists(file.dir(file)))
          system(paste('mkdir -p', file.dir(file)))

      file = paste(normalizePath(file.dir(file)), file.name(file), sep = '/')

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

         dir.create(dirname(normalizePath(file.dir(file))), recursive=TRUE, showWarnings = FALSE)
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
      text = file.name(href)

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
    source('http://bioconductor.org/biocLite.R')
    sapply(pkg, biocLite)
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
brewer.master = function(n, palette = 'Accent')
{
    nms = NULL
    if (is.character(n))
    {
        nms = unique(n)
        n = length(nms)
    }
    
  # library(RColorBrewer)
  palettes = list(
    sequential = c('Blues'=9,'BuGn'=9, 'BuPu'=9, 'GnBu'=9, 'Greens'=9, 'Greys'=9, 'Oranges'=9, 'OrRd'=9, 'PuBu'=9, 'PuBuGn'=9, 'PuRd'=9, 'Purples'=9, 'RdPu'=9, 'Reds'=9, 'YlGn'=9, 'YlGnBu'=9, 'YlOrBr'=9, 'YlOrRd'=9),
    diverging = c('BrBG'=11, 'PiYG'=11, 'PRGn'=11, 'PuOr'=11, 'RdBu'=11, 'RdGy'=11, 'RdYlBu'=11, 'RdYlGn'=11, 'Spectral'=11),
    qualitative = c('Accent'=8, 'Dark2'=8, 'Paired'=12, 'Pastel1'=8, 'Pastel2'=8, 'Set1'=9, 'Set2'=8, 'Set3'=12)
  );

  palettes = unlist(palettes);
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

      col = c(col, RColorBrewer::brewer.pal(max(next.n, 3), names(palettes[i])))
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

    ix = which(x!=zero, arr.ind = T)

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
#' @author Marcin Imielinski9
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
          names(out) = gsub(pattern, rep, file.name(out))
      else
          names(out) = file.name(out)
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

    if (!file.exists(file.dir(filename)))
        system(paste('mkdir -p', file.dir(filename)))

    cat('rendering to', filename, '\n')
    png(filename, height = height, width = width, pointsize = 24*cex.pointsize, ...)

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

    if (!file.exists(file.dir(filename)))
        system(paste('mkdir -p', file.dir(filename)))

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
       
        if (nchar(file.dir(filename))==0)
          filename = paste0('./', filename)
        
        if (!file.exists(file.dir(filename)))
            system(paste('mkdir -p', file.dir(filename)))

        filename = paste(normalizePath(file.dir(filename)), file.name(filename), sep = '/')

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
        
        if (!file.exists(file.dir(filename)))
            system(paste('mkdir -p', file.dir(filename)))
        
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
              if (!file.exists(file.dir(new.fn)))
                  system(paste('mkdir -p', file.dir(new.fn)))

              saveRDS(fn, new.fn)

              return(new.fn)
          }

      new.fn = paste(prefix, file.name(fn), sep = '')

      DEFAULT.OUTDIR = Sys.getenv('PLOP.DIR')
      if (nchar(DEFAULT.OUTDIR)==0)
          DEFAULT.OUTDIR = normalizePath('~/public_html/')

      if (!file.exists(DEFAULT.OUTDIR))
          system(paste('mkdir -p', file.dir(filename)))

      if (!is.null(force))
          {
              if (length(force)==1)
                  force = rep(force, length(new.fn))
              new.fn = force
          }
      
      new.fn = paste(DEFAULT.OUTDIR, new.fn, sep = '/')

      if (any(ix <- !file.exists(file.dir(new.fn))))
          sapply(new.fn[ix], function(x)
              system(paste('mkdir -p', file.dir(x))))
      
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
               
                m = lm(y ~ x, dat)
                abline(m, lwd = 3, lty = 2, col = col.fit)

                eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                                 list(a = format(coef(m)[1], digits = 2),
                                      b = format(coef(m)[2], digits = 2),
                                      r2 = format(summary(m)$r.squared, digits = 3)))


                if (coef(m)[2]>0)
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
#' @param scale scale parameter to geom_vplot (=area)
#' @param log logical flag whether to log transform (=FALSE)
#' @param count logical flag whether to include counts in ylabels (=TRUE)
#' @param xlab character xlabel (=NULL)
#' @param ylab character ylabel (=NULL)
#' @param minsup minimum support to include in a group (=NA)
#' @param scatter logical flag whether to include scatter of points (=FALSE)
#' @param alpha numeric vector between 0 and 1 to specify alpha transparency of points if scatter is TRUE (0.3)
#' @param title character specifying plot title (=NULL)
#' @author Marcin Imielinski
#' @import ggplot2
#' @export
vplot = function(y, group = 'x', facet1 = NULL, facet2 = NULL, transpose = FALSE, flip = FALSE,  mapping = NULL,
    stat = "ydensity",
    position = "dodge",
    trim = TRUE, sample = NA, scale = "area", log = FALSE, count = TRUE, xlab = NULL, ylim = NULL, ylab = NULL, minsup = NA,
    scatter = FALSE,
    text = NULL,
    cex.scatter = 1,
    col.scatter = NULL, alpha = 0.3, title = NULL, legend.ncol = NULL, legend.nrow = NULL, vfilter = TRUE, vplot = TRUE, dot = FALSE, stackratio = 1, binwidth = 0.1, plotly = FALSE, print = TRUE)
    {
        # require(ggplot2)
      if (!is.factor(group))
          group = as.factor(group)
      dat = data.table(y = suppressWarnings(as.numeric(y)), group)
      
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

        ## if (vfilter)
        ##     {
        vgroup = NULL ## NOTE fix
        good = as.data.table(dat)[, list(var = var(y)), keyby = vgroup][var>0, vgroup]
        dat = dat[, vfilter := dat$vgroup %in% as.character(good)]
        ## }

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
                levels(dat$group) = paste(names(tmp)[ix], '\n(', tmp[ix], ')', sep = '')
            }

        if (is.null(mapping))
            mapping = aes(fill=group)

        g = ggplot(dat[vfilter!=0, ], aes(y = y, x = group)) + theme_bw()

        if (vplot)
            g = g + geom_violin(mapping = mapping, stat = stat, position = position, trim = trim, scale = scale)

        if (scatter)
            {
                if (dot)
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
#                    g = g + coord_trans(y = 'log10')
                else
                    g = g+ scale_y_log10(limits = ylim)
#                    g = g + coord_trans(y = 'log10', limits = ylim)
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


      if (flip)
        g = g + coord_flip()

        if (!is.null(dat$facet1))
            {
                if (!is.null(dat$facet2))
                    {
                        if (transpose)
                            g = g + facet_grid(facet2 ~ facet1)
                        else
                            g = g + facet_grid(facet1 ~ facet2)
                    }
                else
                    {
                        if (transpose)
                            g = g + facet_grid(. ~ facet1)
                        else
                            g = g + facet_grid(facet1 ~ .)
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
        x = paste('grep -H', grep, paste(x, collapse = ' '), ' | more')
    system(x)
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

  vcf.paths = dir(oncotator.dir, 'vcf', ful.names = TRUE)
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
      names(bams) = file.name(bams)

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

  out = data.frame(p = rep(NA, length(gene.sets)), fdr = NA, genes = '', stringsAsFactors = F)
  rownames(out) = names(gene.sets);

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

  for (i in 1:nrow(out))
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
          out[i, 'genes'] = paste(paste(names(obs.dat), ' (', obs.dat, ')', sep = ''), collapse = ', ')

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

          out[i, 'leading.edge'] = paste(paste(names(obs.dat[leading.edge]), ' (', obs.dat[leading.edge], ')', sep = ''), collapse = ', ')

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

          out[i, 'p'] = signif((1+sum(perm>obs))/(1+length(perm)), 2);
          out[i, 'obs'] = obs;

                                        #      if (verbose) print(out[i,  setdiff(names(out), 'genes')])
          if (verbose) print(out[i,  ])
        }
    }

  out$fdr = p.adjust(out$p, 'BH');

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
            igv.cmd(paste('snapshotDirectory', gsub('^\\~', '$HOME', file.dir(snapshot))), con)
            igv.cmd(paste('snapshot', file.name(snapshot)), con)
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


#' Minimal overlaps for GRanges/GRangesList
#'
#' Takes any number of GRanges or GRangesList and reduces them to the minimal
#' set of overlapping windows, ignoring strand (optional).  Can also
#' collapse only within levels of a meta data field "by"
#'
#' Will populate output with metadata of first row of input contributing the reduced output range.
#'
#' @param ... \code{GRanges} or \code{GRangesList}
#' @return GRanges
#' @export
gr.reduce <- function(..., by = NULL, ignore.strand = TRUE, span = FALSE) {
    input <- do.call(grbind, list(...))
    if (length(input)==0)
        return(input)
    input.meta = values(input)
    if (is.null(by))
        values(input) = data.frame(i = 1:length(input), bykey = 1)
    else
        values(input) = data.frame(i = 1:length(input), bykey = values(input)[, by])

    if (span)
    {
        if (ignore.strand)
            out = seg2gr(gr2dt(input)[, data.frame(i = i[1], start = min(start), end = max(end)), keyby = list(seqnames, bykey)])
        else
            out = seg2gr(gr2dt(input)[, data.frame(i = i[1], start = min(start), end = max(end)), keyby = list(seqnames, strand, bykey)])
    }
    else
    {
        if (ignore.strand)
            out = seg2gr(gr2dt(input)[, cbind(i = i[1], as.data.frame(reduce(IRanges(start, end)))), keyby = list(seqnames, bykey)])
        else
            out = seg2gr(gr2dt(input)[, cbind(i = i[1], as.data.frame(reduce(IRanges(start, end)))), keyby = list(seqnames, strand, bykey)])
    }

    values(out) = input.meta[out$i, ]

                                        #input = do.call(grbind, input)
    ## for (i in seq_along(input)) {
    ##   if (inherits(input[[i]], 'GRanges'))
    ##     input[[i]] <- reduce(gr.stripstrand(input[[i]]))
    ##   else if (inherits(input[[i]], 'GRangesList'))
    ##     input[[i]] <- reduce(gr.stripstrand(unlist(input[[i]])))
    ##   else
    ##     stop('reduce.window: Need to input GRanges or GRangesList objects')
    ##   seqlengths(input[[i]]) <- seqlengths(input[[i]])*NA
    ## }

    ## output <- do.call('c', input)

    return(out)
                                        #return(sort(reduce(output)))
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
#' Tabulates cluster usage (qstat()) or if full = TRUE flag given will dump out
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
#' Takes string representing time elapsed and returns numeric string representing number of GB
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
qviz = function(res = NULL, queue = NULL, field = "h_vmem", frac = FALSE, all = FALSE, plot = TRUE)
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
                txt = sprintf("detach('package:%s', force = TRUE)", lib)
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


#####################################
#' @name heatmap.plus
#' @title heatmap.plus
#' @description
#' 
#' Additional features:
#'    allows several label tracks on top, bottom, left, and right, with separate top and bottom legend frames tohouse each
#'    allows use of coloredData in tracks
#'
#' @export
###################################
heatmap.plus = function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL,
  bar = FALSE, ## if TRUE will draw barplot for central panel instead of heatmap
  show.rdend = TRUE, # show row dendogram flag
  show.cdend = TRUE,  # show column dendrogram flag
  # these next four args can be coloredData or list of coloredData with category labels
  # for each column / row data point.  If named will be indexed by corresponding row/column names in X
  # otherwise they will be indexed in order.  These args can alternatively can be vector or list of vectors (which
  # will be mapped to default colormaps)
  topColAttr = NULL, # see above
  bottomColAttr = NULL, # ...
  leftRowAttr = NULL, # ...
  rightRowAttr = NULL, # ...
  leg.args = NULL,  ## legend will be populated by color mappings from coloredData and additional args (excluding position) here
  dim.heatmap = c(4, 4), ## use this to make heatmap tall or fat (instead of margins arg)
  distfun = dist, hclustfun = hclust,
  reorderfun = function(d, w) reorder(d, w),
  add.expr,  
  symm = FALSE,
  revC = identical(Colv, "Rowv"), scale = c("row", "column", "none"), na.rm = TRUE,
  margins = c(1, 1), 
  cexRow = 0.2 +  1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL,
  add.grid = F,
  col.grid = 'gray',
  lwd.grid = 1,
  size.legend.panel = 0.4,
  size.feature.panel = 0.2,
  col = topo.colors(100),
  optimal.leaf = T,
  return.clust = F,    
  labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, las.col = 2,
  verbose = getOption("verbose"), ...) 
{
  require(cba)
  hcr = hcc = NULL

  ### set up top and bottom color palettes
  if (!is.null(topColAttr) | !is.null(bottomColAttr) | !is.null(leftRowAttr) | !is.null(rightRowAttr))
    {
      require(RColorBrewer)
      
      brewer.palettes = brewer.pal.info
      brewer.palettes = brewer.palettes[order(brewer.palettes$category), ]
  
      if (!is.null(topColAttr) & !is.list(topColAttr))
        topColAttr = list(topColAttr)
      
      if (!is.null(bottomColAttr) & !is.list(bottomColAttr))
        bottomColAttr = list(bottomColAttr)
      
      if (!is.null(leftRowAttr) & !is.list(leftRowAttr))
        leftRowAttr = list(leftRowAttr)
      
      if (!is.null(rightRowAttr) & !is.list(rightRowAttr))
        rightRowAttr = list(rightRowAttr)

      last.palette = 0;

      .convert.to.cData = function(x) {        
        if (!is.null(x) & class(x) != 'coloredData')
          {
            x = as.vector(x);
            last.palette <<- last.palette + 1;
            if (is.factor(x))
              uval = levels(x)
            else
              uval = unique(x)
            cmap = brewer.pal(min(brewer.palettes[last.palette, 'maxcolors'], length(uval)), rownames(brewer.palettes)[last.palette])

            if (length(uval)>length(cmap))
              {
                warning('Number of colors exceeded for colormap: duplicate colors will be created. ')
                cmap = cmap[((1:length(uval))%%length(cmap))+1]  ## we will repeat colors if the number of unique items is larger than the colormap
              }
            names(cmap) = uval;
            return(coloredData(data = x, colormap = cmap))
          }
        else
          return(x)
      }
      topColAttr = lapply(topColAttr, .convert.to.cData) ## if any attributes are not already colored data, transform them using coloredData using pre-applied palettes
      bottomColAttr = lapply(bottomColAttr, .convert.to.cData) ## if any attributes are not already colored data, transform them using coloredData using pre-applied palettes
      leftRowAttr = lapply(leftRowAttr, .convert.to.cData) ## if any attributes are not already colored data, transform them using coloredData using pre-applied palettes
      bottomColAttr = lapply(bottomColAttr, .convert.to.cData) ## if any attributes are not already colored data, transform them using coloredData using pre-applied palettes      
    }

  if (!is.null(colnames(x)))
    colNames = colnames(x)
  else
    colNames = 1:ncol(x)
  
  if (!is.null(rownames(x)))
    rowNames = rownames(x)
  else
    rowNames = 1:nrow(x)
   
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("'x' must be a numeric matrix")
    nr <- di[1L]
    nc <- di[2L]
    if (nr <= 1 || nc <= 1) 
        stop("'x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2L) 
        stop("'margins' must be a numeric vector of length 2")
    doRdend <- !identical(Rowv, NA)
    doCdend <- !identical(Colv, NA)
    if (!doRdend && identical(Colv, "Rowv")) 
        doCdend <- FALSE
    if (is.null(Rowv)) 
        Rowv <- rowMeans(x, na.rm = na.rm)
    if (is.null(Colv)) 
        Colv <- colMeans(x, na.rm = na.rm)
    if (doRdend) {
        if (inherits(Rowv, "dendrogram")) 
            ddr <- Rowv
        else {
            d = distfun(x);
            hcr = hclustfun(d);
            
            if (optimal.leaf) ## MARCIN ADDED
            {
              tmp <- order.optimal(d, hcr$merge)
              hcr$order = tmp$order;
              hcr$merge = tmp$merge;
              ddr <- as.dendrogram(hcr)
            }
          else
            {
              ddr <- as.dendrogram(hcr)
              if (!is.logical(Rowv) || Rowv) 
                ddr <- reorderfun(ddr, Rowv)
            }
        }
        if (nr != length(rowInd <- order.dendrogram(ddr))) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else rowInd <- 1L:nr
    if (doCdend) {
        if (inherits(Colv, "dendrogram")) 
            ddc <- Colv
        else if (identical(Colv, "Rowv")) {
            if (nr != nc) 
                stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
            ddc <- ddr
        }
        else {
          d = distfun(if (symm) x else t(x));
          hcc = hclustfun(d);
            
          if (optimal.leaf) ## MARCIN ADDED
            {
              tmp <- order.optimal(d, hcc$merge)
              hcc$order = tmp$order;
              hcc$merge = tmp$merge;
              ddc <- as.dendrogram(hcc)
            }
          else
            {
              ddc <- as.dendrogram(hcc)
              if (!is.logical(Colv) || Colv) 
                ddc <- reorderfun(ddc, Colv)
            }
        }
        if (nc != length(colInd <- order.dendrogram(ddc))) 
          stop("column dendrogram ordering gave index of wrong length")
      }
    else colInd <- 1L:nc
    x <- x[rowInd, colInd]
    labRow <- if (is.null(labRow)) 
        if (is.null(rownames(x))) 
            (1L:nr)[rowInd]
        else rownames(x)
    else labRow[rowInd]
    labCol <- if (is.null(labCol)) 
        if (is.null(colnames(x))) 
            (1L:nc)[colInd]
        else colnames(x)
    else labCol[colInd]
    if (scale == "row") {
        x <- sweep(x, 1L, rowMeans(x, na.rm = na.rm), check.margin = FALSE)
        sx <- apply(x, 1L, sd, na.rm = na.rm)
        x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
    }
    else if (scale == "column") {
        x <- sweep(x, 2L, colMeans(x, na.rm = na.rm), check.margin = FALSE)
        sx <- apply(x, 2L, sd, na.rm = na.rm)
        x <- sweep(x, 2L, sx, "/", check.margin = FALSE)
    }
    lmat <- rbind(c(NA, 3), 2:1)
    lwid <- c(if (doRdend) 1 else 0.05, dim.heatmap[1])
    lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0, 
        dim.heatmap[2])
  
  ## COL LABEL LAYOUTS  
  col.panel.height = size.feature.panel;
  row.panel.width = size.feature.panel;
  core.panel.ind = c(2,2);
  
  if (!is.null(topColAttr))  ## add to top
    lapply(topColAttr, function(x)
           {
             lmat <<- rbind(lmat[1, ]+1, c(NA, 1), lmat[2:nrow(lmat), ]+1)             
             lhei <<- c(lhei[1L], col.panel.height, lhei[2:length(lhei)])
             core.panel.ind[1] <<- core.panel.ind[1]+1
           })

  if (!is.null(bottomColAttr)) # add to bottom
    lapply(bottomColAttr, function(x)
           {
             lmat <<- rbind(lmat+1, c(NA, 1))
             lhei <<- c(lhei, col.panel.height)
           })

  ## ROW LABEL LAYOUTS
  if (!is.null(leftRowAttr))
    lapply(leftRowAttr, function(x)
           {             
             lmat <<- cbind(lmat[, 1] + 1, c(rep(NA, core.panel.ind[1]-1),  1, rep(NA, nrow(lmat)-core.panel.ind[1])), lmat[, 2:ncol(lmat)] + 1)
             lwid <<- c(lwid[1L], row.panel.width, lwid[2L])
             core.panel.ind[2] <<- core.panel.ind[2]+1;
           })

  if (!is.null(rightRowAttr))
    lapply(rightRowAttr, function(x)
           {
             lmat <<- cbind(lmat + 1, c(rep(NA, core.panel.ind[1]-1),  1, rep(NA, nrow(lmat)-core.panel.ind[1])))
             lwid <<- c(lwid, row.panel.width)
           })

  ## ADD ROW COL LEGEND PANELS

  # top legend panel
  new.row = rep(NA, ncol(lmat)); new.row[core.panel.ind[2]] = max(lmat, na.rm = T)+1;
  core.panel.ind[1] = core.panel.ind[1]+1
  lhei = c(size.legend.panel, lhei)
  lmat = rbind(new.row, lmat)
  
  # bottom legend panel
  new.row = rep(NA, ncol(lmat)); new.row[core.panel.ind[2]] = max(lmat, na.rm = T)+1;
  lmat = rbind(lmat, new.row)  
  lhei = c(lhei, size.legend.panel)
  
  # left legend panel
  new.col = rep(NA, nrow(lmat)); new.col[core.panel.ind[1]] = max(lmat, na.rm = T)+1;
  core.panel.ind[2] = core.panel.ind[2]+1
  lmat = cbind(new.col, lmat)
  lwid = c(size.legend.panel, lwid)
  
  # right legend panel
  new.col = rep(NA, nrow(lmat)); new.col[core.panel.ind[1]] = max(lmat, na.rm = T)+1;
  lmat = cbind(lmat, new.col)
  lwid = c(lwid, size.legend.panel)
 
  lmat[is.na(lmat)] <- 0
  if (verbose) {
    cat("layout: widths = ", lwid, ", heights = ", lhei, 
        "; lmat=\n")
    print(lmat)
    }
  op <- par(no.readonly = TRUE)
  on.exit(par(op))

  ## LAYOUT 
  layout(lmat, widths = lwid, heights = lhei, respect = TRUE) 

  if (verbose)
    print(lmat)
  
  pad.feature.panel = 0.1
  
  ## DRAW ROW LABEL PANELS
  if (!is.null(rowAttr <- c(leftRowAttr, rightRowAttr)))
    if (length(rowAttr)>0)
      lapply(rev(rowAttr), function(x)
             {
               if (is.null(x)) return
               par(mar = c(margins[1L]/2, pad.feature.panel, margins[1L]/2, pad.feature.panel))
               if (is.null(names(getData(x)))) ## we assume attributes are ordered
                 image(1, 1:nr, z = rbind(1L:nr), col = getColors(x)[rowInd], axes = FALSE, xlab = "", ylab = "")
               else   ## otherwise will use attribute names to properly order (assuming that data rownames are specified)r
                 image(1, 1:nr, z = rbind(1L:nr), col = getColors(x)[rowNames[rowInd]], axes = FALSE, xlab = "", ylab = "")
             })
  
  ## DRAW COLUMN LABEL PANELS 
  if (!is.null(colAttr <- c(topColAttr, bottomColAttr)))
    if (length(colAttr)>0)
      lapply(rev(colAttr), function(x)
             {
               if (is.null(x)) return
               par(mar = c(pad.feature.panel, margins[2L]/2, pad.feature.panel, margins[2L]/2))               
               if (is.null(names(getData(x)))) ## we assume attributes are ordered
                 image(1:nc, 1, z = cbind(1L:nc), col = getColors(x)[colInd], axes = FALSE, xlab = "", ylab = "")
               else ## otherwise will use attribute names to properly order (assuming that data colnames are specified)
                 image(1:nc, 1, z = cbind(1L:nc), col = getColors(x)[colNames[colInd]], axes = FALSE, xlab = "", ylab = "")
             })  
  
  par(mar = c(margins[1L]/2, margins[2L]/2, margins[1L]/2, margins[2L]/2))
  if (!symm || scale != "none") 
    x <- t(x)
  if (revC) {
    iy <- nr:1
    if (doRdend) 
      ddr <- rev(ddr)
    x <- x[, iy]
  }
  else iy <- 1L:nr
  
 ## CENTRAL PANEL (heatmap)
  xlim = 0.5 + c(0, nc);
  ylim = 0.5 + c(0, nr);
  
  if (bar)
    {
      plot.blank()
      par(usr = c(0, nc, 0, max(rowSums(x))))
      barplot(t(x), axes = FALSE, xlab = "", ylab = "", col = col[rowInd], space = 0, names.arg = rep('', nc), add = T)

      if (las.col == 2)
      axis(1, (1L:nc)-0.5, labels = labCol, las = las.col, line = -0.5, tick = 0, cex.axis = cexCol)
    else
      axis(1, (1L:nc)-0.5, labels = labCol, las = las.col, padj = 1, line = -0.5, tick = 0, cex.axis = cexCol)

    }
  else
    {
      image(1L:nc, 1L:nr, x, xlim = xlim, ylim = ylim, axes = FALSE, xlab = "", ylab = "", col = col, ...)
      
      if (las.col == 2)
      axis(1, 1L:nc, labels = labCol, las = las.col, line = -0.5, tick = 0, cex.axis = cexCol)
    else
      axis(1, 1L:nc, labels = labCol, las = las.col, padj = 1, line = -0.5, tick = 0, cex.axis = cexCol)

    }
  
  if (add.grid)
    {
       segments(-.5 + 1:(nc+1), -.5, -.5 + 1:(nc+1), .5 + nr, col= col.grid, lwd=lwd.grid)
       segments(-.5, -.5 + 1:(nr+1), .5 + nc, -.5 + 1:(nr+1), col= col.grid, lwd= lwd.grid)          
     }
  
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1L] - 1.25)

  
   axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, cex.axis = cexRow)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2L] - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))

  ## LEFT DENDROGRAM
  par(mar = c(margins[1L]/2, 0, margins[1L]/2, 0))
    if (doRdend & show.rdend)
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    else
      if (!is.null(rownames(x)))
        {
          plot.blank(ylim = ylim, xlim = c(0, 1))
          text(rep(0.96, length(labRow)), (1:nr), labRow, srt = 0, adj = c(1, 0.5), cex = cexRow)
        }
      else
        frame()
    par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2L]))

  ## TOP DENDROGRAM
  par(mar = c(0, margins[2L]/2, 0, margins[2L]/2))
    if (doCdend & show.cdend) 
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    else
      if (!is.null(colnames(x)))
        {
          plot.blank(xlim = xlim, ylim = c(0, 1))
          par(usr = c(xlim, 0, 1))
          text((1:nc), rep(0.04, length(labCol)), labCol, srt = 90, adj = c(0, 0.5), cex = cexCol)
        }
      else
        frame()
  
    if (!is.null(main)) {
        par(xpd = NA)
        title(main, cex.main = 1.5 * op[["cex.main"]])
    }

  ### LEGENDS

  # set defaults
  if (is.null(leg.args$y)) leg.args$y = 0.5;
  if (is.null(leg.args$x)) leg.args$x = 0.5;
  if (is.null(leg.args$adj)) leg.args$adj = c(0, 0.5);
  if (is.null(leg.args$xjust)) leg.args$xjust = 0.5;
  if (is.null(leg.args$yjust)) leg.args$yjust = 0.5;
  
  ## TOP
  par(mar = c(1,0,0,0))
  par(xpd = NA);
  plot(0:1,0:1, type = "n", axes = FALSE, xlab = "", ylab = "")
  if (length(topColAttr)>0)
    {
      these.leg.args = lapply(1:length(topColAttr),  function(x) {
        out = leg.args;
        out$x = x/(length(topColAttr)+1)
        out$legend = names(getColormap(topColAttr[[x]]))
        out$fill = getColormap(topColAttr[[x]])
        return(out)
      })
      sapply(these.leg.args, function(x) do.call('legend', x)) ## call all the legends
    }

  ## BOTTOM
  par(mar = c(0,0,1,0))
  plot(0:1,0:1, type = "n", axes = FALSE, xlab = "", ylab = "")
  if (length(bottomColAttr)>0)
    {
      these.leg.args = lapply(1:length(bottomColAttr),  function(x) {
        out = leg.args;
        out$x = x/(length(bottomColAttr)+1)
        out$legend = names(getColormap(bottomColAttr[[x]]))
        out$fill = getColormap(bottomColAttr[[x]])
        return(out)
      })
      sapply(these.leg.args, function(x) do.call('legend', x)) ## call all the legends
    }                       

  ## LEFT
  par(mar = c(0,0,0,1))
  plot(0:1,0:1, type = "n", axes = FALSE, xlab = "", ylab = "")
  if (length(leftRowAttr)>0)
    {
      these.leg.args = lapply(1:length(leftRowAttr),  function(x) {
        out = leg.args;
        out$y = x/(length(leftRowAttr)+1)
        out$legend = names(getColormap(leftRowAttr[[x]]))
        out$fill = getColormap(leftRowAttr[[x]])
        return(out)
      })
      sapply(these.leg.args, function(x) do.call('legend', x)) ## call all the legends
    }

  ## RIGHT
  par(mar = c(0,1,0,0))
  plot(0:1,0:1, type = "n", axes = FALSE, xlab = "", ylab = "")
  if (length(rightRowAttr)>0)
    {
      these.leg.args = lapply(1:length(rightRowAttr),  function(x) {
        out = leg.args;
        out$y = x/(length(rightRowAttr)+1)
        out$legend = names(getColormap(rightRowAttr[[x]]))
        out$fill = getColormap(rightRowAttr[[x]])
        return(out)
      })
      sapply(these.leg.args, function(x) do.call('legend', x)) ## call all the legends
    }                       
  
  invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro && 
                                                     doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
  
  if (return.clust)
    return(list(row = hcr, col = hcc))
}


#' @name coloredData
#' @rdname coloredData
#' 
#' @description
#' S4 class for data with colors used by heatmap.plus
#'
#' simple object with data (e.g. vector or matrix of categorical, real numbers) + a colormap
#'
#'  colormap is a (named) vector mapping factor levels / unique values in data into colors,
#'  or otherwise assigning a color range to numeric data.
#'
#' @exportClass coloredData
#' @author Marcin Imielinski
setClass('coloredData', representation(data = 'array', colormap = 'vector', type = 'character', data.names = 'character'),
         prototype(data = matrix(NA), colormap = NA, type = '', data.names = NULL)
         )
setMethod('initialize', 'coloredData', function(.Object, data, colormap, upright = T)
         {
           .Object = callNextMethod()
           .Object@type = class(data)

           if (!is.vector(colormap) | !any(!is.na(colormap)))
             stop('colormap must be a vector with non NA entries')
           
           .Object@colormap = colormap
           if (is.vector(data) & !is.list(data))
             {
               if (!is.null(names(data)))
                 .Object@data.names = names(data)
               if (upright)               
                 data = array(data, dim = c(length(data), 1), dimnames = list(names(data), NULL))
               else
                 data = array(data, dim = c(1, length(data)), dimnames = list(NULL, names(data)))
             }
           else if (is.matrix(data))
               data = array(as.vector(data), dim = dim(data), dimnames = dimnames(data))
           else if (!is.array(data))
             stop('Only vectors, matrices, and arrays supported')
           .Object@data = data
           if (!is.numeric(.Object@data)) {
             if (is.factor(.Object@data))
               uval = levels(.Object@data)
             else if (is.character(.Object@data))
               uval = unique(.Object@data)
             else
               uval = NULL;
             
             if (!is.null(uval))
               {             
                 if (is.null(names(.Object@colormap)))
                   {
                     ix = 1:min(length(uval), length(.Object@colormap))              
                     names(.Object@colormap)[ix] = uval[ix]
                   }
                 if (length(leftover <- setdiff(uval, c(NA, names(.Object@colormap))))>0)
                   warning(sprintf('The following factors are unmapped in the colormap: %s.', paste(leftover, collapse = ",")))
               }
           } 
           .Object
       })


setGeneric('getColormap', function(.Object) standardGeneric('getColormap'))
setGeneric('getData', function(.Object) standardGeneric('getData'))
setGeneric('getColors', function(.Object) standardGeneric('getColors'))
setMethod('getColormap', signature('coloredData'), function(.Object) .Object@colormap)
setMethod('getColors', signature('coloredData'), function(.Object)
          {
            dat = getData(.Object);
            cmap = getColormap(.Object);
            dat[1:length(dat)] = cmap[dat];
            return(dat)
          })
setMethod('getData', signature(.Object = 'coloredData'), function(.Object)
          {
            out = as(.Object@data, .Object@type)
            if (!is.null(.Object@data.names))
              names(out) = .Object@data.names
            return(out)
        })

#' @name coloredData
#' @title coloredData
#' @description
#'
#' Instantiate coloredData
#' 
#' @export
coloredData = function(data, colormap, upright = T) new('coloredData', data = data, colormap = colormap, upright = upright)


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
#' @param camera.res output of camera from limma
#' @param voom.res output of voom from limma
#' @param gene.sets gene set input to camera (named list of indices into the voom.res gene expression matrix)
#' @param design design matrix input to camera
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
cameraplot = function(camera.res, gene.sets, voom.res, design, contrast = ncol(design), 
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
        
        camera.res$name = rownames(camera.res)
        setnames(camera.res, gsub('\\W', '_', names(camera.res), perl = TRUE))
        camera.res = as.data.table(camera.res)[order(-PValue), ]
        my.sets = gene.sets[camera.res$name]
        my.corr = apply(voom.res$E, 1, cor, y = design[,contrast])
        my.rank = rank(-my.corr)
        my.set.rank = sapply(my.sets, function(x) my.rank[x])
        cm = function(x, width = 0.5) rgb(colorRamp(col.ramp)((pmax(-width, pmin(width, x))+width)/(2*width)), maxColorValue = 256)
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
        plot(0, type ="n", xlim = c(0, length(my.rank)), ylim = c(-0.5,0.5), ylab = '', xlab = "", axes = F, main = title)
        par(mar = 0.5*c(5,20, 0, 5))
        axis(2, at = seq(-0.5, 0.5, 0.5), col.axis = col.axis)
        mtext(side = 2, 'Gene Correlations', line = 3, col = col.axis)
        lines(my.rank, my.corr, type = 'h', col = cm(my.corr))
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
             col = cm(gene.corr), cex = cex.glab*0.3, adj = c(0.5, 0))
                                        # draw lines linking gene labels to notches
        segments(gene.coord.top[,1],
                 gene.coord.top[,2],
                 gene.coord.bot[,1],
                 gene.coord.bot[,2]+gtext.shift*0,
                 col = alpha(cm(gene.corr), 0.2), lty = 1)
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
                 notch.coord[,2]+tick.w/2, col = cm(my.corr[rownames(notch.coord)]), lwd = 2*lwd.notch)

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
            out$vartype = ifelse(nchar(out$REF) == nchar(out$ALT), 'SNV',
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

#' @name label.runs
#' @title label.runs
#' @description
#'
#' Labels runs of TRUE in a vector with increasing indices
#' @author Marcin Imielinski
#' @export
label.runs = function(x) as.integer(ifelse(x, cumsum(diff(as.numeric(c(FALSE, x)))>0), NA))    

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
  else if ('lm' %in% class(x))
  {

    coef = as.data.frame(summary(x)$coefficients)
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
    if (!last)
      nm = paste(nm, rownames(coef))
    out = data.frame(name = nm, method = summary(x)$family$family, p = signif(coef$p, 3), estimate = coef$estimate, ci.lower = coef$ci.lower, ci.upper = coef$ci.upper, effect = paste(signif(coef$estimate, 3), ' [',  signif(coef$ci.lower,3),'-', signif(coef$ci.upper, 3), ']', sep = ''))
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
      out = suppressWarnings(read_vcf(x))
      if (is.na(label))
        label = x
    }
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
            vcf$eff = unstrsplit(vcf$ANN)
            vcf$ref = as.character(vcf$REF)
            vcf$alt = as.character(unstrsplit(vcf$ALT))
            vcf$modifier = !grepl('(HIGH)|(LOW)|(MODERATE)', vcf$eff)
            if (!keep.modifier)
              vcf = vcf[!vcf$modifier]
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
        tmp = lapply(out$ANN, function(y) do.call(rbind, strsplit(y, '\\|'))[, 1:15, drop = FALSE])
        tmpix = rep(1:length(out), sapply(tmp, nrow))
        meta = as.data.frame(do.call(rbind, tmp))
        colnames(meta) = fn
        meta$varid = tmpix
        meta$file = x
        out2 = out[tmpix]
        rownames(meta) = NULL
        values(out2) = cbind(values(out2), meta)
        names(out2) = NULL
        out2$ANN = NULL
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
#' @descriptino
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
    kag = karyograph(junctions)
    
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
        tmp = spChain(walks)$y;
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
#' @param ...
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



#' @name ra.overlaps
#' @title ra.overlaps
#' @description
#'
#' Determines overlaps between two piles of rearrangement junctions ra1 and ra2 (each GRangesLists of signed locus pairs)
#' against each other, returning a sparseMatrix that is T at entry ij if junction i overlaps junction j.
#'
#' if argument pad = 0 (default) then only perfect overlap will validate, otherwise if pad>0 is given, then
#' padded overlap is allowed
#'
#' strand matters, though we test overlap of both ra1[i] vs ra2[j] and gr.strandflip(ra2[j])
#'
#' @param ra1 \code{GRangesList} with rearrangement set 1
#' @param ra2 \code{GRangesList} with rearrangement set 2
#' @param pad Amount to pad the overlaps by. Larger is more permissive. Default is exact (0)
#' @param arr.ind Default TRUE
#' @param ignore.strand Ignore rearrangement orientation when doing overlaps. Default FALSE
#' @param ... params to be sent to \code{\link{gr.findoverlaps}}
#' @name ra.overlaps
#' @export
ra.overlaps = function(ra1, ra2, pad = 0, arr.ind = TRUE, ignore.strand=FALSE, ...)
{
    bp1 = grl.unlist(ra1) + pad
    bp2 = grl.unlist(ra2) + pad
    ix = gr.findoverlaps(bp1, bp2, ignore.strand = ignore.strand, ...)


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


    .make_matches = function(ix, bp1, bp2)
    {
        if (length(ix) == 0){
            return(NULL)
        }
        tmp.match = cbind(bp1$grl.ix[ix$query.id], bp1$grl.iix[ix$query.id], bp2$grl.ix[ix$subject.id], bp2$grl.iix[ix$subject.id])
        tmp.match.l = lapply(split(1:nrow(tmp.match), paste(tmp.match[,1], tmp.match[,3])), function(x) tmp.match[x, , drop = F])


        matched.l = sapply(tmp.match.l, function(x) all(c('11','22') %in% paste(x[,2], x[,4], sep = '')) | all(c('12','21') %in% paste(x[,2], x[,4], sep = '')))

        return(do.call('rbind', lapply(tmp.match.l[matched.l], function(x) cbind(x[,1], x[,3])[!duplicated(paste(x[,1], x[,3])), , drop = F])))
    }


    tmp = .make_matches(ix, bp1, bp2)

    if (is.null(tmp)){
        if (arr.ind){
            return(matrix())
        }
        else{
            return(Matrix::sparseMatrix(length(ra1), length(ra2), x = 0))
        }
    }

    rownames(tmp) = NULL

    colnames(tmp) = c('ra1.ix', 'ra2.ix')

    if (arr.ind) {
        ro = tmp[order(tmp[,1], tmp[,2]), ]
        if (class(ro)=='integer'){
            ro <- matrix(ro, ncol=2, nrow=1, dimnames=list(c(), c('ra1.ix', 'ra2.ix')))
        }
        return(ro)
    } 
    else {
        ro = Matrix::sparseMatrix(tmp[,1], tmp[,2], x = 1, dims = c(length(ra1), length(ra2)))
        return(ro)
    }
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
#' @title Simple het "caller" meant to be used at validated het SNP sites for tumor / normal pairs
#' @description 
#'
#' hets() outputs a tsv file of ALT ($alt.count.t, $alt.count.n) and REF ($ref.count.t, $ref.count.n) read counts to out.file
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
                }
            }
        }
    }
    return(out)
}



