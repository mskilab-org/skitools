#############################################################################r
## Marcin Imielinski
## The Broad Institute of MIT and Harvard / Cancer program.
## marcin@broadinstitute.org
##
## Weill-Cornell Medical College
## mai9037@med.cornell.edu
##
## New York Genome Center
## mimielinski@nygenome.org
##ski
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
###############################################################################


#
# General utility functions
#

# used for genome plot
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
lsos <- function(..., n=10) {
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

## lapply(list, dim) shortcut

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
#' @return data.datalbe
#' @export
##############
fready = function(..., pattern = "\\W+", sub = "_")
  {
      tab = fread(...)
      setnames(tab, dedup(gsub(pattern, sub, names(tab), perl = TRUE), suffix = '.'))
      return(tab)
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
#' @author Marcin Imielinski, Eran Hodis
#' @export
qq_pval = function(obs, highlight = c(), exp = NULL, lwd = 1, bestfit=T, col = NULL, col.bg='black', pch=18, cex=1, conf.lines=T, max=NULL, qvalues=NULL, label = NULL,
    subsample = NA, ...)
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
        max <- max(obs,exp) + 0.5
    else
        max <- max

    if (is.exp.null)
        {
            tmp.exp = rev(seq(0, 7, 0.01))
            ix = 10^(-tmp.exp)*N
            c95 <-  qbeta(0.975,ix,N-ix+1)
            c05 <-  qbeta(0.025,ix,N-ix+1)

            if (conf.lines){
                ## plot the two confidence lines
                plot(tmp.exp, -log(c95,10), ylim=c(0,max), xlim=c(0,max), type="l", axes=FALSE, xlab="", ylab="")
                par(new=T)
                plot(tmp.exp, -log(c05,10), ylim=c(0,max), xlim=c(0,max), type="l", axes=FALSE, xlab="", ylab="")
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
            names = names(obs[ord])
            setkey(dat, names)
        }

    if (nrow(dat)>1e5) ## rough guide to subsmapling the lower p value part of the plot
        subsample = 5e4/nrow(dat)

    if (is.na(subsample[1]))
        dat[, plot(x, y, xlab = expression(Expected -log[10](italic(P))), ylab = expression(Observed -log[10](italic(P))), xlim = c(0, max), col = colors, ylim = c(0, max), pch=pch, cex=cex, bg=col.bg, ...)]
    else
        {
            subsample = pmin(pmax(0, subsample[1]), 1)
            dat[ifelse(x<=2, ifelse(runif(length(x))<subsample, TRUE, FALSE), TRUE), plot(x, y, xlab = expression(Expected -log[10](italic(P))), ylab = expression(Observed -log[10](italic(P))), xlim = c(0, max), col = colors, ylim = c(0, max), pch=pch, cex=cex, bg=col.bg, ...)]
        }

    if (!is.null(label))
        {
            if (length(label)>0)
                if (is.null(key(dat)))
                    warning('Need to provide names to input vector to draw labels')
                else
                    dat[list(label), text(x, y, labels=label, pos=3)];
        }

    lines(x=c(0, max), y = c(0, max), col = "black", lwd = lwd);

    if (!is.na(subsample))
        dat = dat[sample(nrow(dat), subsample*nrow(dat)), ]

    lambda = lm(y ~ x-1, dat)$coefficients;

    lines(x=c(0, max), y = c(0, lambda*max), col = "red", lty = 2, lwd = lwd);
    legend('bottomright',sprintf('lambda = %.2f', lambda), text.col='red', bty='n')
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




############
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
#' @name rrbind
#' @title rrbind
#'
#' @description
#' like rbind, but takes the intersecting columns of the dfs
#'
#' if union flag is used then will take union of columns (and put NA's for columns of df1 not in df2 and vice versa)
#'
#' @param ... TODO
#' @export
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

#' grdt
#'
#' Converts gr to data frame
#'
#' and a field grl.iix which saves the (local) index that that gr was in its corresponding grl item
#' @param x \code{GRanges} to convert
#' @name grl.unlist
#' @export
grdt = function(x)
    {
        ## new approach just directly instantiating data table
        cmd = 'data.frame(';
        if (is(x, 'GRanges'))
            {
                was.gr = TRUE
                f = c('seqnames', 'start', 'end', 'strand', 'width')
                f2 = c('as.character(seqnames', 'c(start', 'c(end', 'as.character(strand', 'as.numeric(width')
                cmd = paste(cmd, paste(f, '=', f2, '(x))', sep = '', collapse = ','), sep = '')
                value.f = names(values(x))
            }
        else
            {
                was.gr = FALSE
                value.f = names(x)
            }

        if (length(value.f)>0)
            {
                if (was.gr)
                    cmd = paste(cmd, ',', sep = '')
                class.f = sapply(value.f, function(f) eval(parse(text=sprintf("class(x$'%s')", f))))

                .StringSetListAsList = function(x) ### why do I need to do this, bioconductor peeps??
                    {
                        tmp1 = as.character(unlist(x))
                        tmp2 = rep(1:length(x), elementLengths(x))
                        return(split(tmp1, tmp2))
                    }

                ## take care of annoying S4 / DataFrame / data.frame (wish-they-were-non-)issues
                as.statement = ifelse(grepl('Integer', class.f), 'as.integer',
                    ifelse(grepl('Character', class.f), 'as.character',
                           ifelse(grepl('StringSetList', class.f), '.StringSetListAsList',
                                  ifelse(grepl('StringSet$', class.f), 'as.character',
                                         ifelse(grepl('factor$', class.f), 'as.character',
                                                ifelse(grepl('List', class.f), 'as.list',
                                                       ifelse(grepl('factor', class.f), 'as.character',
                                                              ifelse(grepl('List', class.f), 'as.list', 'c'))))))))
                cmd = paste(cmd, paste(value.f, '=', as.statement, "(x$'", value.f, "')", sep = '', collapse = ','), sep = '')
            }

        cmd = paste(cmd, ')', sep = '')

                                        #      browser()
        return(as.data.table(eval(parse(text =cmd))))
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


#################
#' ra_breaks
#'
#' takes in either file or data frame from dranger or snowman or path to BND / SV type vcf file
#' and returns junctions in VCF format.
#'
#' The default output is GRangesList each with a length two GRanges whose strands point AWAY from the break.  If get.loose = TRUE (only relevant for VCF)
#'
#' @name ra_breaks
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb Seqinfo

#' @export
ra_breaks = function(rafile, keep.features = T, seqlengths = hg_seqlengths(), chr.convert = T, snowman = FALSE, swap.header = NULL,  breakpointer = FALSE, seqlevels = NULL, force.bnd = FALSE, 
    get.loose = FALSE ## if TRUE will return a list with fields $junctions and $loose.ends
  )
  {
      if (is.character(rafile))
          {
              if (grepl('(.bedpe$)', rafile))                  
                  {
                      ra.path = rafile
                      
                      cols = c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'name', 'score', 'str1', 'str2')
                      nh = min(which(!grepl('#', readLines(ra.path))))-1
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
#                      library(VariantAnnotation)
                      
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

                      pos1 = as.logical(strand(vgr.pair1)=='+') ## positive strand junctions shift left by one (ie so that they refer to the base preceding the break for these junctions
                      if (any(pos1))
                          {
                              start(vgr.pair1)[pos1] = start(vgr.pair1)[pos1]-1
                              end(vgr.pair1)[pos1] = end(vgr.pair1)[pos1]-1
                          }

                      pos2 = as.logical(strand(vgr.pair2)=='+') ## positive strand junctions shift left by one (ie so that they refer to the base preceding the break for these junctions
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
                 {
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
     
     
     if (keep.features)
         values(out) = rafile[, ]

      if (!get.loose)
          return(out)
      else                              
          return(list(junctions = out, loose.ends = GRanges()))

     return(out)
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
    title.size = 15, footer.size = 20, header.size = round(1.1*data.size))
  {


    # require(hwriter)
    # require(gplots)

    if (!is.data.frame(tab))
      tab = as.data.frame(tab)

    if (is.null(rownames(tab)))
      row.names = F;

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

###############
# writecols
#
# Takes character vector and dumps into k delimited columns
###############
writeCols = function(v, k = 3, sep = "\t", file = "")
  {
    rows = ceiling(length(v)/k)
    out = matrix(ncol = k, nrow = rows, byrow = TRUE)
    out[1:length(v)] = v;
    out[is.na(out)] = "";
    write.table(as.data.frame(out), file = file, sep = sep, quote = F, row.names = F, col.names = F)
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
    return(gsub('[^\\d]', '', as.character(Sys.time()), perl = T))
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
img.html = function(paths, text = names(paths), height = 1024, width = 768)
    {
        if (is.null(text))
            text = ''

        df = data.frame(paths = paths, text = text)
        return(as.vector(rbind(ifelse(!is.na(df$text), paste("<p> <h1>", df$text, "</h1> <p> "), ""),
            paste("<img src = '", df$paths, "' height = '", height, "' width = '", width, "'>", sep = ''))))
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
html_tag = function(tag, text = NULL, collapse = '\n',  ...)
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
        if (length(dim(tab))==1)
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
dfstring = function(df, oneline = TRUE,  sep1 = '; ', sep2 = ', ')
    {
        if (!class(df)[1]=='data.frame')
            df = as.data.frame(df)

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


####################
#' @name cytoscape
#' @title cytoscape
#'
#' @description
#' shortcut to connect to local cytoscape instance running on LOCAL.COMPUTER (unix environment variable) via RPC call
#'
#' graph must be igraph, or adjacency matrixx
#'
#' @param graph igraph or adjacency matrix
#' @param sessionName character name of cytoscape session to open on local computer
#' @param host character name of host on which cytoscape is running (=Sys.getenv('LOCAL.COMPUTER'))
#' @param port integer port on which cytoscape is running (=9000)
#' @param display logical flag whether to display graph locally (=TRUE)
#' @param layout character specifying layout to display as (=degree-circle)
#' @param verbose logical vector whether to use verbose output
#' @param ... additional arguments to new.CytoscapeWindow
#' @export
#' @author Marcin Imielinski
#' @importFrom igraph V E E<- get.edgelist list.edge.attributes get.edge.attribute get.vertex.attribute
cytoscape = function(graph = NULL, sessionName = 'M-ski', host = Sys.getenv('LOCAL.COMPUTER'), port = 9000, display = T, layout = 'degree-circle', verbose = T, ...)
  {
    # require(RCytoscape)
    # require(igraph)

    if (is(graph, 'matrix'))
      graph = igraph::graph.adjacency(graph, weighted = 'weight');

    if (is(graph, 'igraph'))
      {
        if (!is.null(E(graph)$weight))
          E(graph)$weight = 1

        if (!is.null(E(graph)$arrow.shape))
          E(graph)$arrow.shape = 'ARROW'

        graph.nel = igraph2graph(graph)
      }

    cw = RCytoscape::new.CytoscapeWindow(sessionName, host = host, rpcPort = port, graph = graph.nel,  ...)

    if (display)
      {
      RCytoscape::displayGraph(cw)
      RCytoscape::setDefaultBackgroundColor(cw, gplots::col2hex('white'))

        eG = paste(igraph::get.edgelist(graph)[,1], get.edgelist(graph)[,2], sep = '~')
        ceG = RCytoscape::cy2.edge.names(cw@graph)

        if (verbose)
          cat('Setting line styles\n')

        if ('line.style' %in% igraph::list.edge.attributes(graph))
          {
            uls = setdiff(E(graph)$line.style, NA)
            RCytoscape::setEdgeLineStyleRule(cw, 'line.style', uls, uls)
          }

        if (verbose)
          cat('Setting arrow shape\n')

        if (igraph::is.directed(graph))
          if ('arrow.shape' %in% igraph::list.edge.attributes(graph))
            RCytoscape::setEdgeTargetArrowRule(cw, 'arrow.shape', unique(E(graph)$arrow.shape), unique(E(graph)$arrow.shape))

        if (verbose)
          cat('Setting edge color\n')

        if ('col' %in% igraph::list.edge.attributes(graph))
          {
            uc = setdiff(unique(E(graph)$col), NA);
            RCytoscape::setEdgeColorRule(cw, 'col', uc, uc, mode = 'lookup')
          }

        if (verbose)
          cat('Setting edge width\n')

        if ('width' %in% igraph::list.edge.attributes(graph))
          {
            uw = setdiff(igraph::E(graph)$width, NA)
            RCytoscape::setEdgeLineWidthRule(cw, 'width', as.character(uw), uw)
          }

        if (verbose)
          cat('Setting node size\n')

        if ('size' %in% igraph::list.vertex.attributes(graph))
          {
            us = setdiff(unique(igraph::V(graph)$size), NA)
            RCytoscape::setNodeSizeRule(cw, 'size', us, us, mode = 'lookup')
          }

        if (verbose)
          cat('Setting node color\n')

        if ('col' %in% igraph::list.vertex.attributes(graph))
          {
            uc = setdiff(unique(igraph::V(graph)$col), NA)
            RCytoscape::setNodeColorRule(cw, 'col', uc, uc, mode = 'lookup')
          }

        if (verbose)
          cat('Setting node labels\n')

        if ('label' %in% igraph::list.vertex.attributes(graph))
          RCytoscape::setNodeLabelRule(cw, 'label')

        if (verbose)
          cat('Setting node shapes\n')

        if ('shape' %in% igraph::list.vertex.attributes(graph))
          {
            us = setdiff(unique(V(graph)$shape), NA)
            RCytoscape::setNodeShapeRule(cw, 'shape', us, us, default = 'ELLIPSE')
          }

        if (verbose)
          cat('Setting node width\n')

        if ('border.width' %in% igraph::list.vertex.attributes(graph))
          {
            ubw = setdiff(unique(V(graph)$border.width), NA)
            RCytoscape::setNodeBorderWidthRule(cw, 'border.width', ubw, ubw)
          }

        if (all(c('x', 'y') %in% igraph::list.vertex.attributes(graph)))
          {
            good.ix = !is.na(V(graph)$x) & !is.na(V(graph)$y)
            if (any(good.ix))
              RCytoscape::setNodePosition(cw, V(graph)$name[good.ix], V(graph)$x[good.ix], V(graph)$y[good.ix])
          }
        else
          RCytoscape::layoutNetwork(cw, layout)


        RCytoscape::redraw(cw)
      }

    return(cw)
  }


####################
#' @name igraph2graph
#' @title igraph2graph
#'
#' @description
#' Converts igraph object into janky graphNEL object (for visualization in cytoscape)
#' and populates all edge features both via the edgeL and as NodeAttributes for visualization
#' in cytoscape
#'
#' @importFrom graph edgeData<- nodeData<-
#' @param g igraph object
#' @author Marcin Imielinski
#' @export
#' @return graph object
igraph2graph = function(g)
  {
    # require(igraph)
    # require(graph)
    # require(RCytoscape)

    if (class(V) != 'function' | class(E) != 'function')
      stop('Namespace conflict - either V() or E() no longer mapping to igraph functions')

    if (!is.null(V(g)$name))
      node.labels = V(g)$name
    else
      node.labels = as.character(V(g));

    edge.df = structure(as.data.frame(get.edgelist(g), stringsAsFactors = F), names = c('vertices', 'edges'))
    if (length(igraph::list.edge.attributes(g))>0)
      {
        tmp = do.call('cbind', lapply(igraph::list.edge.attributes(g), function(x) get.edge.attribute(g, x)))
        colnames(tmp) = igraph::list.edge.attributes(g)
        edge.df = cbind(edge.df, as.data.frame(tmp, stringsAsFactors = F))
      }
    if (!is.null(edge.df$weights))
      edge.df$weights = as.numeric(edge.df$weights)
    edge.df[is.na(edge.df)] = "NA"
    edge.df[,1] = as.character(edge.df[,1])
    edge.df[,2] = as.character(edge.df[,2])

    vertex.df = data.frame(vertices = node.labels, stringsAsFactors = F)
    if (length(igraph::list.vertex.attributes(g))>0)
      {
        tmp = do.call('cbind', lapply(igraph::list.vertex.attributes(g), function(x) get.vertex.attribute(g, x)))
        colnames(tmp) = igraph::list.vertex.attributes(g)
        vertex.df = cbind(vertex.df, as.data.frame(tmp, stringsAsFactors = F))
      }

    ## have to reciprocate edges in undirected otherwise graphNEL will barf
    if (!igraph::is.directed(g))
      {
        edge.df.rev = edge.df;
        tmp.col = edge.df[,2]
        edge.df.rev[,2] = edge.df.rev[,1]
        edge.df.rev[,1] = tmp.col;
        edge.df = rbind(edge.df, edge.df.rev)
      }

    edgeL = lapply(split(edge.df, edge.df$vertices), as.list)[node.labels]
    names(edgeL) = node.labels;

    ## retarded GraphNEL object format necessitates these gymnastics
    null.vert = sapply(edgeL, is.null)
    blank.edge.item = list(structure(rep(list(c()), length(igraph::list.edge.attributes(g))+1), names = c('edges', igraph::list.edge.attributes(g))))
    edgeL[null.vert] = blank.edge.item
    edgemode = c('undirected', 'directed')[1 + igraph::is.directed(g)]

    out.g = new('graphNEL', node.labels, edgeL, edgemode = edgemode)

    ## populate edge and node attribute for RCytoscape to access
    if (ncol(edge.df)>2)
      {
        attr.cols = names(edge.df)[3:ncol(edge.df)]
        for (attr in attr.cols)
          {
            if (is.numeric(edge.df[, attr]))
              out.g = RCytoscape::initEdgeAttribute(out.g, attr, 'numeric', NA)
            else if (is.integer(edge.df[, attr]))
              out.g = RCytoscape::initEdgeAttribute(out.g, attr, 'integer', NA)
            else
              {
                cast.numeric = suppressWarnings(as.numeric(edge.df[, attr]))
                cast.character = edge.df[, attr]

                if (any(is.na(cast.numeric) & cast.character != "NA"))
                  {
                    edge.df[, attr] = as.character(edge.df[, attr])
                    out.g = RCytoscape::initEdgeAttribute(out.g, attr, 'char', '')
                  }
                else
                  {
                    cast.numeric[is.na(cast.numeric)] = 0
                    edge.df[, attr] = cast.numeric
                    out.g = RCytoscape::initEdgeAttribute(out.g, attr, 'numeric', '')
                  }
              }
            edgeData(out.g, edge.df[,1], edge.df[,2], attr) = edge.df[,attr]
          }
      }

    if (ncol(vertex.df)>1)
      {
        attr.cols = names(vertex.df)[2:ncol(vertex.df)]
        for (attr in attr.cols)
          {
            if (is.numeric(vertex.df[, attr]))
              out.g = RCytoscape::initNodeAttribute(out.g, attr, 'numeric', NA)
            else if (is.integer(vertex.df[, attr]))
              out.g = RCytoscape::initNodeAttribute(out.g, attr, 'integer', NA)
            else
              {
                cast.numeric = suppressWarnings(as.numeric(vertex.df[, attr]))
                cast.character = suppressWarnings(as.character(vertex.df[, attr]))
                if (any(setdiff(is.na(cast.numeric) & cast.character != "NA", NA)))
                  {
                    vertex.df[, attr] = as.character(vertex.df[, attr])
                    out.g = RCytoscape::initNodeAttribute(out.g, attr, 'char', '')
                  }
                else
                  {
                    cast.numeric[is.na(cast.numeric)] = 0
                    vertex.df[, attr] = cast.numeric
                    out.g = RCytoscape::initNodeAttribute(out.g, attr, 'char', '')
                  }
              }
            nodeData(out.g, vertex.df[,1], attr) = vertex.df[,attr]
          }
      }

    return(out.g)
  }

#' @name cyto2igraph
#' @title cyto2igraph
#'
#' @description
#' Pulls graph in CytoscapeWindow "cw" and returns as igraph object
#'
#' @param cw CytoscapeWindow object to grab from (see Rcytoscape)
#' @author Marcin Imielinski
#' @return igraph object
#' @export
cyto2igraph = function(cw)
{
  node.attr = RCytoscape::getAllNodeAttributes(cw)
  edge.attr = RCytoscape::getAllEdgeAttributes(cw)
  directed = graph::edgemode(cw@graph) == 'directed'

  edge.df = cbind(from = edge.attr$source, to = edge.attr$target,
    edge.attr[, setdiff(names(edge.attr), c('source', 'target'))])

  if ('name' %in% node.attr)
    node.df = node.attr[, c('name', setdiff(names(node.attr), 'name'))]
  else
    node.df = cbind(name = rownames(node.attr), node.attr)

  return(igraph::graph.data.frame(edge.df, directed = directed, vertices = node.df))
}

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
vaggregate = function(...)
  {
    out = aggregate(...);
    return(structure(out[,ncol(out)], names = do.call(paste, lapply(names(out)[1:(ncol(out)-1)], function(x) out[,x]))))
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
#' dbeta: shape1, sheape 2
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
        out = aggregate(prob ~ x, data = out, FUN = sum)
        out$id = 1;
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
            out = matrix(as.integer(rid), as.integer(cid), x = tmp[, 'val'], dimnames = list(levels(rid), levels(cid)))

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
mmatch = function(A, B, dir = 1)
  {
    SEP = ' ';
    Atxt = apply(A, dir, function(x) paste(x, collapse = SEP))
    Btxt = apply(B, dir, function(x) paste(x, collapse = SEP))

    return(match(Atxt, Btxt))
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
#' @name gr.flip
#' @title gr.flip
#'
#' @description
#' alias to gr.flipstrand for backward compatibility
#'
#' @param ... args to gr.flipstrand in gUtils
#' @return flipped GRAnges
#' @author Marcin Imielinski
#' @export
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
#' @param dotsize integer dot size to plot with (= NULL)
#' @param binwidth numeric binwidth of histogram (= NULL)
#' @param title character title of plot (='')
#' @param ylim y limits of plot (= NULL)
#' @param text.size text size of legend (= NULL)
#' @author Marcin Imielinski
#' @export
dplot = function(y, group, ylab = '', xlab = '', log = F, dotsize = NULL, binwidth = NULL, title = NULL, ylim = NULL, text.size = NULL)
  {

    df = data.frame(y = y, group = as.character(group), stringsAsFactors = F)

    if (is.null(binwidth))
        binwidth = as.numeric((quantile(y, c(0.95)) - quantile(y, c(0.05)))/500)
#        binwidth = as.numeric(quantile(diff(sort(y)), 0.3)*10)

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
        names(out) = gsub(pattern, rep, dir(x, pattern))
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
#' @name pad
#' @title pad
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
pad = function(x, k, clip = T)
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

##############################
#' @name affine.map
#' @title affine.map
#'
#' @description
#' affinely maps 1D points in vector x from interval xlim to interval ylim,
#' ie takes points that lie in
#' interval xlim and mapping onto interval ylim using linear / affine map defined by:
#' (x0,y0) = c(xlim(1), ylim(1)),
#' (x1,y1) = c(xlim(2), ylim(2))
#' (using two point formula for line)
#' useful for plotting.
#'
#' if cap.max or cap.min == T then values outside of the range will be capped at min or max
#'
#' @param x input vector to affinely transform
#' @param ylim length 2 vector of y bounds to project data into (=c(0,1))
#' @param xlim length 2 vector of x bounds to project data from (=c(min(x), max(x)))
#' @param cap logical flag whether to cap data below xmin and above xmax of x vector ie return as ymin, ymax respectively (=FALSE)
#' @param cap.min logical flag whether to cap data below xmin  of x vector ie return as ymin (=cap)
#' @param cap.max logical flag whether to cap  data above xmax of x vector ie return as ymax (=cap)
#' @param clip logical flag whether to clip data below xmin and above xmax of x vector ie return as NA (=TRUE)
#' @param clip.min logical flag whether to clip data below xmin  of x vector ie return as NA (=clip)
#' @param clip.max logical flag whether to clip data above xmax  of x vector ie return as NA (=clip)
#' @return x values transformed by affine map defined by xlim and ylim
#' @export
#' @author Marcin Imielinski
 ##############################
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
    png(filename, height = height, width = width, pointsize = 12*cex.pointsize, ...)

    if (!is.null(dim))
        {
            if (length(dim)==1)
                dim = rep(dim, 2)
            dim = dim[1:2]
            layout(matrix(1:prod(dim), nrow = dim[1], ncol = dim, byrow = TRUE))
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
plop = function(fn, prefix = NULL)
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
                dat = dat[!ix, ]

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
vplot = function(y, group, facet1 = NULL, facet2 = NULL, transpose = FALSE, mapping = NULL, stat = "ydensity", position = "dodge", trim = TRUE, scale = "area", log = FALSE, count = TRUE, xlab = NULL, ylim = NULL, ylab = NULL, minsup = NA, scatter = FALSE, cex.scatter = 1, col.scatter = NULL, alpha = 0.3, title = NULL, legend.ncol = NULL, legend.nrow = NULL, vfilter = TRUE, vplot = TRUE, dot = FALSE, stackratio = 1, binwidth = 0.1)
    {
        # require(ggplot2)
        if (!is.factor(group))
            group = as.factor(group)
        dat = data.table(y = suppressWarnings(as.numeric(y)), group)

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
                        g = g + geom_dotplot(data = dat, mapping = aes(x = group, y = y, fill = group), binaxis = 'y', position = 'identity', col = NA, alpha = alpha, method = 'dotdensity', dotsize = cex.scatter, stackratio = stackratio, binwidth = binwidth, stackdir = 'center')
                    }
                else
                    {
                        if (is.null(col.scatter))
                            g = g + geom_jitter(data = dat, mapping = aes(fill = group), shape = 21, size = cex.scatter, alpha = alpha, position = position_jitter(height = 0))
                        else
                            g = g + geom_jitter(data = dat, fill = alpha(col.scatter, shape = 21, alpha), position = position_jitter(height = 0))
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

        print(g)

        'voila'
    }



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


############################
# standardize_segs
#
# (data frame seg function)
#
# Takes and returns segs data frame standardized to a single format (ie $chr, $pos1, $pos2)
#
# if chr = T will ensure "chr" prefix is added to chromossome(if does not exist)
#############################
standardize_segs = function(seg, chr = F)
{
  if (inherits(seg, 'IRangesList'))
    seg = irl2gr(seg);

  if (is(seg, 'matrix'))
    seg = as.data.frame(seg, stringsAsFactors = F)

  if (inherits(seg, 'RangedData') | inherits(seg, 'GRanges') | inherits(seg, 'IRanges'))
  {
    val = as.data.frame(values(seg));
    values(seg) = NULL;
    seg = as.data.frame(seg, row.names = NULL);  ## returns compressed iranges list
    seg$seqnames = as.character(seg$seqnames)
  }
  else
    val = NULL;

  field.aliases = list(
    ID = c('id', 'patient', 'Sample'),
    chr = c('chrom', 'Chromosome', "contig", "seqnames", "seqname", "space", 'chr'),
    pos1 = c('loc.start', 'begin', 'Start', 'start', 'Start.bp', 'Start_position', 'pos', 'pos1', 'left', 's1'),
    pos2 =  c('loc.end', 'End', 'end', "stop", 'End.bp', 'End_position', 'pos2', 'right', 'e1'),
    strand = c('str', 'strand', 'Strand', 'Str')
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
    muts = cov*0;  # compute muts from maf file
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


###########
# pmGSEA "poor man's GSEA " ***
#
# Given a gene.set (character vector) or gene.sets (list of character vectors)
# and given a named vector of significance values or table of significant genes (sig.table)
# (if table then significance column is $p or first column) identifies gene sets that have significant
# negative deviation of a "signed K-S" statistic vs uniform distribution  (ie have p values significantly
# clustering towards zero) ie are significantly enriched in genes showing positive selection.
#
# if positive.selection = F, will identify sets with significantly positive deviation of a "signed K-S" statistic (ie have p values significantly clustering towards 1)
# these are sets showig significant negative selection.
#
# All p-values are computed against a distribution of signed K-S statistic obtained through permutation using random gene sets of the same size chosen from sig.table
#
# Will adaptively perform permutations between minperms and maxperms using following rule of thumb: if there are <PERM.THRESH permutations with
# greater than (lower.tail = F) or less than (lower.tail = T) score than observed score, then will compute additional perms
#
# *** actually not much poorer than the original GSEA, basically a reimplementation of Mootha et al Nat Gen 2002
#
###########
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
    sort.locus = NULL,
    gsub.paths = list(
        c('/seq/picard_aggregation/', '/Volumes/seq_picard_aggregation/'),
        c('/xchip/singtex/', '/Volumes/xchip_singtex/'),
        c('/cga/meyerson/', '/Volumes/cga_meyerson/'),
        c('/seq/tier3b/', '/Volumes/seq_tier3b/'),
        c('/xchip/cle/', '/Volumes/xchip_cle/'),
        c('/xchip/cga/', '/Volumes/xchip_cga/'),
        c('/xchip/beroukhimlab/', '/Volumes/xchip_beroukhimlab/'),
        c('/cgaext/tcga/', '/Volumes/cgaext_tcga/'),
        c('~/', '~/home/'),
        c(Sys.getenv('HOME'), '~/home')),
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
        stop('IGV host field is empty, either specify via host argument to this function or set IGV_HOST environment variable')

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

                        new.paths = normalizePath(paths);

                        if (mac)
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

#' @name dcast.data.table2
#' @title dcast.data.table but allows vector arguments for value.var,
#' @description
#' if value.var is a vector then will combine the right hand side column names with each element of value.var
#' in a merged cast table
#'
#' @export
#' @author Marcin Imielinski
dcast.data.table2 = function(data, formula, ..., value.var = NULL, sep = '_')
    {
        terms = sapply(unlist(as.list(attr(terms(formula), "variables"))[-1]), as.character)
        if (is.null(value.var))
            value.var = setdiff(colnames(data), terms)
        dt = lapply(value.var, function(x)
            {
                d = data.table::dcast.data.table(data, formula,  ..., value.var = x)
                new.cols = setdiff(colnames(d), key(d))
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

#' Read BAM file into GRanges or data.table
#'
#' Wrapper around Rsamtools bam scanning functions,
#' by default, returns GRangesList of read pairs for which <at least one> read lies in the supplied interval
#' @param bam Input bam file. Advisable to make "bam" a BamFile instance instead of a plain string, so that the index does not have to be reloaded.
#' @param bami Input bam index file.
#' @param gr GRanges of intervals to retrieve
#' @param intervals GRanges of intervals to retrieve
#' @param stripstrand Flag to ignore strand information on the query intervals. Default TRUE
#' @param what What fields to pull down from BAM. Default \code{scanBamWhat()}
#' @param unpack.flag Add features corresponding to read flags. Default FALSE
#' @param verbose Increase verbosity
#' @param tag Additional tags to pull down from the BAM (e.g. 'R2')
#' @param isPaired See documentation for \code{scanBamFlag}. Default NA
#' @param isProperPair See documentation for \code{scanBamFlag}. Default NA
#' @param isUnmappedQuery See documentation for \code{scanBamFlag}. Default NA
#' @param hasUnmappedMate See documentation for \code{scanBamFlag}. Default NA
#' @param isNotPassingQualityControls See documentation for \code{scanBamFlag}. Default NA
#' @param isDuplicate See documentation for \code{scanBamFlag}. Default FALSE
#' @param isValidVendorRead See documentation for \code{scanBamFlag}. Default TRUE
#' @param as.grl Return reads as GRangesList. Controls whether \code{get.pairs.grl} does split. Default TRUE
#' @param as.data.table Return reads in the form of a data.table rather than GRanges/GRangesList
#' @param ignore.indels messes with cigar to read BAM with indels removed. Useful for breakpoint mapping on contigs
#' @param size.limit Default 1e6
#' @param ... passed to \code{scanBamFlag}
#' @return Reads in one of GRanges, GRangesList or data.table
#' @export
read.bam = function(bam, intervals = NULL,## GRanges of intervals to retrieve
                    gr = intervals,
                    all = FALSE,
                    bami = NULL,
                    pairs.grl = TRUE, # if TRUE will return GRangesList of read pairs for whom at least one read falls in the supplied interval
                                        #  paired = F, # if TRUE, will used read bam gapped alignment pairs warning: will throw out pairs outside of supplied window
                                        #  gappedAlignment = T, # if false just read alignments using scanbam
                    stripstrand = TRUE,
                    what = scanBamWhat(),
                    unpack.flag = FALSE, # will add features corresponding to read flags
                    verbose = FALSE,
                    tag = NULL,
                    isPaired = NA, ## if these features are NA, then reads satisfying both T and F will be returned
                    isProperPair = NA,
                    isUnmappedQuery = NA,
                    hasUnmappedMate = NA,
                    isNotPassingQualityControls = NA,
                    isDuplicate = F,
                    isValidVendorRead = TRUE,
                    as.grl=TRUE, ## return pairs as grl, rather than GRanges .. controls whether get.pairs.grl does split (t/c rename to pairs.grl.split)
                    as.data.table=FALSE, ## returns reads in the form of a data table rather than GRanges/GRangesList
                    ignore.indels=FALSE, ## messes with cigar to read BAM with indels removed. Useful for breakpoint mapping on contigs
                    size.limit = 1e6,
                    ... # passed to scanBamFlag (
                    )
{
    if (!inherits(bam, 'BamFile'))
    {
        if (is.null(bami))
        {
            if (file.exists(bai <- gsub('.bam$', '.bai', bam)))
                bam = BamFile(bam, bai)
            else if (file.exists(bai <- paste(bam, '.bai', sep = '')))
                bam = BamFile(bam, bai)
            else
                bam = BamFile(bam)
        }
        else
            bam = BamFile(bam, index = bami)
    }


                                        # if intervals unspecified will try to pull down entire bam file (CAREFUL)

    if (length(intervals)==0)
        intervals = NULL

    if (is.null(intervals))
        intervals = gr

    if (is.null(intervals))
    {
        if (all)
            intervals = seqinfo2gr(seqinfo(bam))
        else
            stop('Must provide non empty interval list')
    }

    if (class(intervals) == 'data.frame')
        intervals = seg2gr(intervals);

    if (inherits(intervals, 'GRangesList'))
        intervals = unlist(intervals);

    if (stripstrand)
        strand(intervals) = '*'

    intervals = reduce(intervals);

    now = Sys.time();

    if (pairs.grl)
        paired = F

    flag = scanBamFlag(isPaired = isPaired, isProperPair = isProperPair, isUnmappedQuery = isUnmappedQuery,
                       hasUnmappedMate = hasUnmappedMate, isNotPassingQualityControls = isNotPassingQualityControls,
                       isDuplicate = isDuplicate, ...)

    tag = unique(c('MD', 'MQ', tag))
    param = ScanBamParam(which = gr.fix(intervals, bam, drop = T), what = what, flag = flag, tag = tag)

    if (verbose)
        cat('Reading bam file\n')
    if (class(bam) == 'BamFile')
        out <- scanBam(bam, param=param)
    else
        out <- scanBam(bam, index=bami, param=param)
    if (verbose) {
        print(Sys.time() - now)
        print('BAM read. Making into data.frame')
    }

    out <- out[sapply(out, function(x) length(x$qname)>0)]

    if (length(out)>0)
    {
        if (verbose) {
            print(Sys.time() - now)
            print('combining lists')
        }
        out <- as.data.frame(rbindlist(lapply(out, function(x)
        {
            x <- c(x[-match('tag', names(x))], x$tag)

            x <- x[sapply(x, length)>0]
            conv <- which(!(sapply(x, class) %in% c('integer', 'numeric', 'character')))
            x[conv] <- lapply(x[conv], as.character)

            for (t in tag)
                if (!(t %in% names(x)))
                    x[[t]] = rep(NA, length(x$qname))

            if (!('R2' %in% names(x)) && 'R2' %in% tag)
                x$R2 <- rep(NA, length(x$qname))
            if (!('Q2' %in% names(x)) && 'Q2' %in% tag)
                x$Q2 <- rep(NA, length(x$qname))
            x
        })))

        ## faster CIGAR string parsing with vectorization and data tables
        if (verbose) {
            print(Sys.time() - now)
            print('filling pos2 from cigar')
        }
        if (ignore.indels) {
            cigar <- gsub('[0-9]+D', '', gsub('([0-9]+)I', '\\1M', out$cigar))  ## Remove deletions, turn insertions to matches
            cig <- splitCigar(cigar)
            torun=sapply(cig, function(y) any(duplicated((y[[1]][y[[1]]==M]))))
            M <- charToRaw('M')
            new.cigar <- sapply(cig[torun], function(y) {
                lets <- y[[1]][!duplicated(y[[1]])]
                vals <- y[[2]][!duplicated(y[[1]])]
                vals[lets==M] <- sum(y[[2]][y[[1]]==M])
                lets <- strsplit(rawToChar(lets), '')[[1]]
                paste(as.vector(t(matrix(c(vals, lets), nrow=length(vals), ncol=length(lets)))), collapse='')
            })
            out$cigar[torun] <- new.cigar
        }
        cigs <- countCigar(out$cigar)
        out$pos2 <- out$pos + cigs[, "M"]

        if (verbose) {
            print(Sys.time() - now)
            print('fixing seqdata')
        }
        out$qwidth = nchar(out$seq)
        unm = is.na(out$pos)
        if (any(unm))
        {
            out$pos[unm] = 1
            out$pos2[unm] = 0
            out$strand[unm] = '*'
        }
        gr.fields = c('rname', 'strand', 'pos', 'pos2');
        vals = out[, setdiff(names(out), gr.fields)]

        if (!as.data.table) {
            out <- GRanges(out$rname, IRanges(out$pos, pmax(0, out$pos2-1)), strand = out$strand, seqlengths = seqlengths(intervals))
            values(out) <- vals;
        } else {
            out <- data.table(seqnames=out$rname, start=out$pos, end= pmax(out$pos2-1, 0), strand=out$strand)
            val <- data.table(vals)
            out <- cbind(out, val)
        }
                                        #out$uname = paste(out$qname, ifelse(bamflag(out$flag)[, 'isFirstMateRead'], '_r1', '_r2'), sep = '')
    }
    else {
        if (!as.data.table)
            return(GRanges(seqlengths = seqlengths(intervals)))
        else
            return(data.table())
    }

    if (verbose)
    {
        if (as.data.table)
            cat(sprintf('Extracted %s reads\n', nrow(out)))
        else
            cat(sprintf('Extracted %s reads\n', length(out)))
        print(Sys.time() - now)
    }

    if (pairs.grl)
    {
        if (verbose)
            cat('Pairing reads\n')
        out <- get.pairs.grl(out, as.grl=as.grl)
        if (verbose)
        {
            cat('done\n')
            print(Sys.time() - now)
        }
        if (as.grl && !as.data.table) {
            names(out) = NULL;
            values(out)$col = 'gray';
            values(out)$border = 'gray';
        }
    }

    return(out)
}


#' @name bam.cov.gr
#' @title Get coverage as GRanges from BAM on custom set of GRanges
#' @description
#'
#' gets coverage from bam in supplied ranges using "countBam", returning gr with coverage counts in
#' each of the provided ranges (different from bam.cov above) specified as $file, $records, and $nucleotides
#' columns in the values field
#' basically a wrapper for countBam with some standard settings for ScanBamParams
#'
#' @param bam Input bam file. Advisable to make "bam" a BamFile instance instead of a plain string, so that the index does not have to be reloaded.
#' @param bami Input bam index file.
#' @param gr GRanges of intervals to retrieve
#' @param verbose Increase verbosity
#' @param isPaired See documentation for \code{scanBamFlag}. Default NA
#' @param isProperPair See documentation for \code{scanBamFlag}. Default NA
#' @param isUnmappedQuery See documentation for \code{scanBamFlag}. Default NA
#' @param hasUnmappedMate See documentation for \code{scanBamFlag}. Default NA
#' @param isNotPassingQualityControls See documentation for \code{scanBamFlag}. Default NA
#' @param isDuplicate See documentation for \code{scanBamFlag}. Default FALSE
#' @param isValidVendorRead See documentation for \code{scanBamFlag}. Default TRUE
#' @param mc.cores Number of cores in \code{mclapply} call.
#' @param chunksize How many intervals to process per core. Default 10.
#' @param ... passed to \code{scanBamFlag}
#' @return GRanges parallel to input GRanges, but with metadata filled in.
#' @export
bam.cov.gr = function(bam, gr, bami = NULL, count.all = FALSE, isPaired = T, isProperPair = T, isUnmappedQuery = F, hasUnmappedMate = F, isNotPassingQualityControls = F, isDuplicate = F, isValidVendorRead = T, mc.cores = 1, chunksize = 10, verbose = F, ...)
{
    if (is.character(bam))
        if (!is.null(bami))
            bam = BamFile(bam, bami)
        else
        {
            if (file.exists(paste(bam, 'bai', sep = '.')))
                bam = BamFile(bam, paste(bam, 'bai', sep = '.'))
            else if (file.exists(gsub('.bam$', '.bai', bam)))
                bam = BamFile(bam, paste(bam, 'bai', sep = '.'))
            else
                stop('BAM index not found, please find index and specify bam file argument as valid BamFile object')
        }

    keep = which(seqnames(gr) %in% seqlevels(bam))

    if (length(keep)>0)
    {
        ix = c(keep[c(seq(1, length(keep), chunksize))], keep[length(keep)]+1);  ## prevent bam error from improper chromosomes
        chunk.id = unlist(lapply(1:(length(ix)-1), function(x) rep(x, ix[x+1]-ix[x])))

        gr.chunk = split(gr[keep], chunk.id[keep]);
        if (count.all)
            flag = scanBamFlag()
        else
            flag = scanBamFlag(isPaired = isPaired, isProperPair = isProperPair, isUnmappedQuery = isUnmappedQuery,
                               hasUnmappedMate = hasUnmappedMate, isNotPassingQualityControls = isNotPassingQualityControls,
                               isDuplicate = isDuplicate, ...)
        out = rbindlist(mclapply(1:length(gr.chunk),
                                 function(x) {
                                     if (verbose)
                                         cat(sprintf('Processing ranges %s to %s of %s, extracting %s bases\n', ix[x], ix[x+1]-1, length(keep), sum(width(gr.chunk[[x]]))))
                                     as.data.table(countBam(bam, param = ScanBamParam(which = gr.chunk[[x]], flag = flag)))
                                 }, mc.cores = mc.cores));

        gr.tag = paste(as.character(seqnames(gr)), start(gr), end(gr));
        out.tag = paste(out$space, out$start, out$end);
        ix = match(gr.tag, out.tag);
        values(gr) = cbind(as.data.frame(values(gr)), out[ix, c('file', 'records', 'nucleotides'), with = FALSE])
    }
    else
        values(gr) = cbind(as.data.frame(values(gr)), data.frame(file = rep(gsub('.*\\/([^\\/]+)$', '\\1', path(bam)), length(gr)), records = NA, nucleotides = NA))

    return(gr)
}

#' @name bam.cov.tile
#' @title Get coverage as GRanges from BAM on genome tiles across seqlengths of genome
#' @description
#' Quick way to get tiled coverage via piping to samtools (~10 CPU-hours for 100bp tiles, 5e8 read pairs)
#'
#' Gets coverage for window size "window", pulling "chunksize" records at a time and incrementing bin
#' corresponding to midpoint or overlaps of corresponding (proper pair) fragment (uses TLEN and POS for positive strand reads that are part of a proper pair)
#'
#' @param bam.file character scalar input bam file
#' @param window integer scalar window size (in bp)
#' @param chunksize integer scalar, size of window
#' @param min.mapq integer scalar, minimim map quality reads to consider for counts
#' @param verbose dummy
#' @param max.tlen max paired-read insert size to consider
#' @param st.flag samtools flag to filter reads on [Default: -f 0x02 -F 0x10]
#' @param fragments dummy
#' @param region dummy
#' @param do.gc dummy
#' @param midpoint if TRUE will only use the fragment midpoint, if FALSE will count all bins that overlap the fragment
#' @return GRanges of "window" bp tiles across seqlengths of bam.file with meta data field $counts specifying fragment counts centered
#' in the given bin.
#' @export
bam.cov.tile = function(bam.file, window = 1e2, chunksize = 1e5, min.mapq = 30, verbose = TRUE,
                        max.tlen = 1e4, ## max insert size to consider
                        st.flag = "-f 0x02 -F 0x10",
                        fragments = TRUE,
                        region = NULL,
                        do.gc = FALSE,
                        midpoint = TRUE ## if TRUE will only use the fragment midpoint, if FALSE will count all bins that overlap the fragment
                        )
{
    cmd = 'samtools view %s %s -q %s | cut -f "3,4,9"' ## cmd line to grab the rname, pos, and tlen columns

    sl = seqlengths(BamFile(bam.file))

    counts = lapply(sl, function(x) rep(0, ceiling(x/window)))
    numwin = sum(sapply(sl, function(x) ceiling(x/window)))

    if (!is.null(region))
    {
        cat(sprintf('Limiting to region %s\n', region))
        cmd = 'samtools view %s %s -q %s %s | cut -f "3,4,9"' ## cmd line to grab the rname, pos, and tlen columns
        if (!file.exists(paste(bam.file, '.bam', sep = '')))
            if (file.exists(bai.file <- gsub('.bam$', '.bai', bam.file)))
            {
                .TMP.DIR = '~/temp/.samtools'
                system(paste('mkdir -p', TMP.DIR))
                tmp.fn = paste(normalizePath(TMP.DIR), '/tmp', runif(1), sep = '')
                system(sprintf('ln -s %s %s.bam', bam.file, tmp.fn))
                system(sprintf('ln -s %s %s.bam.bai', bai.file, tmp.fn))
            }

        cat('Calling', sprintf(cmd, st.flag, paste(tmp.fn, 'bam', sep = '.'), min.mapq, region), '\n')
        p = pipe(sprintf(cmd, st.flag, paste(tmp.fn, 'bam', sep = '.'), min.mapq, region), open = 'r')
    }
    else
    {
        cat('Calling', sprintf(cmd, st.flag, bam.file, min.mapq), '\n')
        p = pipe(sprintf(cmd, st.flag, bam.file, min.mapq), open = 'r')
    }

    i = 0
    sl.dt = data.table(chr = names(sl), len = sl)
    counts = sl.dt[, list(start = seq(1, len, window)), by = chr]
    counts = counts[, bin := 1:length(start), by = chr]
    counts[, end := pmin(start + window-1, sl[chr])]
    counts[, count := 0]
    counts[, rowid := 1:length(count)]
    setkeyv(counts, c("chr", "bin")) ## now we can quickly populate the right entries
    totreads = 0

    st = Sys.time()
    if (verbose)
        cat('Starting fragment count on', bam.file, 'with bin size', window, 'and min mapQ', min.mapq, 'and insert size limit', max.tlen, 'with midpoint set to', midpoint, '\n')

    while (length(chunk <- readLines(p, n = chunksize))>0)
    {
        i = i+1

        if (fragments)
        {
            chunk = fread(paste(chunk, collapse = "\n"), header = F)[abs(V3)<=max.tlen, ]
            if (midpoint) ## only take midpionts
                chunk[, bin := 1 + floor((V2 + V3/2)/window)] ## use midpoint of template to index the correct bin
            else ## enumerate all bins containing fragment i.e. where fragments overlap multiple bins  (slightly slower)
            {
                if (verbose)
                    cat('!!!! Counting all overlapping bins !!!\n')
                chunk[, ":="(bin1 = 1 + floor((V2)/window), bin2 = 1 + floor((V2+V3)/window))]
                chunk = chunk[, list(V1, bin = bin1:bin2), by = list(ix = 1:length(V1))]
            }
        }
        else ## just count reads
        {
            cat('counting reads\n')
            chunk = fread(paste(chunk, collapse = "\n"), header = F)
            chunk[, bin := 1 + floor((V2)/window)]
        }

        tabs = chunk[, list(newcount = length(V1)), by = list(chr = as.character(V1), bin)] ## tabulate reads to bins data.table style
        counts[tabs, count := count + newcount] ## populate latest bins in master data.table

        ## should be no memory issues here since we preallocate the data table .. but they still appear
        if (do.gc)
        {
            print('GC!!')
            print(gc())
        }
                                        #          print(tables())

        ## report timing
        if (verbose)
        {
            cat('bam.cov.tile.st ', bam.file, 'chunk', i, 'num fragments processed', i*chunksize, '\n')
            timeelapsed = as.numeric(difftime(Sys.time(), st, units = 'hours'))
            meancov = i * chunksize / counts[tabs[nrow(tabs),], ]$rowid  ## estimate comes from total reads and "latest" bin filled
            totreads = meancov * numwin
            tottime = totreads*timeelapsed/(i*chunksize)
            rate = i*chunksize / timeelapsed / 3600
            cat('mean cov:', round(meancov,1), 'per bin, estimated tot fragments:', round(totreads/1e6,2), 'million fragments, processing', rate,
                'fragments/second\ntime elapsed:', round(timeelapsed,2), 'hours, estimated time remaining:', round(tottime - timeelapsed,2), 'hours', ', estimated total time', round(tottime,2), 'hours\n')
        }
    }

    gr = GRanges(counts$chr, IRanges(counts$start, counts$end), count = counts$count, seqinfo = Seqinfo(names(sl), sl))
    if (verbose)
        cat("Finished computing coverage, and making GRanges\n")
    close(p)

    if (!is.null(region))
        system(sprintf('rm %s.bam %s.bam.bai', tmp.fn, tmp.fn))

    return(gr)
}

#' Compute rpkm counts from counts
#'
#' takes countbam (or bam.cov.gr) output "counts" and computes rpkm by aggregating across "by" variable
#' @param counts GRanges, data.table or data.frame with records, width fields
#' @param by Field to group counts by
#' @note The denominator (ie total reads) is just the sum of counts$records
#' @export
counts2rpkm = function(counts, by)
{
    out = aggregate(1:nrow(counts), by = list(by), FUN = function(x) sum(counts$records[x])/ sum(counts$width[x]/1000));
    out[,2] = out[,2]/sum(counts$records)*1e6;
    names(out) = c('by', 'rpkm');
    return(out);
}

#' Create GRanges of read mates from reads
#'
#' @return \code{GRanges} corresponding to mates of reads
#' @name get.mate.gr
#' @export
get.mate.gr = function(reads)
{

    if (inherits(reads, 'GRanges')) {
        mpos = values(reads)$mpos
        mrnm = as.vector(values(reads)$mrnm)
        mapq = values(reads)$MQ
        bad.chr = !(mrnm %in% seqlevels(reads)); ## these are reads mapping to chromosomes that are not in the current "genome"
        mrnm[bad.chr] = as.character(seqnames(reads)[bad.chr]) # we set mates with "bad" chromosomes to have 0 width and same seqnames (ie as if unmapped)
    } else if (inherits(reads, 'data.table')) {
        mpos <- reads$mpos
        mrnm <- reads$mrnm
        mapq = reads$MQ
        bad.chr <- !(mrnm %in% c(seq(22), 'X', 'Y', 'M'))
        mrnm[bad.chr] <- reads$seqnames[bad.chr]
    }

    if (inherits(reads, 'GappedAlignments'))
        mwidth = qwidth(reads)
    else
    {
        mwidth = reads$qwidth
        mwidth[is.na(mwidth)] = 0
    }

    mwidth[is.na(mpos)] = 0
    mwidth[bad.chr] = 0;  # we set mates with "bad" chromosomes to have 0 width
    mpos[is.na(mpos)] = 1;

    if (inherits(reads, 'GappedAlignments'))
        GRanges(mrnm, IRanges(mpos, width = mwidth), strand = c('+', '-')[1+bamflag(reads)[, 'isMateMinusStrand']], seqlengths = seqlengths(reads), qname = values(reads)$qname, mapq = mapq)
    else if (inherits(reads, 'GRanges'))
        GRanges(mrnm, IRanges(mpos, width = mwidth), strand = c('+', '-')[1+bamflag(reads$flag)[, 'isMateMinusStrand']], seqlengths = seqlengths(reads), qname = values(reads)$qname, mapq = mapq)
    else if (inherits(reads, 'data.table'))
        ab=data.table(seqnames=mrnm, start=mpos, end=mpos + mwidth - 1, strand=c('+','-')[1+bamflag(reads$flag)[,'isMateMinusStrand']], qname=reads$qname, mapq = mapq)
}

#' Takes reads object and returns grl with each read and its mate (if exists)
#'
#' @param reads \code{GRanges} holding reads
#' @param as.grl Default TRUE. Return as a \code{GRangesList}
#' @param verbose Default FALSE
#' @name get.pairs.grl
#' @export
get.pairs.grl = function(reads, as.grl = TRUE, verbose = F)
{

    isdt <- inherits(reads, 'data.table')

    bad.col = c("seqnames", "ranges", "strand", "seqlevels",
                "seqlengths", "isCircular", "genome", "start", "end", "width", "element")

    if (verbose)
        cat('deduping\n')

    if (is(reads, 'GappedAlignmentPairs'))
        reads = unlist(reads)

    if (inherits(reads, 'GRanges')) {
        d <- duplicated(values(reads)$qname) ## duplicates are already paired up
        qpair.ix <- values(reads)$qname %in% unique(values(reads)$qname[d])
    } else if (isdt) {
        d <- duplicated(reads$qname)
        qpair.ix <- reads$qname %in% unique(reads$qname[d])
    }

    if (!inherits(reads, 'GenomicRanges') && !inherits(reads, 'data.table'))
    {
        if (verbose)
            cat('converting to granges\n')
        r.gr = granges(reads)
    }
    else if (!isdt)
        r.gr = reads[, c()]
    else
        r.gr <- reads

    if (verbose)
        cat('grbinding\n')

    m.gr = get.mate.gr(reads[!qpair.ix]);

    if (inherits(reads, 'GRanges')) {
        m.val <- values(m.gr)
        values(m.gr) = NULL;
        r.gr = c(r.gr, m.gr);
        mcols(r.gr) <- rrbind2(mcols(reads)[, setdiff(colnames(values(reads)), bad.col)], m.val)
    } else if (isdt) {
        m.gr <- m.gr[, setdiff(colnames(reads), colnames(m.gr)) := NA, with=FALSE]
        r.gr <- rbind(reads, m.gr, use.names=TRUE)
        setkey(r.gr, qname)
    }

    if (as.grl && !isdt) {
        if (verbose)
            cat('splitting\n')
        return(split(r.gr, as.character(r.gr$qname)))
    } else {
        return(r.gr)
    }
}

#' count.clips
#'
#' takes gr or gappedalignment object and uses cigar field (or takes character vector of cigar strings)
#' and returns data frame with fields (for character input)
#' $right.clips number of "right" soft clips (eg cigar 89M12S)
#' #left.clips number of "left" soft clips (eg cigar 12S89M)
#' or appends these fields to the reads object
#'
#' @param reads GenomicRanges holding the reads
#' @param hard [Default TRUE] option counts hard clips
#' @name count.clips
#' @export
count.clips = function(reads, hard = FALSE)
{
    if (length(reads) == 0)
        return(reads)
    if (inherits(reads, 'GRanges') | inherits(reads, 'GappedAlignments'))
        cigar = values(reads)$cigar
    else
        cigar = reads;

    if (!inherits(cigar, 'character') & !inherits(cigar, 'factor'))
        stop('Input must be GRanges, GappedAlignments, or character vector')

    out = data.frame(left.clips = rep(0, length(cigar)), right.clips = rep(0, length(cigar)));

    re.left = '^(\\d+)S.*'
    re.right = '.*[A-Z](\\d+)S$'

    lclip.ix = grep(re.left, cigar)
    rclip.ix = grep(re.right, cigar)

    if (length(lclip.ix)>0)
    {
        left.clips = gsub(re.left, '\\1', cigar[lclip.ix])
        out$left.clips[lclip.ix] = as.numeric(left.clips)
    }

    if (length(rclip.ix)>0)
    {
        right.clips = gsub(re.right, '\\1', cigar[rclip.ix])
        out$right.clips[rclip.ix] = as.numeric(right.clips)
    }

    if (inherits(reads, 'GRanges') | inherits(reads, 'GappedAlignments'))
    {
        values(reads)$right.clips = out$right.clips
        values(reads)$left.clips = out$left.clips
        out = reads
    }

    return(out)
}

#' varbase
#'
#' takes gr or gappedalignment object "reads" and uses cigar, MD, seq fields
#' to return variant bases and ranges
#'
#' returns grl (of same length as input) of variant base positions with character vector $varbase field populated with variant bases
#' for each gr item in grl[[k]], with the following handling for insertions, deletions, and substitution gr's:
#'
#' substitutions: nchar(gr$varbase) = width(gr) of the corresponding var
#' insertions: nchar(gr$varbase)>=1, width(gr) ==0
#' deletions: gr$varbase = '', width(gr)>=1
#'
#' Each gr also has $type flag which shows the cigar string code for the event ie
#' S = soft clip --> varbase represents clipped bases
#' I = insertion --> varbase represents inserted bases
#' D = deletion --> varbase is empty
#' X = mismatch --> varbase represents mismatched bases
#' @param reads GenomicRanges to extract variants from
#' @param soft [Default TRUE]
#' @param verbose [Default TRUE]
#' @name varbase
#' @export
varbase = function(reads, soft = TRUE, verbose = TRUE)
{
    nreads = length(reads)
    if (inherits(reads, 'GRangesList'))
    {
        was.grl = TRUE
        r.id = as.data.frame(reads)$element
        reads = unlist(reads)
    }
    else if (inherits(reads, 'data.frame'))
    {
        r.id = 1:nrow(reads)
        nreads = nrow(reads)
        was.grl = FALSE
    }
    else
    {
        r.id = 1:length(reads)
        was.grl = FALSE
    }

    if (!inherits(reads, 'GRanges') & !inherits(reads, 'GappedAlignments') & !inherits(reads, 'data.frame'))
        stop('Reads must be either GRanges, GRangesList, or GappedAlignment object')
    else if (inherits(reads, 'data.frame'))
    {
        if (is.null(reads$cigar) | is.null(reads$seq))
            stop('Reads must have cigar and seq fields specified')
    }
    else if (is.null(values(reads)$cigar) | is.null(values(reads)$seq))
        stop('Reads must have cigar and seq fields specified')

    if (is.data.frame(reads))
    {
        sl = NULL
        sn =  reads$seqnames
        cigar = as.character(reads$cigar)
        seq = as.character(reads$seq)
        str = reads$strand

        if (!is.null(reads$MD))
            md = as.character(reads$MD)
        else
            md = rep(NA, length(cigar))
    }
    else
    {
        sl = seqlengths(reads)
        sn =  seqnames(reads)
        cigar = as.character(values(reads)$cigar)
        seq = as.character(values(reads)$seq)
        str = as.character(strand(reads))

        if (!is.null(values(reads)$MD))
            md = as.character(values(reads)$MD)
        else
            md = rep(NA, length(cigar))
    }

    if (!inherits(cigar, 'character') & !inherits(cigar, 'character') & !inherits(md, 'character'))
        stop('Input must be GRanges with seq, cigar, and MD fields populated or GappedAlignments object')

    ix = which(!is.na(cigar))

    if (length(ix)==0)
        return(rep(GRangesList(GRanges()), nreads))

    cigar = cigar[ix]
    seq = seq[ix]
    md = md[ix]
    str = str[ix]

    if (is.data.frame(reads))
    {
        r.start = reads$start[ix]
        r.end = reads$end[ix]
    }
    else
    {
        r.start = start(reads)[ix]
        r.end = end(reads)[ix];
    }

    flip = str == '-'

    seq = strsplit(seq, '')

    cigar.vals = lapply(strsplit(cigar, "\\d+"), function(x) x[2:length(x)])
    cigar.lens = lapply(strsplit(cigar, "[A-Z]"), as.numeric)

    clip.left = sapply(cigar.vals, function(x) x[1] == 'S')
    clip.right = sapply(cigar.vals, function(x) x[length(x)] == 'S')

    if (any(clip.left))
        r.start[clip.left] = r.start[clip.left]-sapply(cigar.lens[which(clip.left)], function(x) x[1])

    if (any(clip.right))
        r.end[clip.right] = r.end[clip.right]+sapply(cigar.lens[which(clip.right)], function(x) x[length(x)])

                                        # split md string into chars after removing "deletion" signatures and also
                                        # any soft clipped base calls (bases followed by a 0
    md.vals = strsplit(gsub('([ATGC])', '|\\1|', gsub('\\^[ATGC]+', '|', md)), '\\|')

                                        # ranges of different cigar elements relative to query ie read-centric coordinates
    starts.seq = lapply(1:length(cigar.lens), function(i)
    {
        x = c(0, cigar.lens[[i]])
        x[which(cigar.vals[[i]] %in% c('D', 'H', 'N'))+1] = 0  ## deletions have 0 width on query
        cumsum(x[1:(length(x)-1)])+1
    })

    ends.seq = lapply(1:length(cigar.lens), function(i)
    {
        x = cigar.lens[[i]];
        x[which(cigar.vals[[i]] %in% c('D', 'H', 'N'))] = 0
        cumsum(x)
    })

                                        # ranges of different cigar elements relative to reference coordinatse
    starts.ref = lapply(1:length(cigar.lens), function(i)
    {
        x = c(0, cigar.lens[[i]]);
        x[which(cigar.vals[[i]] %in% c('I'))+1] = 0 ## insertions have 0 width on reference / subject
        cumsum(x[1:(length(x)-1)]) + r.start[i]
    })

    ends.ref = lapply(1:length(cigar.lens), function(i)
    {
        x = cigar.lens[[i]];
        x[which(cigar.vals[[i]] %in% c('I'))] = 0
        cumsum(x) + r.start[i] - 1
    })

                                        # now using MD find coordinates of mismatched bases (using starts and ends of M regions)

                                        # find coords of subs on genome
    tmp.pos = lapply(1:length(md.vals), function(i)
    {
        x = md.vals[[i]]
        nix = grepl('[ATGC]', x);
        if (!any(nix))
            return(c())
        p = rep(0, length(x))
        p[!nix] = as.numeric(x[!nix])
        p[nix] = 1
        s.pos.m = cumsum(p)[nix] ## position of subs in read
        mix = cigar.vals[[i]]=='M'
        m.st = cumsum(c(1, ends.seq[[i]][mix]-starts.seq[[i]][mix]+1))
        m.st.g = starts.ref[[i]][mix]
        m.st.r = starts.seq[[i]][mix]
        s.match = rep(NA, length(s.pos.m)) ## indices of matching "M" cigar element for each sub
        for (ii in 1:length(s.pos.m))
        {
            j = 0;
            done = FALSE
            for (j in 0:(length(m.st)-1))
                if (s.pos.m[ii] < m.st[j+1])
                    break
            s.match[ii] = j
        }

        s.pos.g = m.st.g[s.match] + s.pos.m-m.st[s.match]
        s.pos.r = m.st.r[s.match] + s.pos.m-m.st[s.match]

        return(rbind(s.pos.g, s.pos.r))
    })
    subs.pos = lapply(tmp.pos, function(x) x[1,])
    subs.rpos = lapply(tmp.pos, function(x) x[2,])
                                        #  subs.base = lapply(md.vals, grep, pattern = '[ATGC]', value = T)
    subs.base = lapply(1:length(seq), function(x) seq[[x]][subs.rpos[[x]]])

                                        # make sure md and cigar are consistent
                                        # (for some reason - sometimes soft clipped mismatches are included in MD leading to a longer MD string)
                                        # also some MD are NA
    mlen.cigar = sapply(1:length(ends.seq), function(x) {mix = cigar.vals[[x]]=='M'; sum(ends.seq[[x]][mix]-starts.seq[[x]][mix]+1)})
    mlen.md = sapply(md.vals, function(x) {ix = grepl('[ATGC]', x); sum(as.numeric(x[!ix])) + sum(nchar(x[ix]))})
    good.md = which(!is.na(md))
                                        #  good.md = which(mlen.md == mlen.cigar & !is.na(md))

    if (any(na <- is.na(md)))
    {
        warning('MD field absent from one or more input reads')
        good.md = which(!na)
    }

                                        # now make granges of subs
    if (length(good.md)>0)
    {
        if (any(mlen.md[good.md] != mlen.cigar[good.md]))
            warning('the lengths of some MD strings do not match the number of M positions on the corresponding CIGAR string: some variants may not be correctly mapped to the genome')

        iix.md = unlist(lapply(good.md, function(x) rep(x, length(subs.pos[[x]]))))
        tmp = unlist(subs.pos[good.md])
        if (!is.null(tmp))
        {
            subs.gr = GRanges(sn[ix][iix.md], IRanges(tmp, tmp), strand = '*')
            values(subs.gr)$varbase = unlist(subs.base[good.md])
            values(subs.gr)$type = 'X'
            values(subs.gr)$iix = ix[iix.md]
        }
        else
            subs.gr = GRanges()
    }
    else
        subs.gr = GRanges()

    iix = unlist(lapply(1:length(cigar.vals), function(x) rep(x, length(cigar.vals[[x]]))))
    cigar.vals = unlist(cigar.vals)
    cigar.lens = unlist(cigar.lens)
    starts.seq = unlist(starts.seq)
    ends.seq = unlist(ends.seq)
    starts.ref = unlist(starts.ref)
    ends.ref = unlist(ends.ref)

    ## pick up other variants (including soft clipped and indel)
    is.var = cigar.vals != 'M'
    iix = iix[is.var]
    cigar.vals = cigar.vals[is.var]
    cigar.lens = cigar.lens[is.var]
    starts.ref = starts.ref[is.var]
    ends.ref = ends.ref[is.var]
    starts.seq = starts.seq[is.var]
    ends.seq = ends.seq[is.var]
    str <- str[iix] # JEREMIAH

    if (length(cigar.vals)>0)
    {
        var.seq = lapply(1:length(cigar.vals),
                         function(i)
                         {
                             if (ends.seq[i]<starts.seq[i])
                                 return('') # deletion
                             else
                                 seq[[iix[i]]][starts.seq[i]:ends.seq[i]] #insertion
                         })
        other.gr = GRanges(sn[ix][iix], IRanges(starts.ref, ends.ref), strand = str, seqlengths = sl)
        values(other.gr)$varbase = sapply(var.seq, paste, collapse = '')
        values(other.gr)$type = cigar.vals
        values(other.gr)$iix = ix[iix];

        out.gr = sort(c(subs.gr, other.gr))
    }
    else
        out.gr = subs.gr


                                        # add default colors to out.gr
    VAR.COL = get.varcol()

    col.sig = as.character(out.gr$type)
    xix = out.gr$type == 'X'
    col.sig[xix] = paste(col.sig[xix], out.gr$varbase[xix], sep = '')
    out.gr$col = VAR.COL[col.sig]
    out.gr$border = out.gr$col

    if (!soft)
        if (any(soft.ix <<- out.gr$type == 'S'))
            out.gr = out.gr[-which(soft.ix)]

    out.grl = rep(GRangesList(GRanges()), nreads)
    out.iix = r.id[values(out.gr)$iix]
    values(out.gr)$iix = NULL

    tmp.grl = split(out.gr, out.iix)
    out.grl[as.numeric(names(tmp.grl))] = tmp.grl
    values(out.grl)$qname[r.id] = reads$qname

    return(out.grl)
}

#' splice.cigar
#'
#' takes gr or gappedalignment object "reads" and parses cigar fields
#' to return grl corresponding to spliced alignments on the genome corresponding to
#' portions of the cigar
#'
#' ie each outputted grl item contains the granges corresponding to all non-N portions of cigar string
#'
#' if grl provided as input (e.g. paired reads) then all of the spliced ranges resulting from each
#' input grl item will be put into the corresponding output grl item
#'
#' NOTE: does not update MD tag
#'
#' if use.D = TRUE, then will treat "D" (deletion) in addition to "N" flags as indicative of deletion event.
#' @param reads \code{GRanges} reads
#' @param verbose Default TRUE
#' @param fast Default TRUE
#' @param use.D Default TRUE
#' @param rem.soft Default TRUE
#' @param get.seq Default FALSE
#' @param return.grl Default TRUE
#' @name splice.cigar
#' @export
splice.cigar = function(reads, verbose = TRUE, fast = TRUE, use.D = TRUE, rem.soft = TRUE, get.seq = FALSE, return.grl = TRUE)
{
    nreads = length(reads)

    if (nreads==0)
        if (return.grl)
            return(GRangesList())
        else
            return(GRanges)

    if (inherits(reads, 'GRangesList'))
    {
        was.grl = TRUE
        r.id = as.data.frame(reads)$element
        reads = unlist(reads)
    }
    else
    {
        r.id = 1:length(reads)
        was.grl = FALSE
    }


    if (is.data.frame(reads))
    {
        sl = NULL
        sn =  reads$seqnames
        cigar = as.character(reads$cigar)
        seq = as.character(reads$seq)
        str = reads$strand

        if (!is.null(reads$MD))
            md = as.character(reads$MD)
        else
            md = rep(NA, length(cigar))
    }
    else
    {
        sl = seqlengths(reads)
        sn =  seqnames(reads)
        cigar = as.character(values(reads)$cigar)
        seq = as.character(values(reads)$seq)
        str = as.character(strand(reads))

        if (!is.null(values(reads)$MD))
            md = as.character(values(reads)$MD)
        else
            md = rep(NA, length(cigar))
    }


    if (!inherits(reads, 'GRanges') & !inherits(reads, 'GappedAlignments') & !inherits(reads, 'data.frame'))
        stop('Reads must be either GRanges, GRangesList, or GappedAlignment object')
    else if (is.null(values(reads)$cigar) | is.null(values(reads)$seq))
        stop('Reads must have cigar and seq fields specified')

    if (!inherits(cigar, 'character') & !inherits(cigar, 'character') & !inherits(md, 'character'))
        stop('Input must be GRanges with seq, cigar, and MD fields populated or GappedAlignments object')


    ix = which(!is.na(reads$cigar))

    if (length(ix)==0)
        if (return.grl)
            return(rep(GRangesList(GRanges()), nreads))
        else
            return(GRanges())

    if (fast)
    {
        ir = cigarRangesAlongReferenceSpace(reads[ix]$cigar, N.regions.removed = FALSE, with.ops = TRUE, reduce.ranges = FALSE)
        irul = unlist(ir)
        out.gr = GRanges(rep(seqnames(reads)[ix], elementLengths(ir)), shift(IRanges(irul), rep(start(reads)[ix]-1, elementLengths(ir))),
                         strand = rep(strand(reads)[ix], elementLengths(ir)), seqlengths = seqlengths(reads))
        out.gr$type = names(irul)
        out.gr$rid = ix[rep(1:length(ir), elementLengths(ir))]
        out.gr$riid = unlist(lapply(elementLengths(ir), function(x) 1:x))
        out.gr$fid = r.id[out.gr$rid]
        out.gr$qname = reads$qname[out.gr$rid]

        if (return.grl)
        {
            out.grl = rep(GRangesList(GRanges()), nreads)
            tmp.grl = split(out.gr, out.gr$fid)
            out.grl[as.numeric(names(tmp.grl))] = tmp.grl
            return(out.grl)
        }
        else
            return(out.gr)
    }
    else
    {

        cigar = cigar[ix]
        str = str[ix]

        if (is.data.frame(reads))
        {
            r.start = reads$start[ix]
            r.end = reads$end[ix]
        }
        else
        {
            r.start = start(reads)[ix]
            r.end = end(reads)[ix];
        }

        flip = str == '-'
        cigar.vals = lapply(strsplit(cigar, "\\d+"), function(x) x[2:length(x)])
        cigar.lens = lapply(strsplit(cigar, "[A-Z]"), as.numeric)

        clip.left = sapply(cigar.vals, function(x) x[1] == 'S')
        clip.right = sapply(cigar.vals, function(x) x[length(x)] == 'S')

                                        # ranges of different cigar elements relative to query ie read-centric coordinates
        starts.seq = lapply(1:length(cigar.lens), function(i)
        {
            x = c(0, cigar.lens[[i]])
            x[which(cigar.vals[[i]] %in% c('D', 'H', 'N'))+1] = 0  ## deletions have 0 width on query
            cumsum(x[1:(length(x)-1)])+1
        })

        ends.seq = lapply(1:length(cigar.lens), function(i)
        {
            x = cigar.lens[[i]];
            x[which(cigar.vals[[i]] %in% c('D', 'H', 'N'))] = 0
            cumsum(x)
        })

                                        # ranges of different cigar elements relative to reference coordinatse
        starts.ref = lapply(1:length(cigar.lens), function(i)
        {
            x = c(0, cigar.lens[[i]]);
            x[which(cigar.vals[[i]] %in% c('I'))+1] = 0 ## insertions have 0 width on reference / subject
            cumsum(x[1:(length(x)-1)]) + r.start[i]
        })

        ends.ref = lapply(1:length(cigar.lens), function(i)
        {
            x = cigar.lens[[i]];
            x[which(cigar.vals[[i]] %in% c('I'))] = 0
            Cumsum(x) + r.start[i] - 1
        })

        iix = unlist(lapply(1:length(cigar.vals), function(x) rep(x, length(cigar.vals[[x]]))))
        cigar.vals = unlist(cigar.vals)
        cigar.lens = unlist(cigar.lens)
        starts.seq = unlist(starts.seq)
        ends.seq = unlist(ends.seq)
        starts.ref = unlist(starts.ref)
        ends.ref = unlist(ends.ref)

        ## pick up splice (including soft clipped and indel)
        splice.char = 'N'

        if (use.D)
            splice.char = c(splice.char, 'D')

        if (rem.soft)
            splice.char = c(splice.char, 'S')

        is.splice = !(cigar.vals %in% splice.char)
        iix = iix[is.splice]
        cigar.vals = cigar.vals[is.splice]
        cigar.lens = cigar.lens[is.splice]
        starts.ref = starts.ref[is.splice]
        ends.ref = ends.ref[is.splice]
        starts.seq = starts.seq[is.splice]
        ends.seq = ends.seq[is.splice]
        str <- str[iix] #

        out.gr = GRanges()
        if (length(cigar.vals)>0)
        {
            other.gr = GRanges(sn[ix][iix], IRanges(starts.ref, ends.ref), strand = str, seqlengths = sl)

            if (get.seq)
            {
                var.seq = lapply(1:length(cigar.vals),
                                 function(i)
                                 {
                                     if (ends.seq[i]<starts.seq[i])
                                         return('') # deletion
                                     else
                                         seq[[iix[i]]][starts.seq[i]:ends.seq[i]] #insertion
                                 })
                values(other.gr)$seq = sapply(var.seq, paste, collapse = '')
            }

            values(other.gr)$type = cigar.vals
            values(other.gr)$iix = ix[iix];

            out.gr = c(other.gr, out.gr)
        }

        out.gr$rid = out.gr$iix
        out.iix = r.id[out.gr$rid]
        values(out.gr)$iix = NULL
        out.gr$fid = out.iix
        out.gr$qname = reads$qname[out.gr$rid]

        if (return.grl)
        {
            out.grl = rep(GRangesList(GRanges()), nreads)
            tmp.grl = split(out.gr, out.iix)
            out.grl[as.numeric(names(tmp.grl))] = tmp.grl
            return(out.grl)
        }
        else
            return(out.gr)
    }
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


#' Improved rbind for intersecting columns of data.frames or data.tables
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


#' Provide transparency to colors
#'
#' takes provided colors and gives them the specified alpha (ie transparency) value
#'
#' @param col color string
#' @param alpha alpha value between 0 and 1
#' @return rgb color like the input, but with transparency added
#' @export
alpha = function(col, alpha)
{
    col.rgb = col2rgb(col)
    return(rgb(red = col.rgb['red', ]/255, green = col.rgb['green', ]/255, blue = col.rgb['blue', ]/255, alpha = alpha))
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

#' More robust and faster implementation of GenomicRangs::setdiff
#'
#' Robust to common edge cases of setdiff(gr1, gr2)  where gr2 ranges are contained inside gr1's (yieldings
#' setdiffs yield two output ranges for some of the input gr1 intervals.
#'
#' @param query \code{GRanges} object as query
#' @param subject \code{GRanges} object as subject
#' @param max.slice Default Inf. If query is bigger than this, chunk into smaller on different cores
#' @param verbose Default FALSE
#' @param mc.cores Default 1. Only works if exceeded max.slice
#' @param ... arguments to be passed to \link{gr.findoverlaps}
#' @return returns indices of query in subject or NA if none found
#' @name gr.match
#' @export
gr.setdiff = function(query, subject, ignore.strand = TRUE, by = NULL,  ...)
{
    if (!is.null(by)) ## in this case need to be careful about setdiffing only within the "by" level
    {
        tmp = grdt(subject)
        tmp$strand = factor(tmp$strand, c('+', '-', '*'))
        sl = seqlengths(subject)
        gp = seg2gr(tmp[, as.data.frame(gaps(IRanges(start, end), 1, sl[seqnames][1])), by = c('seqnames', 'strand', by)], seqinfo = seqinfo(subject))
    }
    else ## otherwise easier
    {
        if (ignore.strand)
            gp = gaps(gr.stripstrand(subject)) %Q% (strand == '*')
        else
            gp = gaps(subject)
    }

    out = gr.findoverlaps(query, gp, qcol = names(values(query)), ignore.strand = ignore.strand, by = by, ...)
    return(out)
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

