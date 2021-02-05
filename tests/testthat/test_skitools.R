library(skitools)
library(gUtils)

library(BSgenome.Hsapiens.UCSC.hg19)
Sys.setenv(DEFAULT_BSGENOME = "BSgenome.Hsapiens.UCSC.hg19::Hsapiens")


context("skitool functions")


test_that("test file.name()", {
    
    example_paths = c("/user/oliviafoofoo/homedirectory/special_project_42/file.bam", "short/path/example.txt", "usr/bin/projects/special_Project3/example.vcf")
    expect_identical(file.name(example_paths), c("file.bam", "example.txt", "example.vcf"))

})

test_that("test file.dir()", {
    
    example_paths = c("/user/oliviafoofoo/homedirectory/special_project_42/file.bam", "short/path/example.txt", "usr/bin/projects/special_Project3/example.vcf")
    expect_identical(file.dir(example_paths[1]), "/user/oliviafoofoo/homedirectory/special_project_42/")
    expect_identical(file.dir(example_paths[2]), "short/path/")
    expect_identical(file.dir(example_paths[3]), "usr/bin/projects/special_Project3/")
    
})


test_that("test gstring()", {
    
    example_igv= "17:7661779-7687550;17:43044295-43170245;12:25204789-25250936"
    expect_equal(width(gstring(example_igv)[1]), 25772)
    expect_equal(width(gstring(example_igv)[2]), 125951)
    expect_equal(width(gstring(example_igv)[3]), 46148)
    
})


test_that("test is.dup()", {
    
    expect_equal(any(is.dup(c(1:10))), FALSE)
    expect_equal(is.dup(c("A", "A", "A", "B", "C", "C")), c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE))
    
})



test_that("test charToDec()", {
    
    expect_equal(charToDec("BRCA2"), c(66, 82, 67, 65, 50))
    expect_equal(charToDec("$Foo.bar!"), c(36, 70, 111, 111, 46, 98, 97, 114, 33))
    
})





test_that('test ra.merge', {

    ## beginning with fake rearrangment data grl1 and grl2
    expect_equal(length(ra.merge(grl1, grl2)), 500)
    expect_true(unique(values(ra.merge(grl1, grl2))[, 1][1:250]))
    expect_false(unique(values(ra.merge(grl1, grl2))[, 1][251:500]))
    expect_false(unique(values(ra.merge(grl1, grl2))[, 2][1:249]))
    expect_true(unique(values(ra.merge(grl1, grl2))[, 2][250:500]))
    ## ignore.strand == TRUE makes no difference...
    expect_equal(ra.merge(grl1, grl2, ignore.strand=TRUE), ra.merge(grl1, grl2))
    ## example in function 'ra.merge()'
    gr1 = GRanges(1, IRanges(1:10, width = 1), strand = rep(c('+', '-'), 5))
    gr2 = GRanges(1, IRanges(4 + 1:10, width = 1), strand = rep(c('+', '-'), 5))
    expect_error(ra.merge(gr1, gr2)) ##   Error: All inputs must be a GRangesList
    ## create GRangesLists
    ra1 = split(gr1, rep(1:5, each = 2))
    ra2 = split(gr2, rep(1:5, each = 2))
    ram = ra.merge(ra1, ra2)
    expect_warning(ra.merge(ra1, ra2))  ## warning: GRanges object contains 10 out-of-bound ranges located on sequence 1.
    expect_equal(length(ram), 7)
    ## 'values(ram)' shows the metadata with TRUE / FALSE flags
    expect_equal(values(ram)[, 1], c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)) 
    expect_equal(values(ram)[, 2], c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE)) 
    ## ram2 = ra.merge(ra1, ra2, pad = 5) # more inexact matching results in more merging
    ram2 = ra.merge(ra1, ra2, pad = 500) ## more inexact matching results in more merging
    ##values(ram2)
    expect_equal(length(ram2), 5)    
    expect_equal(values(ram2)[, 1], c(TRUE, TRUE, TRUE, TRUE, TRUE)) 
    expect_equal(values(ram2)[, 2], c(TRUE, TRUE, TRUE, TRUE, TRUE)) 
    expect_error(ra.merge(ra1, ra2, pad =-1)) ## adjustment would result in ranges with negative widths
    ## ram3
    ram3 = ra.merge(ra1, ra2, ind = TRUE) ## indices instead of flags
    ##values(ram3)
    expect_equal(length(ram3), 7)
    expect_equal(values(ram3)[, 1], c(1, 2, 3, 4, 5, NA, NA)) 
    expect_equal(values(ram3)[, 2], c(NA, NA, 3, 4, 5, 4, 5)) 
    ## test both 'pad', 'ind'
    ram4 = ra.merge(ra1, ra2, pad = 500, ind = TRUE) 
    expect_equal(values(ram4)[, 1], c(1, 2, 3, 4, 5))
    expect_equal(values(ram4)[, 2], c(1, 2, 3, 4, 5))
    ## ignore.strand == TRUE
    expect_error(ra.merge(ra1, ra2, ignore.stand = TRUE)) ##  unable to find an inherited method for function ‘values’ for signature ‘"logical"’
    ### all args
    expect_error(ra.merge(ra1, ra2, pad = 500, ind = TRUE, ignore.stand = TRUE)) ## unable to find an inherited method for function ‘values’ for signature ‘"logical"’
 
})




