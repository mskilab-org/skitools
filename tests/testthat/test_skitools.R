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

## NOTE: the old "test ra.merge" block that lived here was dropped (mskilab-org/skitools#11).
## ra.merge was moved out of skitools into gGnome in commit 1d04255 ("removed ra.merge
## (duplicated with gGnome::ra.merge)"), so this package no longer defines it, and the block
## could not pass. It was stale in a second way too: it referenced grl1 / grl2, which are
## never defined in this file, and its expectations index values(ram)[, 1:2], which under
## gGnome::ra.merge are tmp.ix / merged.ix, not the seen.by.ra* flags. Any test of ra.merge
## belongs in gGnome, which now owns the function.
