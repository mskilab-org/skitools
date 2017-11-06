library(skitools)

library(BSgenome.Hsapiens.UCSC.hg19)
Sys.setenv(DEFAULT_BSGENOME = "BSgenome.Hsapiens.UCSC.hg19::Hsapiens")


context("skitool functions")


test_that("test file.name()", {
    
    example_paths = c("/user/olivia/homedirectory/special_project_42/file.bam", "short/path/example.txt", "usr/bin/projects/special_Project3/example.vcf")
    expect_identical(file.name(example_paths), c("file.bam", "example.txt", "example.vcf"))

})

test_that("test file.dir()", {
    
    example_paths = c("/user/olivia/homedirectory/special_project_42/file.bam", "short/path/example.txt", "usr/bin/projects/special_Project3/example.vcf")
    expect_identical(file.dir(example_paths[1]), "/user/olivia/homedirectory/special_project_42/")
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







