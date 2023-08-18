blackspotted-rougheye-SST
================
Diana Baetscher
2023-08-18

Walkthrough of species id for unknown samples:

Using the rds file output from microhaplot: 1. read in rds files 2.
apply read depth filters 3. apply allele balance filter

``` r
library(tidyverse)
```

    ## -- Attaching packages --------------------------------------- tidyverse 1.3.1 --

    ## v ggplot2 3.4.1      v purrr   0.3.4 
    ## v tibble  3.1.2      v dplyr   1.0.10
    ## v tidyr   1.2.0      v stringr 1.4.0 
    ## v readr   1.4.0      v forcats 0.5.1

    ## Warning: package 'ggplot2' was built under R version 4.1.3

    ## Warning: package 'tidyr' was built under R version 4.1.3

    ## Warning: package 'dplyr' was built under R version 4.1.3

    ## -- Conflicts ------------------------------------------ tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(readxl)
library(stringr)
library(lubridate)
```

    ## 
    ## Attaching package: 'lubridate'

    ## The following objects are masked from 'package:base':
    ## 
    ##     date, intersect, setdiff, union

``` r
library(rubias)
```

    ## Warning: package 'rubias' was built under R version 4.1.3

``` r
library(ggpattern)
```

    ## Warning: package 'ggpattern' was built under R version 4.1.3

``` r
source("../R/rockfish-funcs2.R")


#### Call genos from the microhaplot rds files ####

# the directory with the rds file
dir <- "../data_AFSC/rds_files/"


# cycle over them, read them and add the gtseq_run column on each.
# at the end, bind them together.
genos_long <- call_genos_from_haplotRDS(path = file.path(dir, "blackspot_rougheye_SST.rds"))
```

    ## Warning: `tbl_df()` was deprecated in dplyr 1.0.0.

    ## Warning: Please use `tibble::as_tibble()` instead.

    ## Joining, by = c("id", "locus", "rank")

``` r
#### In the end, let us get a data frame that includes genotypes for all the individuals  ####
# and which explicitly has NAs in places where data are missing
genos_long_explicit_NAs <- genos_long %>%
  select(id) %>%
  unique() %>%
  unlist() %>%
  unname() %>%
  expand.grid(id = ., locus = unique(genos_long$locus), gene_copy = 1:2, stringsAsFactors = FALSE) %>%
  tbl_df() %>% 
  left_join(., genos_long) %>%
  arrange(id, locus, gene_copy)
```

    ## Joining, by = c("id", "locus", "gene_copy")

``` r
genos_long_explicit_NAs %>%
  ggplot(aes(x = locus, y = id, fill = depth)) +
  geom_tile()
```

![](08-blackspotted-rougheye-SST_files/figure-gfm/overall-pattern-of-missing-data-1.png)<!-- -->
There are five loci that drop-out in blackspotted and rougheye (and
shortspine thornyhead).

Using those genotypes…

``` r
genos_long_explicit_NAs %>%
  group_by(id) %>%
  tally()
```

    ## # A tibble: 287 x 2
    ##    id        n
    ##    <chr> <int>
    ##  1 s1      180
    ##  2 s10     180
    ##  3 s100    180
    ##  4 s101    180
    ##  5 s102    180
    ##  6 s103    180
    ##  7 s104    180
    ##  8 s105    180
    ##  9 s106    180
    ## 10 s107    180
    ## # ... with 277 more rows

287 samples

Look at missing data: 180 gene copies total (90 loci x2)

``` r
ind_to_toss <- genos_long_explicit_NAs %>%
  group_by(id) %>%
  filter(is.na(allele)) %>% # missing data
  tally() %>%
  arrange(desc(n)) %>% # remove samples with >20% missing data
  filter(n > 36) 

# remove those from the df
genos_ind_filtered <- genos_long_explicit_NAs %>%
  anti_join(., ind_to_toss)
```

    ## Joining, by = "id"

105 samples removed due to \>20% missing data.

``` r
ind_to_toss %>%
  inner_join(., genos_long_explicit_NAs) %>%
  ggplot(aes(x = locus, y = id, fill = depth)) +
  geom_tile()
```

    ## Joining, by = "id"

![](08-blackspotted-rougheye-SST_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

Load baseline data

``` r
# baseline data - curated, 997 indivs
baseline <- readRDS("../new_baseline_data/processed/sebastes_spp_id_baseline_haplotypes.rds")

# remove the 6 loci that had HWE and other issues
to_remove <- read_csv("../data/loci_to_remove.csv")
```

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   locus = col_character()
    ## )

``` r
baseline90 <- baseline %>%
  anti_join(., to_remove)
```

    ## Joining, by = "locus"

``` r
# remind myself which species are in the baseline:
baseline90 %>%
  select(collection) %>%
  unique() %>%
  arrange()
```

    ## # A tibble: 54 x 1
    ##    collection   
    ##    <chr>        
    ##  1 aleutianus   
    ##  2 alutus       
    ##  3 auriculatus  
    ##  4 aurora       
    ##  5 babcocki     
    ##  6 borealis     
    ##  7 caurinus     
    ##  8 chlorostictus
    ##  9 constellatus 
    ## 10 crameri      
    ## # ... with 44 more rows

``` r
tossers <- baseline90 %>%
  select(indiv, gtseq_run, id) %>%
  unique() %>%
  group_by(indiv) %>%
  tally() %>%
  filter(n >1)

baseline90_one_each <- baseline90 %>%
  anti_join(., tossers)
```

    ## Joining, by = "indiv"

``` r
# baseline data - curated, 997 indivs
baseline_spp_info <- baseline90_one_each %>%
  select(sample_type, repunit, collection, indiv, gtseq_run, id, species) %>%
  unique()
baseline_spp_info$gtseq_run <- as.character(baseline_spp_info$gtseq_run)
```

``` r
# slim that down to just the matching field with the unknowns
for_alleidx <- baseline90_one_each %>%
  select(-indiv, -c(1:3, 12:13), -species)
  

for_alleidx$gtseq_run <- as.character(for_alleidx$gtseq_run)
```

``` r
# merge the two dataframes
merged_df <- bind_rows(for_alleidx, genos_ind_filtered)

# first make integers of the alleles
alle_idxs <- merged_df %>% 
  dplyr::select(gtseq_run, id, locus, gene_copy, allele) %>%
  group_by(locus) %>%
  mutate(alleidx = as.integer(factor(allele, levels = unique(allele)))) %>%
  ungroup() %>%
  arrange(gtseq_run, id, locus, alleidx) # rubias can handle NA's, so no need to change them to 0's

  
# and spread the alleles
two_col <- alle_idxs %>%
  #group_by(indiv, locus) %>%
  unite(loc, locus, gene_copy, sep = ".") %>%
  #ungroup() %>%
  select(-allele) %>%
  pivot_wider(names_from = loc, values_from = alleidx) 
```

add back on info for reference and make two-column format for rubias

``` r
# baseline
reference <- two_col %>%
  left_join(., baseline_spp_info) %>%
  filter(!is.na(species)) %>%
  select(-gtseq_run, -id, -species) %>%
  select(sample_type, repunit, collection, indiv, everything())
```

    ## Joining, by = c("gtseq_run", "id")

``` r
# mixture
rubias_mix <- two_col %>%
  anti_join(., baseline_spp_info) %>%
  select(-gtseq_run) %>%
  mutate(sample_type = "mixture", collection = "blackspotted_rougheye", repunit = NA) %>%
  select(sample_type, repunit, collection, everything()) %>%
  rename(indiv = id)
```

    ## Joining, by = c("gtseq_run", "id")

## Mixture assignment with rubias

``` r
rubias_output <- infer_mixture(reference = reference, mixture = rubias_mix, gen_start_col = 5)
```

    ## All gene copies missing at locus tag_id_1471.1 in mixture data. Ploidy indeterminate in that data frame.

    ## All allelic data missing for the following loci in the mixture. Ploidy inferred from reference. tag_id_1471.1

    ## Collating data; compiling reference allele frequencies, etc.   time: 1.66 seconds
    ## Computing reference locus specific means and variances for computing mixture z-scores   time: 0.17 seconds
    ## Working on mixture collection: blackspotted_rougheye with 182 individuals
    ##   calculating log-likelihoods of the mixture individuals.   time: 0.08 seconds
    ##   performing 2000 total sweeps, 100 of which are burn-in and will not be used in computing averages in method "MCMC"   time: 0.25 seconds
    ##   tidying output into a tibble.   time: 0.08 seconds

``` r
# mixing proportions
rubias_output$mixing_proportions %>%
  ggplot(aes(x = collection, y = pi)) +
  geom_bar(stat = "identity") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.95)
  )
```

![](08-blackspotted-rougheye-SST_files/figure-gfm/run-rubias-1.png)<!-- -->

``` r
# take the top output for each sample
top_assign <- rubias_output$indiv_posteriors %>%
  group_by(indiv) %>%
  slice_max(., order_by = PofZ)


top_assign 
```

    ## # A tibble: 182 x 10
    ## # Groups:   indiv [182]
    ##    mixture~1 indiv repunit colle~2  PofZ log_l~3 z_score n_non~4 n_mis~5 missi~6
    ##    <chr>     <chr> <chr>   <chr>   <dbl>   <dbl>   <dbl>   <int>   <int> <list> 
    ##  1 blackspo~ s1    melano~ melano~     1   -64.1  2.37        80      10 <int>  
    ##  2 blackspo~ s10   melano~ melano~     1   -94.7 -0.0824      83       7 <int>  
    ##  3 blackspo~ s100  melano~ melano~     1   -73.9  1.64        82       8 <int>  
    ##  4 blackspo~ s101  melano~ melano~     1   -75.4  1.59        82       8 <int>  
    ##  5 blackspo~ s104  aleuti~ aleuti~     1   -64.3  0.128       83       7 <int>  
    ##  6 blackspo~ s105  aleuti~ aleuti~     1   -68.7 -0.245       82       8 <int>  
    ##  7 blackspo~ s106  aleuti~ aleuti~     1   -45.3  1.93        81       9 <int>  
    ##  8 blackspo~ s107  aleuti~ aleuti~     1   -71.8 -1.07        80      10 <int>  
    ##  9 blackspo~ s108  aleuti~ aleuti~     1   -67.6 -0.548       81       9 <int>  
    ## 10 blackspo~ s109  aleuti~ aleuti~     1   -50.8  1.14        82       8 <int>  
    ## # ... with 172 more rows, and abbreviated variable names 1: mixture_collection,
    ## #   2: collection, 3: log_likelihood, 4: n_non_miss_loci, 5: n_miss_loci,
    ## #   6: missing_loci

``` r
top_assign %>%
  filter(z_score < -3 | z_score > 3)
```

    ## # A tibble: 1 x 10
    ## # Groups:   indiv [1]
    ##   mixture_~1 indiv repunit colle~2  PofZ log_l~3 z_score n_non~4 n_mis~5 missi~6
    ##   <chr>      <chr> <chr>   <chr>   <dbl>   <dbl>   <dbl>   <int>   <int> <list> 
    ## 1 blackspot~ s88   melano~ melano~  1.00   -153.   -4.92      82       8 <int>  
    ## # ... with abbreviated variable names 1: mixture_collection, 2: collection,
    ## #   3: log_likelihood, 4: n_non_miss_loci, 5: n_miss_loci, 6: missing_loci

``` r
top_assign_filtered <- top_assign %>%
  filter(z_score > -3 | z_score < 3) %>%
  filter(PofZ > 0.9) 

top_assign_filtered %>%
  group_by(repunit) %>%
  tally()
```

    ## # A tibble: 3 x 2
    ##   repunit           n
    ##   <chr>         <int>
    ## 1 aleutianus       83
    ## 2 aurora            1
    ## 3 melanostictus    98

Map those results:

``` r
samplesheet <- read_csv("../data_AFSC/samplesheets/20230807_rockfish_SST.csv", skip = 18)
```

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   Sample_ID = col_character(),
    ##   Sample_Plate = col_character(),
    ##   Sample_Well = col_character(),
    ##   I7_Index_ID = col_character(),
    ##   index = col_character(),
    ##   I5_Index_ID = col_logical(),
    ##   index2 = col_character(),
    ##   Description = col_logical(),
    ##   SampleProject = col_logical(),
    ##   id = col_character()
    ## )

``` r
meta <- read_csv("../data_AFSC/ABLG_database_20230724_export.csv")
```

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   .default = col_character(),
    ##   ABLG = col_double(),
    ##   CollectionYear = col_double(),
    ##   CollectionMonth = col_double(),
    ##   CollectionDay = col_double(),
    ##   StartLatitudeDD = col_double(),
    ##   StartLongitudeDD = col_double(),
    ##   EstimatedLatLongYesNo = col_logical(),
    ##   DepthM = col_double(),
    ##   StartDepthM = col_logical(),
    ##   EndDepthM = col_logical(),
    ##   TempC = col_double(),
    ##   LengthMM = col_double(),
    ##   LengthType = col_logical(),
    ##   FishDepthMM = col_logical(),
    ##   WeightG = col_double(),
    ##   CaptureTime = col_logical(),
    ##   ReleaseTime = col_logical(),
    ##   ExtractionMethod = col_logical(),
    ##   DNAngul = col_double(),
    ##   Comments = col_logical()
    ##   # ... with 12 more columns
    ## )
    ## i Use `spec()` for the full column specifications.

    ## Warning: 20160 parsing failures.
    ##  row      col           expected                   actual                                             file
    ## 1927 Comments 1/0/T/F/TRUE/FALSE date range - march 23-26 '../data_AFSC/ABLG_database_20230724_export.csv'
    ## 1928 Comments 1/0/T/F/TRUE/FALSE date range - march 23-26 '../data_AFSC/ABLG_database_20230724_export.csv'
    ## 1929 Comments 1/0/T/F/TRUE/FALSE date range - march 23-26 '../data_AFSC/ABLG_database_20230724_export.csv'
    ## 1930 Comments 1/0/T/F/TRUE/FALSE date range - march 23-26 '../data_AFSC/ABLG_database_20230724_export.csv'
    ## 1931 Comments 1/0/T/F/TRUE/FALSE date range - march 23-26 '../data_AFSC/ABLG_database_20230724_export.csv'
    ## .... ........ .................. ........................ ................................................
    ## See problems(...) for more details.

``` r
meta$ABLG <- as.character(meta$ABLG)

results_for_map <- samplesheet %>%
  select(Sample_ID, id) %>%
  left_join(., top_assign, by = c("id" = "indiv")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "ABLG", "")) %>%
  left_join(., meta, by = c("Sample_ID" = "ABLG"))
```

``` r
library(ggplot2)
library(sf)  
```

    ## Warning: package 'sf' was built under R version 4.1.3

    ## Linking to GEOS 3.10.2, GDAL 3.4.1, PROJ 7.2.1; sf_use_s2() is TRUE

``` r
library(rnaturalearth)
```

    ## Warning: package 'rnaturalearth' was built under R version 4.1.3

``` r
library(rnaturalearthdata)
```

    ## Warning: package 'rnaturalearthdata' was built under R version 4.1.3

    ## 
    ## Attaching package: 'rnaturalearthdata'

    ## The following object is masked from 'package:rnaturalearth':
    ## 
    ##     countries110

``` r
world <- ne_countries(scale = "medium", returnclass = "sf")

goa_map <- ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-140, -120), 
           ylim = c(55, 48), 
           expand = FALSE) 

goa_map +
  geom_jitter(data = results_for_map, aes(x = -1*StartLongitudeDD, y = StartLatitudeDD, color = repunit), size = 1, alpha = 0.5, width = 0.3, height = 0.3) +
  scale_color_manual(values = c("darkred", "gold", "salmon", "darkorange", "dodgerblue", "lightblue", "midnightblue")) +
  theme_bw() +
  labs(x = "Longitude",
       y = "Latitude",
       color = "Species (genetic)",
       shape = "Species (morphology)") +
  theme(
    legend.text = element_text(face = "italic"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
  ) +
  scale_y_continuous(breaks = c(40, 45, 50, 55, 60, 65)) +
  scale_x_continuous(breaks = c(-180, -160, -150, -140, -130, -120)) +
  guides(color = guide_legend(
    override.aes = list(alpha = 1, size = 3), 
    label.theme = element_text(size = 12, face = "italic"),
    title.theme = element_text(size = 14)),
    shape = guide_legend(override.aes = list(alpha = 1, size = 3)))
```

    ## Warning: Removed 98 rows containing missing values (`geom_point()`).

![](08-blackspotted-rougheye-SST_files/figure-gfm/map-indiv-data-1.png)<!-- -->

``` r
ggsave("pdf_outputs/blackspotted_rougheye_OwensLab_map.pdf", width = 6, height = 5)
```

    ## Warning: Removed 98 rows containing missing values (`geom_point()`).

``` r
blackspottedRougheye_sppID_results_20230818 <- results_for_map %>%
  select(-id, -mixture_collection, -log_likelihood, -z_score, -n_miss_loci, -n_non_miss_loci, -missing_loci, -PofZ) %>%
  rename(ABLG = Sample_ID)
  
write_csv(blackspottedRougheye_sppID_results_20230818, "csv_outputs/blackspottedRougheye_sppID_results_20230818.csv")
```
