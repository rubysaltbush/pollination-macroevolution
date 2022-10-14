# reads in data on pollination syndromes, makes names consistent with RB2020
# tree, separates pollination data into columns categorising it differently for
# different ASR models (i.e. animal+abiotic, wind, water + animal)

pollination1209 <- cache_RDS("data_output/pollination1209.csv", read_function = readr::read_csv,
  save_function = write_csv, function() {
    # read in original data
    pollination1209 <- readr::read_csv("data_input/qryeFLOWER-1209taxa_AusTraits_TRY_2021_data_entry.csv")

    # remove austraits and TRY columns from data, condense to columns of
    # interest for now
    pollination1209 <- pollination1209 %>%
      dplyr::select(taxon_name:ParCladeTax, pollination:wind_pollination_explicitly_tested,
        pollination_data_source:scent, rank_scored_at) %>%
      dplyr::distinct()
    # replace spaces in species names with underscores and get rid of hyphens
    # to match tree
    pollination1209$taxon_name <- gsub(" ", "_", pollination1209$taxon_name)
    pollination1209$taxon_name <- gsub("-", "", pollination1209$taxon_name)

    # adding tip numbers (position in tree) to pollination1209 df
    treetips <- data.frame(taxon_name = tree$tip.label, position = c(1:1201))
    pollination1209 <- pollination1209 %>%
      dplyr::left_join(treetips, by = "taxon_name")
    rm(treetips)

    # reduce confidence score with definition written out to just the score
    pollination1209 <- pollination1209 %>%
      dplyr::mutate(conf_score = gsub(" .*", "", pollination1209$confidence_score)) %>%
      dplyr::select(-confidence_score)
    table(pollination1209$conf_score)

    # new columns separating pollination values out in different ways
    
    # abiotic and animal
    pollination1209 <- pollination1209 %>%
      dplyr::mutate(abiotic_animal = gsub("wind|water", "abiotic", pollination1209$wind_water_animal_pollination))
    pollination1209$abiotic_animal <- gsub("abiotic&abiotic", "abiotic", pollination1209$abiotic_animal)
    # replace NAs with ? for corHMM functions
    pollination1209$abiotic_animal <- tidyr::replace_na(pollination1209$abiotic_animal, "?")
    table(pollination1209$abiotic_animal)
    
    # then abiotic and animal with polymorphisms as a third state
    pollination1209$abiotic_animal_nopoly <- gsub("abiotic&animal", "abiotic_animal", pollination1209$abiotic_animal)
    table(pollination1209$abiotic_animal_nopoly)

    # wind, water, vertebrate, insect
    # 4 states with polymorphisms
    pollination1209$wind_water_vert_insect <- pollination1209$pollination
    pollination1209$wind_water_vert_insect <- gsub(" and ", "&", pollination1209$wind_water_vert_insect)
    pollination1209$wind_water_vert_insect <- gsub("bat|hummingbird|bird|lizard|mammal|rodent|small mammal",
      "vertebrate", pollination1209$wind_water_vert_insect)
    pollination1209$wind_water_vert_insect <- gsub("bee|moth|beetle|butterfly|fly|thrip|gall midge|hymenoptera|lepidoptera|wasp|cricket&cockroach",
      "insect", pollination1209$wind_water_vert_insect)
    pollination1209$wind_water_vert_insect <- gsub("abiotic", "wind&water", pollination1209$wind_water_vert_insect)
    pollination1209$wind_water_vert_insect <- gsub("animal", "vertebrate&insect",
      pollination1209$wind_water_vert_insect)
    pollination1209$wind_water_vert_insect <- gsub("insect&insect",
      "insect", pollination1209$wind_water_vert_insect)
    pollination1209$wind_water_vert_insect <- gsub("insect&vertebrate|vertebrate&insect&vertebrate|vertebrate, vertebrate&insect",
      "vertebrate&insect", pollination1209$wind_water_vert_insect)
    pollination1209$wind_water_vert_insect <- gsub("insect&wind", "wind&insect",
      pollination1209$wind_water_vert_insect)
    # replace NAs with ? for corHMM functions
    pollination1209$wind_water_vert_insect <- tidyr::replace_na(pollination1209$wind_water_vert_insect, "?")
    table(pollination1209$wind_water_vert_insect)
    # NOT coding this with polymorphisms as new states as this results in 9 states
    # which is way too many

    # wind, water & animal, 3 states with polymorphisms
    # replace NAs with ? for corHMM functions
    pollination1209$wind_water_animal <- tidyr::replace_na(pollination1209$wind_water_animal_pollination, "?")
    table(pollination1209$wind_water_animal)

    # code wind, water, animal polymorphic taxa as new states (6 states total)
    pollination1209$wind_water_animal_nopoly <- gsub("&", "_", pollination1209$wind_water_animal)
    table(pollination1209$wind_water_animal_nopoly)
    
    # abiotic, vertebrate & insect, 3 states with polymorphisms
    pollination1209$abiotic_vert_insect <- gsub("wind|water|wind&water", "abiotic", pollination1209$wind_water_vert_insect)
    table(pollination1209$abiotic_vert_insect)
    
    # code vertebrate, insect, abiotic polymorphic taxa as new states (6 states total)
    pollination1209$abiotic_vert_insect_nopoly <- gsub("&", "_", pollination1209$abiotic_vert_insect)
    table(pollination1209$abiotic_vert_insect_nopoly)
    
    # binary characters for stochastic character mapping
    # wind and animal (drop water)
    pollination1209$wind_animal <- gsub("water&animal|wind&water|water", "?", pollination1209$wind_water_animal)
    table(pollination1209$wind_animal)
    
    # vertebrate and insect (drop abiotic)
    pollination1209$vert_insect <- gsub("abiotic&insect|abiotic&vertebrate|abiotic", "?", pollination1209$abiotic_vert_insect)
    table(pollination1209$vert_insect)
    
    # for colouring polymorphic tips in trees (but not ASR as too many states)
    # add column with no poly for wind_water_vert_insect
    pollination1209$wind_water_vert_insect_nopoly <- gsub("&", "_", pollination1209$wind_water_vert_insect)
    table(pollination1209$wind_water_vert_insect_nopoly)
    
    # replace named character states with integers (necessary for plotting later)
    pollination1209 <- pollination1209 %>%
      dplyr::mutate(dplyr::across(19:28, ~gsub("water_insect", "14", .))) %>%
      dplyr::mutate(dplyr::across(19:28, ~gsub("wind_vertebrate", "13", .))) %>%
      dplyr::mutate(dplyr::across(19:28, ~gsub("wind_insect", "12", .))) %>%
      dplyr::mutate(dplyr::across(19:28, ~gsub("abiotic_vertebrate", "11", .))) %>%
      dplyr::mutate(dplyr::across(19:28, ~gsub("abiotic_insect", "10", .))) %>%
      dplyr::mutate(dplyr::across(19:28, ~gsub("water_animal", "9", .))) %>%
      dplyr::mutate(dplyr::across(19:28, ~gsub("wind_animal", "8", .))) %>%
      dplyr::mutate(dplyr::across(19:28, ~gsub("abiotic_animal", "7", .))) %>%
      dplyr::mutate(dplyr::across(19:28, ~gsub("insect", "5", .))) %>%
      dplyr::mutate(dplyr::across(19:28, ~gsub("vertebrate", "6", .))) %>%
      dplyr::mutate(dplyr::across(19:28, ~gsub("abiotic", "1", .))) %>%
      dplyr::mutate(dplyr::across(19:28, ~gsub("animal", "4", .))) %>%
      dplyr::mutate(dplyr::across(19:28, ~gsub("water", "3", .))) %>%
      dplyr::mutate(dplyr::across(19:28, ~gsub("wind", "2", .))) %>%
      #leaves some "no poly" states as e.g. 2_3 wind_water, code as 1 abiotic
      dplyr::mutate(dplyr::across(19:28, ~gsub("2_3", "1", .))) %>%
      #leaves some "no poly" states as e.g. 6_5 vert_insect, code as 4 animal
      dplyr::mutate(dplyr::across(19:28, ~gsub("6_5", "4", .)))
    
    #### TABLES ####
    # syndromes vs systems?
    table(pollination1209$syndrome_or_system)
    # 734 syndromes, 425 systems, makes sense looking at confidence scores
    
    # how many families covered by final data?
    table(pollination1209$wind_water_vert_insect)
    # number of families in final data? 434?!
    nrow(pollination1209 %>% dplyr::select(ParFamTax) %>% distinct())
    # number of families WITH data 433!
    nrow(pollination1209 %>% dplyr::filter(wind_water_vert_insect != "?") %>% dplyr::select(ParFamTax) %>% distinct())
    # number of families withOUT data 35!
    nrow(pollination1209 %>% dplyr::filter(wind_water_vert_insect == "?") %>% dplyr::select(ParFamTax) %>% distinct())
    # families missing from final data
    without <- pollination1209 %>% dplyr::filter(wind_water_vert_insect == "?") %>% dplyr::select(ParFamTax) %>% distinct()
    with <- pollination1209 %>% dplyr::filter(wind_water_vert_insect != "?") %>% dplyr::select(ParFamTax) %>% distinct()
    without[!(without$ParFamTax %in% with$ParFamTax),]
    # Hoplestigmataceae
    rm(with, without)
    
    # how often has wind pollination been explicitly tested?
    table(pollination1209$wind_pollination_explicitly_tested)
    # only 47 out of 1201! not often explicitly tested, even (or especially?)
    # in taxa assumed to be wind pollinated
    table(pollination1209[pollination1209$abiotic_animal == "1",]$wind_pollination_explicitly_tested)
    # only explicitly tested in 10 out of 130 abiotically pollinated taxa

    # most common data source?
    table(pollination1209$pollination_data_source)
    # most explicit from text, then guessed from illustration, then interpreted
    # from text

    # ranks scored at? ONLY for those WITH data!!
    table(pollination1209[pollination1209$abiotic_animal != "?",]$rank_scored_at)

    # how many taxa have information about pollination rewards? what are they
    # all?
    table(pollination1209$pollination_reward)
    # nectar is the most commonly described reward, followed by nectar and
    # pollen, then pollen

    # how many taxa have nectar? 249 with, 70 defs without
    table(pollination1209$nectar)
    # how many have scent? 141 with, 26 defs without
    table(pollination1209$scent)
    
    #### BIOME DATA ####
    # add biome data from RB et al 2020
    # read in biomes file - shows biome by cleaned GBIF records of contemporary
    # species distribution, with different thresholds for inclusion in a particular
    # biome. 7 species with arctic/antarctic distributions are excluded. e.g. SB50 = >50%
    # of species records are in a particular superbiome (temperate = 0, tropical = 1
    # arid = 2)
    pollination1209_biome <- readr::read_csv("data_input/ATT2_superbiomes_v1.2_NoCold.csv")
    # some names in biome data have been updated in my pollination data
    # but are not updated in (older) biome data
    pollination1209_biome$taxon_name <- gsub("Anagallis_tenella", "Lysimachia_tenella", pollination1209_biome$taxon_name)
    pollination1209_biome$taxon_name <- gsub("Sedum_rubrotinctum", "Sedum_x_rubrotinctum", pollination1209_biome$taxon_name)
    pollination1209_biome$taxon_name <- gsub("Tetracarpaea_tasmanica", "Tetracarpaea_tasmannica", pollination1209_biome$taxon_name)
    
    # join onto pollination 1209
    pollination1209 <- pollination1209 %>%
      dplyr::left_join(pollination1209_biome, by = "taxon_name")
    rm(pollination1209_biome)

    #### GBIF location data ####
    # read in GBIF location data downloaded and cleaned by Will Cornwell
    gbif <- read_csv("data_input/filtered_obs_ruby.csv")
    
    # filter to preserved specimens only, human obs records a bit too iffy
    gbif_specimens <- gbif %>%
      dplyr::filter(basisofrecord == "PRESERVED_SPECIMEN")
    # still have matches for most species
    
    # add _ to species name (and rename taxon_name) in gbif to match pollination1209
    gbif_specimens$taxon_name <- gsub(" ", "_", gbif_specimens$species)
    # remove - from taxon_name in gbif to match pollination1209
    gbif_specimens$taxon_name <- gsub("-", "", gbif_specimens$taxon_name)
    # taxonomic updates and mismatches between gbif and pollination data
    # ONLY updating names for taxa that are recognised as synonyms in POWO
    # some GBIF names don't match and this is legitimate based on e.g. their distributions
    # in some cases our taxonomy is out of date, in some cases GBIF taxonomy is out of date
    gbif_specimens$taxon_name <- gsub("Tetracarpaea_tasmanica", "Tetracarpaea_tasmannica", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Citrus_aurantium", "Citrus_x_aurantium", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Tarenaya_rosea", "Cleome_rosea", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Hypopitys_monotropa", "Monotropa_hypopitys", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Glischrocaryon_racemosum", "Haloragodendron_racemosum", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Nama_demissum", "Nama_demissa", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Bdallophyton_americanum", "Bdallophytum_americanum", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Grahamia_frutescens", "Talinopsis_frutescens", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Sophronanthe_pilosa", "Gratiola_pilosa", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Morella_cerifera", "Myrica_cerifera", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Aphyllon_fasciculatum", "Orobanche_fasciculata", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Staphylea_japonica", "Euscaphis_japonica", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Isotrema_durius", "Aristolochia_macrophylla", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Hilairanthus_bicolor", "Avicennia_bicolor", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Octomeles_sumatranum", "Octomeles_sumatrana", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Diplacus_aurantiacus", "Mimulus_aurantiacus", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Sanicula_odorata", "Sanicula_gregaria", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Craterostigma_nummulariifolium", "Lindernia_nummulariifolia", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Miconia_petiolaris", "Clidemia_petiolaris", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Hololachna_soongarica", "Reaumuria_soongarica", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Liquidambar_excelsa", "Altingia_excelsa", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Rhodoleia_championiae", "Rhodoleia_championii", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Procris_repens", "Pellionia_repens", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Actinostemma_paniculatum", "Bolbostemma_paniculatum", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Manekia_incurva", "Manekia_sydowii", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Durandea_jenkinsii", "Durandea_pentagyna", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Muraltia_pauciflora", "Polygala_pauciflora", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Sedum_rubrotinctum", "Sedum_x_rubrotinctum", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Morus_alba", "Morus_indica", gbif_specimens$taxon_name) # can't believe only two records available for Morus indica??? Weird
    gbif_specimens$taxon_name <- gsub("Marathrum_foeniculaceum", "Marathrum_rubrum", gbif_specimens$taxon_name)
    gbif_specimens$taxon_name <- gsub("Fagus_grandiflora", "Fagus_grandifolia", gbif_specimens$taxon_name)
    # all remaining unmatched taxa from GBIF records can be discarded as they
    # are taxonomically suspect (i.e. distribution does not match species distribution)
    
    # check species coverage if gbif "species" matched to taxa in pollination1209
    sum(pollination1209$taxon_name %in% gbif_specimens$taxon_name)
    # 1140 of 1201 pollination1209 taxa with gbif herbarium records
    
    # clean GBIF records (already cleaned by Will Cornwell but left sea records in)
    # convert country code from ISO2c to ISO3c for cleaning function
    gbif_specimens$countrycode <- countrycode::countrycode(gbif_specimens$countrycode, origin = 'iso2c', destination = 'iso3c')
    gbif_specimens <- data.frame(gbif_specimens)
    # then test out flagging records for potential errors
    flags <- CoordinateCleaner::clean_coordinates(x = gbif_specimens,
                                                  countries = "countrycode",
                                                  tests = c("capitals", 
                                                            "centroids", 
                                                            "equal","gbif", 
                                                            "institutions",
                                                            "outliers", "seas",
                                                            "zeros", "countries"))
    # Flagged 32063 of 337450 records, EQ = 0.1.
    
    # remove all records flagged as not on land (this will include several water 
    # pollinated taxa but these will not be considered in analyses of transitions 
    # between wind and animal pollination so this is okay)
    gbif_spec_filtered <- flags %>%
      dplyr::filter(.sea == TRUE) %>%
    # and remove all outliers as many of these outside native range according to POWO
      dplyr::filter(.otl == TRUE)
    
    paste("There are", 
          length(gbif_spec_filtered$.summary) - sum(gbif_spec_filtered$.summary),
          "records remaining with a flagged issue, but these have been checked 
          and are okay")
    
    # check species coverage after cleaning
    sum(pollination1209$taxon_name %in% gbif_spec_filtered$taxon_name)
    # 1137 of 1201 pollination1209 taxa with gbif herbarium records, lost 3 species
    
    # work out number of GBIF records per taxon
    no_records <- gbif_spec_filtered %>%
      dplyr::group_by(taxon_name) %>%
      dplyr::summarise(no_records = length(decimallatitude))
    pollination1209 <- pollination1209 %>%
      dplyr::left_join(no_records, by = "taxon_name")
    rm(no_records)
    paste("There are", mean(pollination1209$no_records, na.rm = TRUE), "GBIF records per taxon on average")
    # 271.6 records per taxon
    # calculate standard error for number of records per taxon
    sqrt(var(pollination1209$no_records, na.rm = TRUE) / length(pollination1209$no_records))
    # 12.9 standard error - minimum 1 (for 26 taxa), maximum 5568
    # if we had a cutoff of e.g. 10 how many would we lose?
    test_cutoff <- pollination1209[!(is.na(pollination1209$no_records)),]
    sum(test_cutoff$no_records < 10)
    # 162 with records fewer than ten
    sum(test_cutoff$no_records < 5)
    # 100 with records fewer than 5 BUT these are likely genuine records of
    # rarely encountered/obscure taxa, will leave for now
    rm(test_cutoff)
    
    rm(gbif, gbif_specimens, flags)
    
    # calculate average latitude for each species
    latitude <- gbif_spec_filtered %>%
      dplyr::group_by(taxon_name) %>%
      dplyr::summarise(meanlat = mean(decimallatitude))
    # and average ABSOLUTE latitude
    abs_latitude <- gbif_spec_filtered %>%
      dplyr::group_by(taxon_name) %>%
      dplyr::summarise(abs_meanlat = mean(abs(decimallatitude)))
    # join onto pollination1209 data
    pollination1209 <- pollination1209 %>%
      dplyr::left_join(latitude, by = "taxon_name") %>%
      dplyr::left_join(abs_latitude, by = "taxon_name")
    rm(latitude, abs_latitude)
    
    # time to get spatial!
    #first duplicate latitude/longitude columns so I don't lose these
    gbif_filt_spatial <- gbif_spec_filtered %>%
      dplyr::mutate(long = decimallongitude, lat = decimallatitude)
    #convert gbif_spec_filtered to an sf spatial object and assign CRS of WGS84
    gbif_filt_spatial <- sf::st_as_sf(gbif_filt_spatial, coords = c("long","lat"))
    sf::st_crs(gbif_filt_spatial) <- 4326
    
    # try plotting e.g. taxon with most records
    gbif_filt_spatial %>%
      dplyr::filter(taxon_name == "Casearia_sylvestris") %>%
      dplyr::select(countrycode, geometry) %>%
      plot()
    # looks like distribution on POWO! yay
    
    # read in rasters of LAI, mean temperature, mean rainfall
    # use geodata package to access WorldClim data for mean temp and precipitation
    # takes a while to download first time - reset R options so it doesn't timeout
    options(timeout = 10000)
    worldclim <- geodata::worldclim_global(var = "bio", res = 2.5,
                                           path = "cache/")
    
    # BIO1 = Annual Mean Temperature
    meantemp <- worldclim$wc2.1_2.5m_bio_1
    plot(meantemp)
    # BIO12 = Annual Precipitation
    totprec <- worldclim$wc2.1_2.5m_bio_12
    plot(totprec)
    
    rm(worldclim)
    
    # might not look at wind speed for now as not sure how reliable a coarse
    # measure of wind speed will be at this scale, given windiness probably 
    # varies a lot with topography, humidity and habitat openness
    
    # read in global raster of mean LAI for 1981-2020, calculated by taking
    # mean of GLOBMAP V3 (https://doi.org/10.5281/zenodo.4700264) ~8km res
    # averaging method documented in prep_LAI.R
    globalmeanLAI1981_2020 <- terra::rast("data_output/globalmeanLAI1981_2020.tif")
    plot(globalmeanLAI1981_2020)
    # slightly different extent to WorldClim rasters but otherwise lines up
    
    # intersect records with rasters, calculate average LAI, temp and rain for each
    record_meantemp <- terra::extract(meantemp, vect(gbif_filt_spatial))
    record_totprec <- terra::extract(totprec, vect(gbif_filt_spatial))
    record_meanLAI <- terra::extract(globalmeanLAI1981_2020, vect(gbif_filt_spatial))
    rm(meantemp, totprec, globalmeanLAI1981_2020)
    
    # corrupts gbif_filt_spatial if not calculated separately then joined as below
    gbif_filt_spatial <- gbif_filt_spatial %>%
      dplyr::mutate(ID = as.numeric(rownames(gbif_filt_spatial))) %>%
      dplyr::left_join(record_meantemp, by = "ID") %>%
      dplyr::left_join(record_totprec, by = "ID") %>%
      dplyr::left_join(record_meanLAI, by = "ID") %>%
      dplyr::select(-ID)
    rm(record_meantemp, record_totprec, record_meanLAI)
    
    # remove sf spatial element from gbif_filt_spatial otherwise VERY slow
    gbif_spec_filtered <- sf::st_drop_geometry(gbif_filt_spatial)
    # and output this data just in case
    write_csv(gbif_spec_filtered, "data_output/gbif_specimens_filtered.csv")
    rm(gbif_filt_spatial)
    
    # calculate mean MATemperature per species
    species_meantemp <- gbif_spec_filtered %>%
      dplyr::select(taxon_name, wc2.1_2.5m_bio_1) %>%
      dplyr::group_by(taxon_name) %>%
      dplyr::summarise(meanMAT = mean(wc2.1_2.5m_bio_1, na.rm = TRUE))
    
    # calculate mean TAPrecipitation per species
    species_totprec <- gbif_spec_filtered %>%
      dplyr::select(taxon_name, wc2.1_2.5m_bio_12) %>%
      dplyr::group_by(taxon_name) %>%
      dplyr::summarise(meanTAP = mean(wc2.1_2.5m_bio_12, na.rm = TRUE))
    
    # calculate mean LeafAI per species
    species_meanLAI <- gbif_spec_filtered %>%
      dplyr::select(taxon_name, globalmeanLAI1981_2020) %>%
      dplyr::group_by(taxon_name) %>%
      dplyr::summarise(meanLAI = mean(globalmeanLAI1981_2020, na.rm = TRUE))
    # 2 species with 0 mean LAI occur on small islands, LAI measure likely in 
    # error for these as may be averaging with ocean LAI (which 0)
    # remove 2 taxa with 0 LAI
    species_meanLAI <- species_meanLAI %>%
      dplyr::filter(meanLAI > 0)
    
    # per Rachael's request, also calculate species median and SEM LAI
    # Standard Error of the Mean
    species_SEMLAI <- gbif_spec_filtered %>%
      dplyr::select(taxon_name, globalmeanLAI1981_2020) %>%
      dplyr::group_by(taxon_name) %>%
      dplyr::summarise(SEMLAI = sd(globalmeanLAI1981_2020, na.rm = TRUE)/sqrt(length(globalmeanLAI1981_2020)))
    
    # species median LAI
    species_medianLAI <- gbif_spec_filtered %>%
      dplyr::select(taxon_name, globalmeanLAI1981_2020) %>%
      dplyr::group_by(taxon_name) %>%
      dplyr::summarise(medianLAI = median(globalmeanLAI1981_2020, na.rm = TRUE))
    # again taxa with 0 median LAI occur on small islands, likely 0 in error
    # remove 5 taxa with 0 LAI
    species_medianLAI <- species_medianLAI %>%
      dplyr::filter(medianLAI > 0)
    
    # join species means to pollination data
    pollination1209 <- pollination1209 %>%
      dplyr::left_join(species_meantemp, by = "taxon_name") %>%
      dplyr::left_join(species_totprec, by = "taxon_name") %>%
      dplyr::left_join(species_meanLAI, by = "taxon_name") %>%
      dplyr::left_join(species_medianLAI, by = "taxon_name") %>%
      dplyr::left_join(species_SEMLAI, by = "taxon_name")
    rm(species_meantemp, species_totprec, species_meanLAI, species_medianLAI,
       species_SEMLAI, gbif_spec_filtered)
    
    # add category of tropical/nontropical based on species' mean absolute latitude
    pollination1209 <- pollination1209 %>%
      dplyr::mutate(trop_nontrop = if_else(abs_meanlat <= 23.43635, "tropical", "nontropical"))
    
    # work out appropriate cutoff for LAI - from mapping about 40 looks appropriate??
    sum(pollination1209$meanLAI >= 40, na.rm = TRUE) # 280 above 40, maybe too high?
    sum(pollination1209$meanLAI >= 30, na.rm = TRUE) # 437 above 30
    
    # try kmeans clustering to determine best cutoff point for LAI
    # subset just taxon and LAI
    forcluster <- pollination1209 %>%
      dplyr::select(meanLAI, taxon_name) %>%
      dplyr::filter(!is.na(meanLAI)) %>%
      as.data.frame()
    rownames(forcluster) <- forcluster$taxon_name
    forcluster[,2] <- NULL
    
    # cluster using kmeans to put into two categories
    kmeansLAI <- kmeans(forcluster, centers = 2)
    open_closed <- data.frame(forcluster, open_closed = kmeansLAI$cluster)
    open_closed$taxon_name <- rownames(open_closed)
    rownames(open_closed) <- NULL
    # kmeans clustering makes cutoff at 30 LAI
    open_closed <- open_closed %>%
      dplyr::select(open_closed, taxon_name)
    
    # join clusters back onto pollination data
    # 1 is closed vegetation, 2 is open vegetation
    pollination1209 <- pollination1209 %>%
      dplyr::left_join(open_closed, by = "taxon_name")
    rm(open_closed, forcluster, kmeansLAI)
    table(pollination1209$open_closed)
    
    # scatter plot of LAI vs absolute latitude coloured by superbiome
    # out of interest
    pdf(file="figures/meanLAI_meanlat_superbiome.pdf", width = 8, height = 5)
    par(las = 1, bty = "l") # remove plot outline
    palette(c("black", "#91bfdb", "#fc8d59", "#d7191c"))# set colour palette
    plot(meanLAI ~ abs_meanlat, 
         data = pollination1209, 
         pch = 18, cex = 1, 
         col = factor(SB50), #(temperate = 0, tropical = 1, arid = 2)
         xlab = "Species mean latitude (absolute)", 
         ylab = "Species mean Leaf Area Index")
    legend("topright",
           legend = c("temperate", "tropical", "arid", "no biome data"),
           col = c("#91bfdb", "#fc8d59", "#d7191c", "black"),
           pch = 18, bty = "n")
    dev.off()
    palette("default") # set colour palette back to default
    
    #FAMILY POLLINATION BREAKDOWN
    # how many families in data overall?
    families <- pollination1209 %>% 
      dplyr::filter(wind_water_vert_insect != "?") %>%
      dplyr::select(ParFamTax) %>%
      dplyr::distinct()
    # "There are 433 families with pollination data"
    
    # how many families contain wind pollinated taxa?
    windpollfams <- pollination1209 %>% 
      dplyr::filter(wind_water_animal %in% c("2", "2&3", "2&4")) %>%
      dplyr::select(ParFamTax) %>% 
      dplyr::distinct()
    # 100, a fair few! 100/433 = 23 % of families
    
    # how many families contain animal pollinated taxa?
    animalpollfams <- pollination1209 %>% 
      dplyr::filter(wind_water_animal %in% c("4", "2&4", "3&4")) %>%
      dplyr::select(ParFamTax) %>% 
      dplyr::distinct()
    # 389, majority by far - 389/433 = 89.8% of families
    # 373 if excluding ambophily - 373/433 = 86% of families
    
    # how many families water pollinated?
    waterpollfams <- pollination1209 %>% 
      dplyr::filter(wind_water_animal %in% c("3", "2&3", "3&4")) %>%
      dplyr::select(ParFamTax) %>% 
      dplyr::distinct()
    # 8 families with water pollination! 8/433 = 1.9%
    # 6 families if excluding polymorphies, 6/433 = 1.4%
    
    # how many families ONLY wind or ONLY animal pollinated?
    # first work out how many families overlap
    windanimalfams <- animalpollfams %>%
      dplyr::filter(ParFamTax %in% windpollfams$ParFamTax|
                      ParFamTax %in% waterpollfams$ParFamTax)
    # and also add polymorphic water pollinated families
    waterpolyfams <- pollination1209 %>% 
      dplyr::filter(wind_water_animal %in% c("2&3", "3&4")) %>%
      dplyr::select(ParFamTax) %>% 
      dplyr::distinct()
    windwateranimalfams <- rbind(windanimalfams, waterpolyfams)
    windwateranimalfams <- dplyr::distinct(windwateranimalfams)
    # 64 families with both wind, water and animal pollination - 64/433 = 15%
    # now animal only families?
    animalonly <- animalpollfams %>%
      dplyr::filter(!(ParFamTax %in% windwateranimalfams$ParFamTax))
    # 327 families entirely animal pollinated - 76%
    # wind only families?
    windonly <- windpollfams %>%
      dplyr::filter(!(ParFamTax %in% windwateranimalfams$ParFamTax))
    # 37 families entirely wind pollinated - 37/433 = 9%
    # water only families?
    wateronly <- waterpollfams %>%
      dplyr::filter(!(ParFamTax %in% windwateranimalfams$ParFamTax))
    # 5 families entirely water pollinated - 5/433 = 1.2%
    rm(windpollfams, animalpollfams, animalonly, windonly, windanimalfams,
       wateronly, windwateranimalfams, waterpollfams, waterpolyfams, families)
    
    # write pollination1209 to data output folder so it can be cached!
    write_csv(pollination1209, "data_output/pollination1209.csv")
  })
