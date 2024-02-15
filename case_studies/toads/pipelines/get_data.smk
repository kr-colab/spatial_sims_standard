rule all:
    input:
        aus_map = "data/AUS_2021_AUST_GDA2020.shp",
        gbif = "data/cane_toads_gbif.csv",
        bioclim = "data/wc2.1_10m_bio_12.tif",
        biomap = "data/australia_bio_12.png"

rule get_aus:
    output:
        rules.all.input.aus_map
    shell:
        """
        wget -O temp https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs-edition-3/jul2021-jun2026/access-and-downloads/digital-boundary-files/AUS_2021_AUST_SHP_GDA2020.zip
        unzip temp
        mv AUS_2021_AUST_GDA2020* data/
        rm temp
        """

rule get_gbif:
    output:
        rules.all.input.gbif
    shell:
        """
        wget -O temp https://api.gbif.org/v1/occurrence/download/request/0229932-220831081235567.zip
        unzip temp
        mv 0229932-220831081235567.csv {output}
        rm temp
        """

rule get_bioclim:
    output:
        rules.all.input.bioclim
    shell:
        """
        mkdir -p data/
        wget -O wc2.1_10m_bio.zip https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_bio.zip
        unzip wc2.1_10m_bio.zip
        mv wc2.1_10m_bio_*.tif data/
        rm wc2.1_10m_bio.zip
        """

rule biomap:
    input:
        aus_map = rules.all.input.aus_map,
        bioclim = "data/wc2.1_10m_bio_{i}.tif"
    output:
        "data/australia_bio_{i}.png"
    shell:
        """
        python src/make_biomap.py \
          --bioclim_tiff {input.bioclim} \
          --shapefile {input.aus_map} \
          --outfile {output}
        """
        
