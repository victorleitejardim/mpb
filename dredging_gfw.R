library(gfwr)
# Vessel search by MMSI
gfw_vessel_info(query = 224224000,
                search_type = "search")
# Encounter events for a date range
gfw_event(event_type = "ENCOUNTER",
          start_date = "2020-01-01",
          end_date = "2020-01-02")


#polygone bretagne
eu_sig <- st_read("data/FRA_adm_shp/FRA_adm0.shp")

brit2 <- st_crop(eu_sig, y = c(xmin = -4.6, xmax = -4.2, ymin = 48.249, ymax =48.4))

test <- gfw_ais_fishing_hours(spatial_resolution = "HIGH",
                      temporal_resolution = "HOURLY",
                      start_date = "2023-10-01",
                      end_date = "2024-03-01",
                      region_source = "USER_SHAPEFILE",
                      region = brit2
                      )
