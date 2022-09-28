

extract_raster_prop <- function(
    my_points,
    location_column,
    my_buffer,
    my_raster_data,
    lulc_cats = NULL
){
  
  sites <- my_points[,location_column]
  
  cli::cli_h1("Reprojecting my_points to map projection")
  points_RP <- sf::st_transform(
    my_points,
    sf::st_crs(my_raster_data)
  )
  
  cli::cli_alert_success("my_points reprojected")
  
  # Step 2.
  # use the raster::extract function to extract the mean value of the raster data
  # within a particular buffer around each site.
  # We need a sub-function to calculate the proportion of each category
  spatial_summary <- function(
    x,
    ncats = my_raster_data@data@max,
    ...
  ){
    return(
      prop.table(
        tabulate(x, ncats)
      )
    )
  }
  
  cli::cli_h1("Extracting spatial data")
  prop_extract <- suppressWarnings(
    raster::extract(
      my_raster_data,
      points_RP,
      fun=spatial_summary,
      buffer= my_buffer,
      na.rm = TRUE
    )
  )
  
  
  if(is.numeric(prop_extract)){
    prop_extract <- matrix(
      prop_extract,
      ncol = my_raster_data@data@max,
      nrow = nrow(points_RP),
      byrow = TRUE
    )
  }
  
  
  # if lulc_cats is a list
  if(is.list(lulc_cats)){
    prop_extract <- apply(
      prop_extract,
      1,
      function(x){
        sapply(
          lulc_cats,
          function(y)
            ifelse(
              length(y) == 1,
              x[y],
              sum(x[y])
            )
        )
      }
    )
    if(length(lulc_cats) == 1){
      prop_extract <- t(t(prop_extract))
    }else{
      prop_extract <- t(prop_extract)
    }
  }
  
  # if it is a numeric
  if(is.numeric(lulc_cats)){
    prop_extract <- prop_extract[,lulc_cats]
  }
  # if it's a names list
  if(!is.null(names(lulc_cats)) & is.list(lulc_cats)){
    colnames(prop_extract) <- names(lulc_cats)
  }
  
  # create dataframe matching the sites with the extracted data
  df <- data.frame(
    LocationName = sites,
    prop_extract,
    row.names = NULL
  )
  # give the site column the same name as my_points
  if(is.character(location_column)){
    colnames(df)[1] <- location_column
  } else {
    colnames(df)[1] <- colnames(my_points)[location_column]
  }
  cli::cli_alert_success("Spatial data extracted")
  return(df[order(df[,1]),])
}