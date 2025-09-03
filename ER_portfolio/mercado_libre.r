library(httr)
library(jsonlite)

# Function to fetch all pages of listings by category with optional location filter
fetch_all_listings <- function(category_id, location_id = NULL, page_size = 50) {
  base_url <- "https://api.mercadolibre.com/sites/MLU/search"
  all_listings <- list()
  
  offset <- 0
  while(TRUE) {
    query_params <- list(category = category_id, limit = page_size, offset = offset)

    # Add location filter if provided
    if (!is.null(location_id)) {
      query_params$location <- location_id
    }

    response <- GET(url = base_url, query = query_params)
    
    # Error handling
    if (http_error(response)) {
      stop("API request failed: ", http_status(response)$message)
    }

    data <- content(response, as = "parsed", simplifyDataFrame = TRUE)
    
    # Check if no results are returned or no more results available
    if (is.null(data$results) || length(data$results) == 0) {
      break
    }
    
    # Append the result page to the listings list
    all_listings <- append(all_listings, data$results)
    
    # Increase offset to fetch the next page
    offset <- offset + page_size
    
    # Optional: Pause after each request to avoid being rate-limited
    Sys.sleep(0.5)
  }
  
  # Convert list of results to a data frame
  if (length(all_listings) > 0) {
    listings_df <- do.call(rbind.data.frame, lapply(all_listings, function(x) data.frame(t(unlist(x)))))

    # Select relevant columns
    desired_columns <- c("id", "title", "price", "permalink", "condition", "available_quantity", "sold_quantity")
    listings_cleaned <- listings_df[, desired_columns, drop = FALSE]

    return(listings_cleaned)
  } else {
    message("No results found for category ", category_id)
    return(NULL)
  }
}

# Fetch categories and define necessary IDs
url_categories <- "https://api.mercadolibre.com/sites/MLU/categories"
categories_response <- GET(url_categories)
categories <- content(categories_response, as = "parsed", simplifyDataFrame = TRUE)

print("Available categories:")
print(categories)

# Define category IDs
real_estate_category <- "MLU1459"  # Replace with the correct category ID
vehicles_category <- "MLU1744"  # Replace with the correct category ID
location_montevideo <- "TUxVQ09VTmJhYjU0"  # Montevideo location ID

# Get and save all real estate listings (country-wide)
real_estate_listings <- fetch_all_listings(real_estate_category)

# Get and save all vehicle listings (country-wide)
vehicle_listings <- fetch_all_listings(vehicles_category)

# Get and save all real estate listings specifically in Montevideo
montevideo_real_estate <- fetch_all_listings(real_estate_category, location_montevideo)

# Save results to CSV (only if listings exist)
if (!is.null(real_estate_listings)) {
  write.csv(real_estate_listings, "real_estate_listings.csv", row.names = FALSE)
}
if (!is.null(vehicle_listings)) {
  write.csv(vehicle_listings, "vehicle_listings.csv", row.names = FALSE)
}
if (!is.null(montevideo_real_estate)) {
  write.csv(montevideo_real_estate, "montevideo_real_estate.csv", row.names = FALSE)
}

print("Data fetching completed. CSV files saved.")