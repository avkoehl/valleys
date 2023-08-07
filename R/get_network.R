# gets NDH HR data from the web
# splits the lines into 1km segments
# and saves it as dataframe to the data folder
library('nhdplusTools')

dl_network = function() {
  # download the network
  # and save it to the data folder
  # as a dataframe
  # returns the dataframe
  # requires nhdplusTools
  #
}

split_network = function() {
  # splits the network into 1km segments
  # and saves it to the data folder
  # as a dataframe
  # returns the dataframe
  #
}

get_attributes = function() {
  # gets the attributes for each segment
  # and saves it to the data folder
  # as a dataframe
  # returns the dataframe
  #

    # attributes are for stratified sample
    # confinement
    # slope
    # stream order
    # stability ?
}

get_network = function() {
    # dl_network()
    # split_network()
    # get_attributes()
    # return the dataframe
}

# for all of california
hr_urls = download_hr_data("../data/ndhdplus_data/"
