Path.ACA.border = "Import/SIG/Extern_border/Out_border_ACA.shp"

if(file.exists(Path.ACA.border) == T){
  ACA.bo = sf::as_Spatial(sf::st_read(Path.ACA.border, quiet = TRUE))
  proj4string(ACA.bo) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  ACA.bo.proj = fortify(ACA.bo)

  Eurasia_map <- c("Mongolia", "Russia", "Spai n", "France", "Italy", "Greece", "Germany", "Finland", "Bhutan", "Bangladesh",
                   "Czech Republic", "Denmark", "Kosovo", "Vietnam", "Laos", "Japan", "Nepal", "India", "Myanmar",
                   "Latvia","Lithuania", "Estonia", "Belarus", "Romania", "Bulgaria", "Hungary", "Austria", "Croatia",
                   "Albania", "Serbia", "Slovenia", "Slovakia", "Bosnia", "Montenegro", "UK", "Ireland", "Moldova", "Macedonia",
                   "Norway", "Sweden", "Turkey", "Belgium", "Uzbekistan", "Tajikistan", "Syria", "Israel", "Jordan", "Pakistan",
                   "Kazakhstan","Turkmenistan", "Ukraine", "Poland", "Portugal", "Switzerland", "Kyrgyzstan", "Morocco",
                   "Kuwait", "Oman", "Saudi Arabia", "United Arab Emirates", "South Korea", "North Korea", "Qatar", "Bahrein",
                   "China", "Iran", "Armenia", "Georgia", "Afghanistan", "Iraq", "Azerbaijan")
  Eurasia_map <- map_data("world", region = Eurasia_map)
  World_map <- map_data("world")
   
  Map_theme_grid <- theme(
    axis.title = element_blank(), title = element_text(size = 20),
    legend.title = element_text(size = 12), legend.key = element_blank(), legend.justification = c("center"), legend.text = element_text(size = 10),
    legend.position = c(0.95, 0.35), legend.background = element_rect(fill = "#EBEBEB"),
    panel.background = element_rect(fill = "#EBEBEB", colour = "black"), panel.border = element_rect(fill = NA, colour = "black"),
    plot.background = element_blank(), plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm")
  )

  map_empty  <-  ggplot()+
    geom_polygon(data = Eurasia_map, aes(x=long, y=lat, group = group), alpha = 1, fill = "grey85", color = "grey30", size = 0.3)+
    coord_quickmap(xlim = c(-3,175), ylim = c(10, 80)) +
    guides(colour = guide_legend(override.aes = list(size = 5.5))) + Map_theme_grid

  map_empty_world  <-  ggplot()+
    geom_polygon(data = World_map, aes(x=long, y=lat, group = group), alpha = 1, fill = "grey85", color = "grey30", size = 0.3)+
    coord_quickmap(xlim = c(-3,175), ylim = c(10, 80)) +
    guides(colour = guide_legend(override.aes = list(size = 5.5))) + Map_theme_grid
}