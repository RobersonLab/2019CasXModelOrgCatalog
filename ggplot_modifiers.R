gg_bigger_texts = theme(
    axis.title = element_text( size=22 ),
    axis.text = element_text( size=20 ),
    legend.text = element_text( size=14 ),
    legend.title = element_text( size=15 ),
    plot.title = element_text( size=22 ),
	#strip.text.x = element_text( size=15, face = "bold", margin = margin( t=2, b=2, l=2, r=2 ),
	strip.text.x = element_text( size=17, margin=margin( b=5, t=5 ) ),
	strip.text.y = element_text( size=15 )
)

gg_no_legend = theme(
	legend.position='none'
)

gg_no_grid = theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

gg_no_x_grid = theme( 
  panel.grid.major.x = element_blank() )

gg_no_y_grid = theme( 
  panel.grid.major.y = element_blank() )

gg_center_title = theme(
	plot.title = element_text( hjust = 0.5 )
)

gg_no_x_label = theme( 
	axis.title.x = element_blank()
)

gg_no_y_label = theme(
	axis.title.y = element_blank()
)

gg_angled_x_text = theme (
	axis.text.x = element_text( angle = 45, vjust = 1, hjust = 1, color = 'black' )
)
###############################
# color blind friendly colors #
###############################
colorBlindPalette <- c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )
blues5Palette <- c( '#ffffcc', '#a1dab4', '#41b6c4', '#2c7fb8', '#253494' )
greens5Palette <- c( '#ffffcc', '#c2e699', '#78c679', '#31a354', '#006837' )
purples5Palette <- c( '#feebe2', '#fbb4b9', '#f768a1', '#c51b8a', '#7a0177' )
reds5Palette <- c( '#ffffb2', '#f3cc5c', '#fd8d3c', '#f03b20', '#bd0026' )
reds4Palette <- c( '#fef0d9', '#fdcc8a', '#fc8d59', '#d7301f' )
reds3Palette <- c( '#fee0d2', '#fc9272', '#de2d26' )
