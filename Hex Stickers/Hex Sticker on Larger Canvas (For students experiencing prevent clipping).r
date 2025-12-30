## If any of the packages are missing, run this code first:
## install.packages(c("hexSticker","ggplot2","ragg"))  

library(hexSticker)
library(ggplot2)
library(ragg)

## Set working directory or give full file path:
my_image = "plot_transparent.png"

sticker_plot = sticker(
  subplot  = my_image,
  package  = "UTM FSU \nBiostats",
  lineheight = 0.035,
  p_size = 100,
  p_x = 1, 
  p_y = 1.50,
  s_x = 1, 
  s_y = 0.84, 
  s_width = 0.65,
  p_color = "#FFFFFF",
  h_color = "#baf2e5", h_fill = "#275d95",
  h_size  = 1.50, # THICK BORDER
  url     = "https://nishanmudalige.github.io/UTM_FSU_Biostats/",
  u_color = "#FFFFFF", 
  u_size = 20,
  white_around_sticker = FALSE,
  filename = "R0_hex_sticker.png"
)


## Add padding around Hex sticker
extra_padding = 0.05    # adjust this value as required


## Add extra padding to sticker_plot
sticker_plot = sticker_plot +
  coord_fixed(clip = "off") +
  scale_x_continuous(expand = expansion(mult = 0, add = extra_padding)) +
  scale_y_continuous(expand = expansion(mult = 0, add = extra_padding)) +
  theme(plot.margin = margin(6, 6, 6, 6, "mm"))  


# Create a larger canvas
canvas_scale = 1.12 # Adjust this value as required
w_mm = 43.9 * canvas_scale # NOT recommended to adjust this value. Maintains aspect raio
h_mm = 50.8 * canvas_scale # NOT recommended to adjust this value. Maintains aspect raio

## Overwrite over 
ggplot2::ggsave(
  "R0_hex_sticker.png", 
  sticker_plot,
  bg = "transparent",
  width = w_mm, height = h_mm, units = "mm",
  dpi = 2500, limitsize = FALSE
)