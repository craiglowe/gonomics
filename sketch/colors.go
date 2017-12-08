package draw

import (
	"image/color"
)

/*
K. Kelly (1965): Twenty-two colors of maximum contrast. Color Eng., 3(6)
*/

var KellyPalette = color.Palette{
	color.RGBA{0xF2, 0xF3, 0xF4, 0xff},
	color.RGBA{0x22, 0x22, 0x22, 0xff},
	color.RGBA{0xF3, 0xC3, 0x00, 0xff},
	color.RGBA{0x87, 0x56, 0x92, 0xff},
	color.RGBA{0xF3, 0x84, 0x00, 0xff},
	color.RGBA{0xA1, 0xCA, 0xF1, 0xff},
	color.RGBA{0xBE, 0x00, 0x32, 0xff},
	color.RGBA{0xC2, 0xB2, 0x80, 0xff},
	color.RGBA{0x84, 0x84, 0x82, 0xff},
	color.RGBA{0x00, 0x88, 0x56, 0xff},
	color.RGBA{0xE6, 0x8F, 0xAC, 0xff},
	color.RGBA{0x00, 0x67, 0xA5, 0xff},
	color.RGBA{0xF9, 0x93, 0x79, 0xff},
	color.RGBA{0x60, 0x4E, 0x97, 0xff},
	color.RGBA{0xF6, 0xA6, 0x00, 0xff},
	color.RGBA{0xB3, 0x44, 0x6C, 0xff},
	color.RGBA{0xDC, 0xD3, 0x00, 0xff},
	color.RGBA{0x88, 0x2D, 0x17, 0xff},
	color.RGBA{0x8D, 0xB6, 0x00, 0xff},
	color.RGBA{0x65, 0x45, 0x22, 0xff},
	color.RGBA{0xE2, 0x58, 0x22, 0xff},
	color.RGBA{0x2B, 0x3D, 0x26, 0xff},
}

/*
https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
*/

var TrubetskoyPalette = color.Palette{
	color.RGBA{0xe6, 0x19, 0x4b, 0xff}, /* 0  red */
	color.RGBA{0x3c, 0xb4, 0x4b, 0xff}, /* 1  green */
	color.RGBA{0xff, 0xe1, 0x19, 0xff}, /* 2  yellow */
	color.RGBA{0x00, 0x82, 0xc8, 0xff}, /* 3  blue */
	color.RGBA{0xf5, 0x82, 0x31, 0xff}, /* 4  orange */
	color.RGBA{0x91, 0x1e, 0xb4, 0xff}, /* 5  purple */
	color.RGBA{0x46, 0xf0, 0xf0, 0xff}, /* 6  cyan */
	color.RGBA{0xf0, 0x32, 0xe6, 0xff}, /* 7  magenta */
	color.RGBA{0xd2, 0xf5, 0x3c, 0xff}, /* 8  lime */
	color.RGBA{0xfa, 0xbe, 0xbe, 0xff}, /* 9  pink */
	color.RGBA{0x00, 0x80, 0x80, 0xff}, /* 10 teal */
	color.RGBA{0xe6, 0xbe, 0xbe, 0xff}, /* 11 lavender */
	color.RGBA{0xaa, 0x6e, 0x28, 0xff}, /* 12 brown */
	color.RGBA{0xff, 0xfa, 0xc8, 0xff}, /* 13 beige */
	color.RGBA{0x80, 0x00, 0x00, 0xff}, /* 14 maroon */
	color.RGBA{0xaa, 0xff, 0xcf, 0xff}, /* 15 mint */
	color.RGBA{0x80, 0x80, 0x00, 0xff}, /* 16 olive */
	color.RGBA{0xff, 0xd8, 0xb1, 0xff}, /* 17 coral */
	color.RGBA{0x00, 0x00, 0x80, 0xff}, /* 18 navy */
	color.RGBA{0x80, 0x80, 0x80, 0xff}, /* 19 grey */
	color.RGBA{0xff, 0xff, 0xff, 0xff}, /* 20 white */
	color.RGBA{0x00, 0x00, 0x00, 0xff}, /* 21 black */
}