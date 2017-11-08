package draw

import (
	"image"
	"image/color"
)

func HLine(img *image.RGBA, xStart int, xEnd int, y int, col color.Color) {
	for x := xStart; x <= xEnd; x++ {
		img.Set(x, y, col)
	}
}

func VLine(img *image.RGBA, x int, yStart int, yEnd int, col color.Color) {
	for y := yStart; y <= yEnd; y++ {
		img.Set(x, y, col)
	}
}

func Rectangle(img *image.RGBA, xOne int, yOne int, xTwo int, yTwo int, col color.Color) {
	HLine(img, xOne, xTwo, yOne, col)
	HLine(img, xOne, xTwo, yTwo, col)
	VLine(img, xOne, yOne, yTwo, col)
	VLine(img, xTwo, yOne, yTwo, col)
}

func FilledRectangle(img *image.RGBA, xOne int, yOne int, xTwo int, yTwo int, col color.Color) {
	for x := xOne; x <= xTwo; x++ {
		for y := yOne; y <= yTwo; y++ {
			img.Set(x, y, col)
		}
	}
}
