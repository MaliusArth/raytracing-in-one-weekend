package main

import "core:fmt"
import "core:os"

color :: [3]f64

write_color :: proc (dst : os.Handle, pixel_color : color)
{
	ir := i32(255.999 * pixel_color.x)
	ig := i32(255.999 * pixel_color.y)
	ib := i32(255.999 * pixel_color.z)

	fmt.fprintfln(dst, "%v %v %v", ir, ig, ib)
}

main :: proc ()
{
	// Image

	image_width : i32 = 256
	image_height : i32 = 256

	// Render

	fmt.println("P3")
	fmt.printfln("%v %v", image_width, image_height)
	fmt.println("255")

	for j : i32; j < image_height; j+=1
	{
		fmt.eprintf("\rScanlines remaining: %v", image_height - j)
		for i : i32; i < image_width; i+=1
		{
			pixel_color := color{f64(i) / f64(image_width-1), f64(j) / f64(image_height-1), 0}
			write_color(os.stdout, pixel_color)
		}
	}

	fmt.eprintln("\rDone.                   ")

}