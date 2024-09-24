package main

import "core:fmt"

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
			r := f64(i) / f64(image_width-1)
			g := f64(j) / f64(image_height-1)
			b := 0.0

			ir := i32(255.999 * r)
			ig := i32(255.999 * g)
			ib := i32(255.999 * b)

			fmt.printfln("%v %v %v", ir, ig, ib)
		}
	}

	fmt.eprintln("\rDone.                   ")

}