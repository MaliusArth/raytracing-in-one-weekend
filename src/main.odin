package main

import "core:fmt"
import "core:os"
import "core:math"

color :: distinct [3]f64
vec3 :: distinct [3]f64
point3 :: vec3

magnitude :: proc(vec: vec3) -> f64
{
	return math.sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z)
}

normalize :: proc(vec: vec3) -> vec3
{
	vec := vec // explicit mutation via shadowing
	vec = vec / magnitude(vec)
	return vec
}

dot :: proc(a: vec3, b: vec3) -> (result: f64)
{
	result = a.x*b.x + a.y*b.y + a.z*b.z
	return
}

substract :: proc(a: vec3, b: vec3) -> (result: vec3)
{
	result.x = a.x-b.x
	result.y = a.y-b.y
	result.z = a.z-b.z
	return
}

ray :: struct
{
	origin: point3,
	direction: vec3,
}

hit_sphere :: proc(center: ^point3, radius: f64, r: ^ray) -> f64
{
	o_c: vec3 = substract(center^, r.origin)
	a := dot(r.direction, r.direction)
	b := -2.0 * dot(r.direction, o_c)
	c := dot(o_c, o_c) - radius*radius
	discriminant := b*b - 4*a*c
	if discriminant < 0 { return -1 }

	return (-b - math.sqrt(discriminant)) / (2.0*a)
}

write_color :: proc (dst: os.Handle, pixel_color: color)
{
	ir := i32(255.999 * pixel_color.x)
	ig := i32(255.999 * pixel_color.y)
	ib := i32(255.999 * pixel_color.z)

	fmt.fprintfln(dst, "%v %v %v", ir, ig, ib)
}

ray_color :: proc(r: ^ray) -> color
{
	// linear gradient between a and b
	a := color{1.0, 1.0, 1.0}
	b := color{0.5, 0.7, 1.0}
	t := 0.5 * (r.direction.y + 1.0)
	return (1 - t) * a + t * b
}

main :: proc ()
{
	// Image

	aspect_ratio :: 16.0 / 9.0
	image_width :: 400

	// Calculate the image height, and ensure that it's at least 1.
	image_height :: max(1, i64(f64(image_width)/aspect_ratio))

	// Camera

	focal_length :: 1.0
	viewport_height :: 2.0
	viewport_width :: viewport_height * (f64(image_width)/f64(image_height))
	camera_center :: point3{} // default zero initialized

	// Calculate the vectors across the horizontal and down the vertical viewport edges.
	viewport_u :: vec3{viewport_width, 0, 0}
	viewport_v :: vec3{0, -viewport_height, 0}

	// Calculate the horizontal and vertical delta vectors from pixel to pixel.
	pixel_delta_u := viewport_u / f64(image_width)
	pixel_delta_v := viewport_v / f64(image_height)

	// Calculate the location of the upper left pixel.
	viewport_top_left := camera_center - vec3{0, 0, focal_length} - viewport_u/2.0 - viewport_v/2.0
	pixel00_loc := viewport_top_left + 0.5 * (pixel_delta_u + pixel_delta_v)

	// Render

	fmt.println("P3")
	fmt.printfln("%v %v", image_width, image_height)
	fmt.println("255")

	for j in 0..<image_height
	{
		fmt.eprintf("\rScanlines remaining: %v", image_height - j)
		for i in 0..<image_width
		{
			pixel_center := pixel00_loc + (f64(i) * pixel_delta_u) + (f64(j) * pixel_delta_v)
			ray_direction := pixel_center - camera_center
			ray_direction = normalize(ray_direction)
			r := ray{camera_center, ray_direction}

			pixel_color := ray_color(&r)

			sphere_center := point3{0,0,-1}
			sphere_radius :: 0.5
			t := hit_sphere(&sphere_center, sphere_radius, &r)
			if t > 0
			{
				hit_point := (r.origin + t * r.direction)
				normal := normalize(hit_point - sphere_center)
				pixel_color = 0.5*color{normal.x+1, normal.y+1, normal.z+1}
			}

			write_color(os.stdout, pixel_color)
		}
	}

	fmt.eprintln("\rDone.                   ")

}