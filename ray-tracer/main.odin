package main

import "core:fmt"
import "core:os"
import "core:math"

vec2 :: distinct [2]f64
vec2i :: distinct [2]i64
vec3 :: distinct [3]f64
color :: vec3
point3 :: vec3

dot :: proc(a: vec3, b: vec3) -> f64
{
	return a.x*b.x + a.y*b.y + a.z*b.z
}

magnitude_squared :: proc(vec: vec3) -> f64
{
	return dot(vec, vec)
}

magnitude :: proc(vec: vec3) -> f64
{
	return math.sqrt(magnitude_squared(vec))
}

normalize :: proc(vec: vec3) -> vec3
{
	return vec / magnitude(vec)
}

ray :: struct
{
	origin : point3,
	direction : vec3,
}

hit_record :: struct
{
	p : point3,
	normal : vec3,
	t : f64,
	front_face : b8,
}

hit_sphere_ranged :: proc(center: ^point3, radius: f64, r: ^ray, t_range: struct{min, max: f64}) -> (value: hit_record, ok: bool) #optional_ok
{
	o_c := center^ - r.origin
	a := magnitude_squared(r.direction)
	h := dot(r.direction, o_c)
	c := magnitude_squared(o_c) - radius*radius
	discriminant := h*h - a*c
	if discriminant < 0 do return

	sqrtd := math.sqrt(discriminant)

	// Find the nearest root that lies in the acceptable range
	root := (h - sqrtd) / a
	if root < t_range.min || t_range.max < root
	{
		root = (h + sqrtd) / a
		if root < t_range.min || t_range.max < root do return
	}

	record : hit_record
	record.t = root
	record.p = r.origin + root*r.direction
	record.normal = (record.p - center^) / radius
	record.front_face = dot(r.direction, record.normal) < 0

	// adjust normal to ray being inside/outside the sphere
	record.normal = record.front_face ? record.normal : -record.normal

	return record, true
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

sphere :: struct
{
	center : vec3,
	radius : f64,
}

camera :: struct
{
	position : point3,

	focal_length : f64,
	aspect_ratio : f64,
	image_size : vec2i,
	viewport : vec2,
}

create_camera :: proc(
	position := point3{},
	focal_length : f64 = 1.0,
	target_aspect_ratio : f64 = 1.0,
	image_width : i64 = 100,
	// viewport : vec2,
) -> (camera : camera)
{
	// Image
	camera.position = position

	camera.focal_length = focal_length
	camera.image_size.x = image_width

	// Calculate the image height, and ensure that it's at least 1.
	camera.image_size.y = max(1, i64(f64(camera.image_size.x)/target_aspect_ratio))
	// TODO(viktor): instead of adjusting the ratio, adjust width?
	// adjust ratio in case height had to be overwritten to 1
	camera.aspect_ratio = f64(camera.image_size.x)/f64(camera.image_size.y)

	// TODO(viktor): why is this hardcoded and why height?
	viewport_height := 2.0
	camera.viewport = { viewport_height * camera.aspect_ratio, viewport_height }

	return
}

render :: proc(camera : camera, spheres : []sphere)
{
	// 1280 × 720
	// 16/9 = 1.778

	// rasterization

	// Calculate the horizontal and vertical delta vectors from pixel to pixel.
	pixel_deltas : vec3 = {
		camera.viewport.x / f64(camera.image_size.x), // 3.556/1280=0.002778
		-camera.viewport.y / f64(camera.image_size.y), // -2/720=-0.002778
		0,
	}

	// Calculate the location of the top left pixel.
	// {0, 0, -1} - {1.778, 0, 0} - {0, -1, 0} -> {-1.778, -1, -1}
	// NOTE(viktor): this is used to calculate the ray_direction which substracts the camera.position away again so we might as well leave it out altogether
	viewport_top_left := /* camera.position */ - vec3{camera.viewport.x*0.5, -camera.viewport.y*0.5, camera.focal_length}
	pixel00_center_in_3d := viewport_top_left + 0.5 * pixel_deltas

	fmt.printfln("P3\n%v %v\n255", camera.image_size.x, camera.image_size.y)

	for j in 0..<camera.image_size.y
	{
		fmt.eprintf("\rScanlines remaining: %v", camera.image_size.y - j)
		for i in 0..<camera.image_size.x
		{
			pixel_center := pixel00_center_in_3d + pixel_deltas * vec3{f64(i), f64(j), 0}
			ray_direction := pixel_center /* - camera.position */
			// ray_direction = normalize(ray_direction)
			r := ray{camera.position, ray_direction}

			closest_t := math.inf_f64(0)
			// NOTE(viktor): only needs to be normalized for the background gradient
			r.direction = normalize(r.direction)
			pixel_color := ray_color(&r) // NOTE(viktor): background color (gradient)

			for &sphere in spheres
			{
				if record, ok := hit_sphere_ranged(&sphere.center, sphere.radius, &r, {0, closest_t}); ok
				{
					if record.t < closest_t
					{
						closest_t = record.t
						pixel_color = 0.5*color(record.normal+1)
					}
				}
			}
			write_color(os.stdout, pixel_color)
		}
	}

	fmt.eprintln("\rDone.                   ")
}

main :: proc ()
{
	camera := create_camera(target_aspect_ratio=16.0/9.0, image_width=400)

	// Scene
	spheres :: []sphere{
		{center={0, 0, -1}, radius=0.5},
		// {center={0.2, 0.5, -1}, radius=0.2},
		// {center={0.3, -0.1, -0.7}, radius=0.2},
		// {center={-0.5, -0.25, -0.5}, radius=0.2},
		// {center={0.65, 0.6, -0.8}, radius=0.3},
		{center={0, -100.5, -1}, radius=100},
	}

	render(camera, spheres)

}