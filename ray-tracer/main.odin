package main

import "core:fmt"
import "core:os"
import "core:math"
import "core:math/rand"
import "core:strings"

/// vec

vec2 :: distinct [2]f64
vec2i :: distinct [2]i64
vec3 :: distinct [3]f64
color :: vec3
point3 :: vec3

dot :: proc(a: vec3, b: vec3) -> f64 {
	return a.x*b.x + a.y*b.y + a.z*b.z
}

magnitude_squared :: proc(vec: vec3) -> f64 {
	return dot(vec, vec)
}

magnitude :: proc(vec: vec3) -> f64 {
	return math.sqrt(magnitude_squared(vec))
}

normalize :: proc(vec: vec3) -> vec3 {
	return vec / magnitude(vec)
}

is_near_zero :: proc(vec: vec3) -> bool {
	EPSILON :: 1e-8
	return (math.abs(vec.x) <= EPSILON &&
	        math.abs(vec.y) <= EPSILON &&
	        math.abs(vec.z) <= EPSILON)
}

// ![reflection](https://raytracing.github.io/images/fig-1.15-reflection.jpg|width=200)
reflect :: proc(vec, normal: vec3) -> vec3 {
	return vec-2*dot(vec, normal)*normal
}

// src_refractive_index: refractive index of the medium the light is exiting
// dst_refractive_index: refractive index of the medium the light is entering
refract_with_reference_medium :: proc(vec, normal: vec3, src_refractive_index, dst_refractive_index: f64) -> vec3 {
	return refract_with_relative_refractive_index(vec, normal, src_refractive_index / dst_refractive_index)
}

// ![refraction](https://raytracing.github.io/images/fig-1.17-refraction.jpg|width=200)
refract_with_relative_refractive_index :: proc(vec, normal: vec3, relative_refractive_index: f64) -> vec3 {
	// vec & normal are expected to be normalized

	// min() mitigates floating-point errors, valid range: [-1..1]
	cos_theta := math.min(dot(-vec, normal), 1.0)
	ray_out_perpendicular := relative_refractive_index * (vec + cos_theta*normal)

	// max() mitigates floating-point errors, valid range:  [0..1]
	ray_out_parallel := -math.sqrt(math.max(1.0-magnitude_squared(ray_out_perpendicular), 0.0))*normal
	return ray_out_perpendicular + ray_out_parallel
}

refract :: proc { refract_with_reference_medium, refract_with_relative_refractive_index }

/// random

random_vec3 :: proc() -> vec3 {
	return {rand.float64(), rand.float64(), rand.float64()}
}

random_vec3_range :: proc(min, max : f64) -> vec3 {
	return {rand.float64_range(min, max), rand.float64_range(min, max), rand.float64_range(min, max)}
}

// Returns the vector to a random point in the [min,min]-[max,max] square.
random_vec2_range :: proc(min, max : f64) -> vec3 {
	return {rand.float64_range(min, max), rand.float64_range(min, max), 0}
}

random_unit_vector :: proc() -> vec3 {
	// rejection sampling to ensure normal distribution
	for {
		p := random_vec3_range(-1, 1 + math.F64_EPSILON) // excl. max
		length_squared := magnitude_squared(p)
		if 0 < length_squared && length_squared <= 1 {
			return p / math.sqrt(length_squared)
		}
	}
}

random_point_on_hemisphere :: proc(normal : vec3) -> vec3 {
	on_unit_sphere := random_unit_vector()
	if dot(on_unit_sphere, normal) > 0.0 { // In the same hemisphere as the normal
		return on_unit_sphere
	} else {
		return -on_unit_sphere
	}
}

/// raycast

ray :: struct {
	origin : point3,
	direction : vec3,
}

// material_proc :: #type proc(data: rawptr)
material_proc :: #type proc(
	data: rawptr,
	input_ray : ^ray,
	hit : ^hit_record,
) -> (output_ray : ray, attenuation : color)

material :: struct {
	procedure : material_proc,
	data : rawptr,
}

lambertian_data :: struct {
	albedo : color,
}

lambertian_proc :: proc(
	data: rawptr,
	input_ray : ^ray,
	hit : ^hit_record,
) -> (output_ray : ray, attenuation : color) {
	// Lambertian (diffuse) reflectance can either
	//  * always scatter and attenuate light according to its reflectance R,
	//  * or it can sometimes scatter (with probability 1−R) with no attenuation (where a ray that isn't scattered is just absorbed into the material).
	// It could also be a mixture of both those strategies. We will choose to always scatter

	// normal_offset :: 1.001
	output_direction := hit.normal /* * normal_offset */ + random_unit_vector()
	// TODO(viktor): can't we just add some epsilon to the normal? this is just shadow acne all over again
	if is_near_zero(output_direction) {
		output_direction = hit.normal
	}
	output_ray = ray{hit.p, output_direction}

	material := cast(^lambertian_data)data
	attenuation = material.albedo
	// Note the third option: we could scatter with some fixed probability p and have attenuation be albedo/p.
	// p := fixed probability
	// attenuation = material.albedo/p
	return
}

metallic_data :: struct {
	albedo : color,
}

metallic_proc :: proc(
	data: rawptr,
	input_ray : ^ray,
	hit : ^hit_record,
) -> (output_ray : ray, attenuation : color) {
	output_direction := reflect(input_ray.direction, hit.normal)
	output_ray = ray{hit.p, output_direction}

	material := cast(^metallic_data)data
	attenuation = material.albedo
	return
}

hit_record :: struct {
	p : point3,
	normal : vec3,
	t : f64,
}

hit_sphere_ranged :: proc(
	center: ^point3, radius: f64,
	r: ^ray, t_range: struct{min, max: f64},
) -> (value: hit_record, ok: bool) #optional_ok {
	o_c := center^ - r.origin
	a := magnitude_squared(r.direction)
	h := dot(r.direction, o_c)
	c := magnitude_squared(o_c) - radius*radius
	discriminant := h*h - a*c
	if discriminant < 0 do return

	sqrtd := math.sqrt(discriminant)

	// Find the nearest root that lies in the acceptable range
	root := (h - sqrtd) / a
	if root < t_range.min || t_range.max < root {
		root = (h + sqrtd) / a
		if root < t_range.min || t_range.max < root do return
	}

	record : hit_record
	record.t = root
	record.p = r.origin + root*r.direction
	record.normal = (record.p - center^) / radius
	front_face := dot(r.direction, record.normal) < 0
	record.normal = front_face ? record.normal : -record.normal

	return record, true
}

background_color :: proc(r: ^ray) -> color {
	// linear gradient between a and b
	a := color{1.0, 1.0, 1.0}
	b := color{0.5, 0.7, 1.0}
	t := 0.5 * (r.direction.y + 1.0)
	return (1 - t) * a + t * b
}

ray_color :: proc(r: ^ray, bounces : i64, spheres : []sphere) -> color {
	// < 0 NOT <= 0 or do it in the if below
	// if bounces < 0 do return color{0,0,0}

	// raycast
	closest_hit := hit_record{t=math.F64_MAX}
	material : ^material
	RAY_OFFSET :: 0.001 // prevent shadow acne
	for &sphere in spheres {
		if hit, ok := hit_sphere_ranged(&sphere.center, sphere.radius, r, {RAY_OFFSET, closest_hit.t}); ok {
			if hit.t < closest_hit.t {
				closest_hit = hit
				material = sphere.material
			}
		}
	}

	output_color : color
	if closest_hit.t < math.F64_MAX {
		if bounces <= 0 { // don't recurse
			output_color = color{0,0,0}
		} else {
			reflected_ray, attenuation := material.procedure(material.data, r, &closest_hit)
			output_color = attenuation * ray_color(&reflected_ray, bounces-1, spheres)
		}
	} else {
		// NOTE(viktor): only needs to be normalized for the background gradient
		r.direction = normalize(r.direction)
		output_color = background_color(r)
	}
	return output_color
}

///

linear_to_gamma2 :: proc(linear_component : f64) -> f64 {
	return linear_component > 0.0 ? math.sqrt(linear_component) : 0.0
}

write_color :: proc (dst: ^strings.Builder, pixel_color: color) {
	// Apply a linear to gamma transform for gamma 2
	r := linear_to_gamma2(pixel_color.r)
	g := linear_to_gamma2(pixel_color.g)
	b := linear_to_gamma2(pixel_color.b)

	// Translate the [0,1] component values to the byte range [0,255].
	intensity :: vec2{0.000, 0.999}
	ir := i32(256 * clamp(r, intensity.x, intensity.y))
	ig := i32(256 * clamp(g, intensity.x, intensity.y))
	ib := i32(256 * clamp(b, intensity.x, intensity.y))

	fmt.sbprintfln(dst, "%v %v %v", ir, ig, ib)
}

sphere :: struct {
	center : vec3,
	radius : f64,
	material : ^material,
}

camera :: struct {
	position : point3,
	focal_length : f64,
	aspect_ratio : f64,
	image_size : vec2i,
	viewport : vec2,
	samples_per_pixel : i64,
	max_ray_bounces : i64,
}

camera_init :: proc(
	camera : ^camera,
	position := point3{},
	focal_length : f64 = 1.0,
	target_aspect_ratio : f64 = 1.0,
	image_width : i64 = 100,
	// viewport : vec2,
	samples_per_pixel : i64 = 10,
	max_ray_bounces : i64 = 10,
) {
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
	camera.viewport = {viewport_height * camera.aspect_ratio, viewport_height}

	camera.samples_per_pixel = samples_per_pixel
	camera.max_ray_bounces = max_ray_bounces
}

render :: proc(str : ^strings.Builder, camera : camera, spheres : []sphere, $print_progress : bool) {
	// rasterization

	// Calculate the horizontal and vertical delta vectors from pixel to pixel.
	pixel_deltas : vec3 =
	{
		 camera.viewport.x / f64(camera.image_size.x),
		-camera.viewport.y / f64(camera.image_size.y),
		0,
	}

	// Calculate the location of the top left pixel.
	/*
	  NOTE(viktor): this is used to calculate the ray_direction which substracts the camera.position away again
	  so we might as well leave it out altogether
	*/
	viewport_top_left := /* camera.position */ - vec3{camera.viewport.x*0.5, -camera.viewport.y*0.5, camera.focal_length}
	pixel00_center_in_3d := viewport_top_left + 0.5 * pixel_deltas

	fmt.sbprintfln(str, "P3\n%v %v\n255", camera.image_size.x, camera.image_size.y)
	pixel_samples_scale := 1.0 / f64(camera.samples_per_pixel)
	for j in 0..<camera.image_size.y {
		when print_progress do fmt.eprintf("\rScanlines remaining: %v ", camera.image_size.y - j)
		for i in 0..<camera.image_size.x {
			pixel_color : color
			for _ in 0..<camera.samples_per_pixel {
				// get ray

				offset := random_vec2_range(-0.5, 0.5)
				pixel_sample := pixel00_center_in_3d + pixel_deltas * (vec3{f64(i), f64(j), 0} + offset)
				ray_direction := pixel_sample /* - camera.position */
				r := ray{camera.position, ray_direction}

				pixel_color += ray_color(&r, camera.max_ray_bounces, spheres)
			}
			write_color(str, pixel_samples_scale * pixel_color)
		}
	}

	when print_progress do fmt.eprintln("\rDone.                   ")
}

main :: proc () {
	rand.reset(seed=1)

	camera : camera
	camera_init(
		&camera,
		target_aspect_ratio=16.0/9.0,
		image_width=400,
		samples_per_pixel=100,
		max_ray_bounces=50,
	)

	// Scene
	ground_data := lambertian_data{albedo={0.8, 0.8, 0.0}}
	blue_data := lambertian_data{albedo={0.1, 0.2, 0.5}}

	ground := material{procedure=lambertian_proc, data=&ground_data}
	blue   := material{procedure=lambertian_proc, data=&blue_data}
	silver := material{procedure=metallic_proc,   data=  &metallic_data{albedo={0.8, 0.8, 0.8}}}
	gold   := material{procedure=metallic_proc,   data=  &metallic_data{albedo={0.8, 0.6, 0.2}}}

	spheres := []sphere{
		{center={ 0.0, -100.5, -1.0}, radius=100, material=&ground},
		{center={ 0.0,    0.0, -1.2}, radius=0.5, material=&blue},
		{center={-1.0,    0.0, -1.0}, radius=0.5, material=&silver},
		{center={ 1.0,    0.0, -1.0}, radius=0.5, material=&gold},
	}

	PPM_HEADER_SIZE :: 3 + 2 * 4 + 3

	str: strings.Builder
	strings.builder_init(&str, 0, cast(int)(camera.image_size.x * camera.image_size.y * 3 * 4 + PPM_HEADER_SIZE))
	defer strings.builder_destroy(&str)

	// diff : time.Duration
	// {
	// 	time.SCOPED_TICK_DURATION(&diff)
	render(&str, camera, spheres, true)
	// }

	fmt.fprintln(os.stdout, strings.to_string(str))

	// fmt.eprintfln("render took %v", diff)
}
