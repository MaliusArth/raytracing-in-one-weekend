package main

import "core:fmt"
import "core:os"
import "core:math"
import "core:math/rand"
import "core:strings"

/// vec

f64x2  :: [2]f64
f64x3  :: [3]f64
i64x2  :: [2]i64
intx2  :: [2]int
v3     :: f64x3
vec3   :: f64x3
point3 :: f64x3
color  :: f64x3
v2     :: f64x2
vec2   :: f64x2
turns  :: f64 // [0, 1]
// turns :: distinct f64 // [0, 1]
// radians :: distinct f64 // [0, math.TAU]
turns_to_radians :: proc(θ: turns) -> f64 { return θ * math.TAU }


/// linalg

dot :: proc(a, b: vec3) -> f64 {
	return a.x*b.x + a.y*b.y + a.z*b.z
}

cross :: proc(a, b: vec3) -> vec3 {
	result: vec3
	result[0] = a[1]*b[2] - b[1]*a[2]
	result[1] = a[2]*b[0] - b[2]*a[0]
	result[2] = a[0]*b[1] - b[0]*a[1]
	return result
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

is_normalized :: proc(vec: vec3) -> bool {
	EPSILON :: 1e-10
	mag2 := magnitude_squared(vec)
	return 1-EPSILON <= mag2 && mag2 <= 1+EPSILON
}

lookat :: proc(position: point3 = {}, target: point3 = {0, 0, -1}, axis_up: vec3 = {0, 1, 0}) -> (right, up, forward: vec3) {

	// view cartesian coordinate system
	forward = normalize(target - position)
	right   = normalize(cross(forward, axis_up))
	up      = cross(right, forward)
	return
}

// ![reflection](https://raytracing.github.io/images/fig-1.15-reflection.jpg|width=200)
reflect :: proc(vec, normal: vec3) -> vec3 {
	// vec & normal don't need to be normalized
	return vec-2*dot(vec, normal)*normal
}

// * src_refractive_index: refractive index of the medium the light is exiting, aka. incident refractive index
// * dst_refractive_index: refractive index of the medium the light is entering, aka. transmitted refractive index
refract_with_reference_medium :: proc(vec, normal: vec3, src_refractive_index, dst_refractive_index: f64) -> vec3 {
	return refract_with_relative_refractive_index(vec, normal, src_refractive_index / dst_refractive_index)
}

// ![refraction](https://raytracing.github.io/images/fig-1.17-refraction.jpg|width=200)
refract_with_relative_refractive_index :: proc(vec, normal: vec3, relative_refractive_index: f64) -> vec3 {
	// contract: vec & normal are expected to be normalized
	// contract: can_refract := (rel_refractive_index * sin_theta) <= 1.0
	fmt.assertf(is_normalized(vec),    "vector needs to be normalized! %v", magnitude_squared(vec))
	fmt.assertf(is_normalized(normal), "vector needs to be normalized! %v", magnitude_squared(normal))
	when !ODIN_DISABLE_ASSERT {
	{
		cos_theta := math.min(dot(-vec, normal), 1.0)
		sin_theta := math.sqrt(1.0 - cos_theta*cos_theta)
		can_refract := relative_refractive_index * sin_theta <= 1.0
		assert(can_refract)
	}
	}

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

hit_record :: struct {
	p : point3,
	normal : vec3,
	t : f64,
	front_face : bool,
}

lambertian_data :: struct {
	albedo : color,
}

lambertian_data_init :: proc(data : ^lambertian_data, albedo : color) {
	data^ = {albedo}
}

lambertian_data_make :: proc(albedo : color) -> (result : lambertian_data) {
	lambertian_data_init(&result, albedo)
	return
}

// ![](https://raytracing.github.io/images/fig-1.14-rand-unitvec.jpg|width=200)
lambertian_proc :: proc(data: rawptr, ray_in: ^ray, hit: ^hit_record) ->
                       (ray_out: ray, attenuation: color, ok: bool) {
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
	ray_out = ray{hit.p, output_direction}

	material := cast(^lambertian_data)data
	attenuation = material.albedo
	// Note the third option: we could scatter with some fixed probability p and have attenuation be albedo/p.
	// p := fixed probability
	// attenuation = material.albedo/p
	ok = true
	return
}

metallic_data :: struct {
	albedo : color,
	fuzz : f64,
}

metallic_data_init :: proc(data: ^metallic_data, albedo: color, fuzz: f64) {
	data^ = {albedo, math.clamp(fuzz, 0.0, 1.0)}
}

metallic_data_make :: proc(albedo: color, fuzz: f64) -> (result: metallic_data) {
	metallic_data_init(&result, albedo, fuzz)
	return
}

// ![](https://raytracing.github.io/images/fig-1.16-reflect-fuzzy.jpg|width=200)
metallic_proc :: proc(data: rawptr, ray_in: ^ray, hit: ^hit_record) ->
                     (ray_out: ray, attenuation: color, ok: bool) {
	material := cast(^metallic_data)data
	output_direction := normalize(reflect(ray_in.direction, hit.normal))
	output_direction += material.fuzz * random_unit_vector()
	ray_out = ray{hit.p, output_direction}

	attenuation = material.albedo
	ok = dot(output_direction, hit.normal) > 0
	return
}

dielectric_data :: struct {
	refractive_index: f64,
}

dielectric_data_init :: proc(data: ^dielectric_data, refraction_index: f64) {
	data^ = {refraction_index}
}

dielectric_data_make :: proc(refraction_index: f64) -> (result: dielectric_data) {
	dielectric_data_init(&result, refraction_index)
	return
}

// ![](https://raytracing.github.io/images/fig-1.17-refraction.jpg|width=200)
dielectric_proc :: proc(data: rawptr, ray_in: ^ray, hit: ^hit_record) -> (ray_out: ray, attenuation: color, ok: bool) {
	material := cast(^dielectric_data)data

	// materials with a refractive_index < 1 (air/vaccuum) are treated as air/vaccuum materials
	// inside a denser material with the effective refractive index of 1/material.refractive_index
	src_refractive_index := 1.0                       if material.refractive_index >= 1.0 else 1.0 / material.refractive_index
	dst_refractive_index := material.refractive_index if material.refractive_index >= 1.0 else 1.0

	if !hit.front_face {
		// swap
		tmp := src_refractive_index
		src_refractive_index = dst_refractive_index
		dst_refractive_index = tmp
	}

	rel_refractive_index := src_refractive_index / dst_refractive_index

	unit_direction := normalize(ray_in.direction)
	cos_theta := dot(-unit_direction, hit.normal)
	// math.saturate/math.clamp [0..1]
	cos_theta  = math.min(cos_theta, 1.0)
	// assert(cos_theta >= 0.0)
	sin_theta := math.sqrt(1.0-cos_theta*cos_theta)

	reflectance_fresnel :: proc(cos_i, sin_i, src_refractive_index, dst_refractive_index: f64) -> f64 {
		n1 := src_refractive_index
		n2 := dst_refractive_index

		x := (n1 / n2) * sin_i
		x *= x
		x_invpow := 1.0 - x
		subexp := math.sqrt(x_invpow)

		x1  := n1 * cos_i
		x2  := n2 * subexp
		R_s := (x1 - x2)
		R_s /= (x1 + x2)
		R_s *= R_s

		x1   = n1 * subexp
		x2   = n2 * cos_i
		R_p := (x1 - x2)
		R_p /= (x1 + x2)
		R_p *= R_p

		return (R_s + R_p) * 0.5 // avg
	}

	reflectance_schlicks_approximation :: proc(cos_i, rel_refractive_index: f64) -> f64 {
		r0 := (1.0 - rel_refractive_index) / (1.0 + rel_refractive_index)
		r0 *= r0

		a := 1.0 - cos_i
		a = math.pow(a, 5)
		return math.lerp(a, 1.0, r0)
	}
	reflectance :: proc{reflectance_schlicks_approximation, reflectance_fresnel}

	// NOTE: total internal reflection:
	// When a ray interfaces with the surface of a less dense external medium at an angle beyond a certain critical angle
	// it is completely reflected inside the internal medium.

	// determine if angle is beyond critical angle for total internal reflection
	must_reflect := (rel_refractive_index * sin_theta) > 1.0
	output_direction: vec3

	if must_reflect || reflectance(cos_theta, rel_refractive_index) > rand.float64() {
		output_direction = reflect(unit_direction, hit.normal)
	} else {
		output_direction = refract(unit_direction, hit.normal, rel_refractive_index)
	}
	ray_out = ray{hit.p, output_direction}

	attenuation = {1,1,1}
	ok = true
	return
}

material_proc :: #type proc(data: rawptr, ray_in: ^ray, hit: ^hit_record) -> (ray_out: ray, attenuation: color, ok: bool)

material :: struct {
	procedure : material_proc,
	data : rawptr,
}

material_make_lambertian :: proc(data: ^lambertian_data) -> material { return {lambertian_proc, data} }
material_make_metallic   :: proc(data: ^metallic_data  ) -> material { return {metallic_proc,   data} }
material_make_dielectric :: proc(data: ^dielectric_data) -> material { return {dielectric_proc, data} }
material_make :: proc {
	material_make_lambertian,
	material_make_metallic,
	material_make_dielectric,
}

hit_sphere_ranged :: #force_inline proc(center: ^point3, radius: f64, r: ^ray, t_range: struct{min, max: f64}) ->
                                       (record: hit_record, ok: bool) {
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

	record.t = root
	record.p = r.origin + root*r.direction
	record.normal = (record.p - center^) / radius
	record.front_face = dot(r.direction, record.normal) < 0
	record.normal = record.front_face ? record.normal : -record.normal

	return record, true
}

background_color :: proc(r: ^ray) -> color {
	// linear gradient between a and b
	a := color{1.0, 1.0, 1.0}
	b := color{0.5, 0.7, 1.0}
	t := 0.5 * (r.direction.y + 1.0)
	return math.lerp(a, b, t)
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
			output_color = {0,0,0}
		} else {
			if reflected_ray, attenuation, ok := material.procedure(material.data, r, &closest_hit); ok {
				output_color = attenuation * ray_color(&reflected_ray, bounces-1, spheres)
			} else {
				output_color = {0,0,0}
			}
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
	// orientation : quaternion256,
	right: vec3,
	up: vec3,
	forward: vec3,
	focal_length : f64,
	aspect_ratio : f64,
	image_size : vec2,
	viewport_size : vec2,
	samples_per_pixel : i64,
	max_ray_bounces : i64,
	vfov : turns,
}

camera_init :: proc(
	camera : ^camera,
	position : point3 = {},
	right : vec3 = {1, 0, 0},
	up : vec3 = {0, 1, 0},
	forward : vec3 = {0, 0, -1},
	// orientation : quaternion256 = {},
	focal_length : f64 = 1.0,
	target_aspect_ratio : f64 = 1.0,
	image_width : i64 = 100,
	// viewport : vec2,
	samples_per_pixel : i64 = 10,
	max_ray_bounces : i64 = 10,
	vfov : turns = 0.25, // 1/4 turn
) {
	camera.position = position
	// camera.orientation = orientation
	camera.right = right
	camera.up = up
	camera.forward = forward
	camera.focal_length = focal_length

	camera.image_size.x = cast(f64)image_width
	camera.image_size.y = max(1, cast(f64)image_width / target_aspect_ratio)
	// TODO(viktor): instead of adjusting the ratio, adjust width?
	// adjust ratio in case height had to be overwritten to 1
	camera.aspect_ratio = camera.image_size.x / camera.image_size.y

	// TODO(viktor): I don't like the use of vertical fov, the world is mostly horizontal, hfov is more intuitive
	// ![](https://raytracing.github.io/images/fig-1.18-cam-view-geom.jpg|width=200)
	unit_height := 2 * math.tan(turns_to_radians(vfov / 2))
	viewport_height := focal_length * unit_height
	viewport_width  := viewport_height * camera.aspect_ratio
	camera.viewport_size = {viewport_width, viewport_height}

	camera.samples_per_pixel = samples_per_pixel
	camera.max_ray_bounces = max_ray_bounces
}

render :: proc(str : ^strings.Builder, camera : camera, spheres : []sphere, $print_progress : bool) {
	// view-space
	pixel_deltas_vs := camera.viewport_size / camera.image_size * vec2{1, -1} // vertical flip

	// world-space
	pixel_delta_u := camera.right * pixel_deltas_vs.x
	pixel_delta_v := camera.up    * pixel_deltas_vs.y

	// view-space
	view_plane_top_left_pixel_center := camera.viewport_size * 0.5 * vec2{-1, 1} + pixel_deltas_vs * 0.5

	// world-space
	view_plane_top_left_pixel_center_ws :=
		camera.position +
		camera.forward  * camera.focal_length +
		camera.right    * view_plane_top_left_pixel_center.x +
		camera.up       * view_plane_top_left_pixel_center.y

	image_width  := cast(int)camera.image_size.x
	image_height := cast(int)camera.image_size.y
	pixel_sample_contribution_scale := 1.0 / f64(camera.samples_per_pixel)
	for v in 0..<image_height {
		when print_progress do fmt.eprintf("\rScanlines remaining: %v ", image_height - v)
		for u in 0..<image_width {
			pixel_color : color
			for _ in 0..<camera.samples_per_pixel {
				offset := random_vec2_range(-0.5, 0.5)
				pixel_sample :=
					view_plane_top_left_pixel_center_ws +
					pixel_delta_u * (cast(f64)u + offset.x) +
					pixel_delta_v * (cast(f64)v + offset.y)

				ray_direction := pixel_sample - camera.position
				r := ray{camera.position, ray_direction}
				pixel_color += ray_color(&r, camera.max_ray_bounces, spheres)
			}
			write_color(str, pixel_sample_contribution_scale * pixel_color)
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
		position={-2, 2, 1},
		vfov=20.0/360.0,
	)

	// Scene
	ground := material_make(&lambertian_data{albedo={0.8, 0.8, 0.0}})
	blue   := material_make(&lambertian_data{albedo={0.1, 0.2, 0.5}})
	// silver := material_make(  &metallic_data{albedo={0.8, 0.8, 0.8}, fuzz=0.3})
	gold   := material_make(  &metallic_data{albedo={0.8, 0.6, 0.2}, fuzz=1.0})
	glass  := material_make(&dielectric_data{refractive_index=1.5})
	air_bubble := material_make(&dielectric_data{refractive_index=1.0/1.5})

	// black  := material_make(&lambertian_data{albedo={0.1, 0.1, 0.1}})
	// red    := material_make(&lambertian_data{albedo={0.5, 0.2, 0.1}})
	// green  := material_make(&lambertian_data{albedo={0.2, 0.5, 0.1}})

	spheres := []sphere{
		{center={ 0.0, -100.5, -1.0}, radius=100, material=&ground},
		{center={ 0.0,    0.0, -1.2}, radius=0.5, material=&blue},
		{center={-1.0,    0.0, -1.0}, radius=0.5, material=&glass},
		{center={-1.0,    0.0, -1.0}, radius=0.4, material=&air_bubble},
		{center={ 1.0,    0.0, -1.0}, radius=0.5, material=&gold},

		// {center={ 0.0,    0.0,  0.0}, radius=0.1,  material=&black},
		// {center={ 0.0,    0.0,  0.0}, radius=0.05, material=&red},
		// {center={ 0.0,    0.0,  0.0}, radius=0.05, material=&green},
		// {center={ 0.0,    0.0,  0.0}, radius=0.05, material=&blue},
	}

	// black_ball   := &spheres[5]
	// right_ball   := &spheres[6]
	// up_ball      := &spheres[7]
	// forward_ball := &spheres[8]

	// right_ball.center, up_ball.center, forward_ball.center =
	// 	lookat(position=black_ball.center, target=spheres[4].center, axis_up={0, 1, 0})
	// right_ball.center *= 0.1
	// up_ball.center *= 0.1
	// forward_ball.center *= 0.1
	camera.right, camera.up, camera.forward = lookat(position=camera.position, target={0, 0, -1}, axis_up={0, 1, 0})

	str: strings.Builder
	header := fmt.aprintfln("P3\n%v %v\n255", cast(int)camera.image_size.x, cast(int)camera.image_size.y)
	defer delete_string(header)
	strings.builder_init(&str, 0, len(header) + cast(int)camera.image_size.x * cast(int)camera.image_size.y * 3 * 4)
	fmt.sbprint(&str, header)
	defer strings.builder_destroy(&str)

	render(&str, camera, spheres, true)

	fmt.fprint(os.stdout, strings.to_string(str))
}
