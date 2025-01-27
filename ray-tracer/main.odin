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

// cite: Ray Tracing Gems II (pp.109-114): Majercik, Zander. (2021). The Schlick Fresnel Approximation. 10.1007/978-1-4842-7185-8_9.
// ref: https://www.researchgate.net/publication/354065225_The_Schlick_Fresnel_Approximation
reflectance_fresnel :: proc(cos_i, sin_i, src_refractive_index, dst_refractive_index: f64) -> f64 {
	n1 := src_refractive_index
	n2 := dst_refractive_index

	subexp := (n1 / n2) * sin_i
	subexp  = 1.0 - subexp * subexp
	subexp  = math.sqrt(subexp)
	// assert(!math.is_nan(subexp))

	// parallel polarization
	x1  := n1 * cos_i
	x2  := n2 * subexp
	R_s := (x1 - x2)
	R_s /= (x1 + x2)
	R_s *= R_s

	// perpendicular polarization
	x1   = n1 * subexp
	x2   = n2 * cos_i
	R_p := (x1 - x2)
	R_p /= (x1 + x2)
	R_p *= R_p

	// we ignore light polarization and average the equations
	return (R_s + R_p) * 0.5
}

// cite: Schlick, C. An inexpensive BRDF model for physically-based rendering.
// ref: https://onlinelibrary.wiley.com/doi/10.1111/1467-8659.1330233
reflectance_schlick_approximation :: proc(cos_i, rel_refractive_index: f64) -> f64 {
	// calculate reflectance at normal incidence
	r0 := (1.0 - rel_refractive_index) / (1.0 + rel_refractive_index)
	r0 *= r0

	a  := 1.0 - cos_i
	return r0 + (1.0 - r0) * a*a*a*a*a
	// this is essentially a lerp: a^5*(1-r0) + 1.0*r0
	// return math.lerp(a*a*a*a*a, 1.0, r0)
	// r0 & a^5 are interchangeable since these terms are the same
	// a^5 * 1.0 + (1.0 - a^5) * r0 == r0 * 1.0 + (1.0 - r0) * a^5
	// so we arrive at a lerp by a function of cos_i
	// return math.lerp(r0, 1.0, a*a*a*a*a)
}

// cite: Lazányi, István & Szirmay-Kalos, László. (2005). Fresnel Term Approximations for Metals.. 77-80.
// ref: https://www.researchgate.net/publication/221546550_Fresnel_Term_Approximations_for_Metals
reflectance_schlick_lazanyi_approximation :: proc(cos_i, rel_refractive_index: f64, a: f64, alpha: f64) -> f64 {
	// return reflectance_schlick_approximation(r0, cos_i) - a * cos_i * math.pow(1 - cos_i, alpha)
	return reflectance_schlick_approximation(rel_refractive_index, cos_i) - a * cos_i * math.pow(1 - cos_i, alpha)
}

// reflectance_schlick_lazanyi_approximation :: proc(r0: vec3, cos_i: float, a: vec3, alpha: float) -> f64 {
// 	return schlickFresnel(r0, cos_i) - a * cos_i * math.pow(1 - cos_i, alpha);
// }

// cite: Hoffman, N. Fresnel equations considered harmful. In Eurographics Workshop on Material Appearance Modeling, pages 7–11, 2019. DOI: 10.2312/mam.20191305.
// ref: https://diglib.eg.org/server/api/core/bitstreams/726dc384-d7dd-4c0e-8806-eadec0ff3886/content
reflectance_hoffman_approximation :: proc(cos_i, rel_refractive_index: f64, h: f64) -> f64 {
	// HACK(viktor): we calculate r0 twice, extract it out of the schlicks impl instead, this would also match the papers
	// calculate reflectance at normal incidence
	r0 := (1.0 - rel_refractive_index) / (1.0 + rel_refractive_index)
	r0 *= r0

	a := 823543 / 46656 * (r0 - h) + 49 / 6 * (1 - r0)
	return reflectance_schlick_lazanyi_approximation(cos_i, rel_refractive_index, a, alpha = 6)
}
// reflectance :: proc{reflectance_fresnel, reflectance_schlick_approximation, reflectance_schlick_lazanyi_approximation, reflectance_hoffman_approximation}

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

random_point_on_disk :: proc() -> point3 {
	// rejection sampling to ensure normal distribution
	for {
		p := random_vec2_range(-1, 1 + math.F64_EPSILON) // excl. max
		length_squared := magnitude_squared(p)
		if length_squared <= 1 {
			return p
		}
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
	data^ = {albedo=albedo}
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
	assert(0 <= fuzz && fuzz <= 1.0)
	data^ = {albedo=albedo, fuzz=fuzz}
	// data^ = {albedo, math.clamp(fuzz, 0.0, 1.0)}
}

// ![](https://raytracing.github.io/images/fig-1.16-reflect-fuzzy.jpg|width=200)
metallic_proc :: proc(data: rawptr, ray_in: ^ray, hit: ^hit_record) ->
                     (ray_out: ray, attenuation: color, ok: bool) {
	material := cast(^metallic_data)data
	output_direction := normalize(reflect(ray_in.direction, hit.normal))
	output_direction += material.fuzz * random_unit_vector()
	ray_out = ray{hit.p, output_direction}

	USE_METAL_FRESNEL :: #config(USE_METAL_FRESNEL, false)
	when USE_METAL_FRESNEL {
		cos_theta := dot(-normalize(ray_in.direction), hit.normal)
		// math.saturate/math.clamp [0..1]
		cos_theta  = math.min(cos_theta, 1.0)

		HARDCODED_REFRACTIVE_INDEX :: 1.27035
		// GOLD_REFRACTIVE_INDEX :: complex(cast(f64)0.27035, cast(f64)2.7790)
		METAL_FRESNEL_KIND :: 2
		when METAL_FRESNEL_KIND == 0 {
			reflection_factor := reflectance_schlick_approximation(cos_theta, 1.0/HARDCODED_REFRACTIVE_INDEX)
		} else when METAL_FRESNEL_KIND == 1 {
			H :: 0.5
			reflection_factor := reflectance_hoffman_approximation(cos_theta, 1.0/HARDCODED_REFRACTIVE_INDEX, H)
		} else {
			sin_theta_squared := 1.0-cos_theta*cos_theta
			sin_theta := math.sqrt(sin_theta_squared)
			reflection_factor := reflectance_fresnel(cos_theta, sin_theta, 1.0, HARDCODED_REFRACTIVE_INDEX)
		}

		// refraction_color*refraction_factor + reflection_color*reflection_factor
		attenuation = math.lerp(material.albedo, color{1, 1, 1}, reflection_factor)
	} else {
		attenuation = material.albedo
	}
	ok = dot(output_direction, hit.normal) > 0
	return
}

dielectric_data :: struct {
	refractive_index: f64,
}

dielectric_data_init :: proc(data: ^dielectric_data, refractive_index: f64) {
	data^ = {refractive_index=refractive_index}
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
	sin_theta_squared := 1.0-cos_theta*cos_theta
	// sin_theta := math.sqrt(sin_theta_squared)

	// NOTE: total internal reflection:
	// When a ray interfaces with the surface of a less dense external medium at an angle beyond a certain critical angle
	// it is completely reflected inside the internal medium.

	// determine if angle is beyond critical angle for total internal reflection
	must_reflect := (rel_refractive_index*rel_refractive_index * sin_theta_squared) > 1.0
	// must_reflect := (rel_refractive_index * sin_theta) > 1.0
	output_direction: vec3

	if must_reflect || reflectance_schlick_approximation(cos_theta, rel_refractive_index) > rand.float64() {
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

material_make_lambertian :: proc(albedo: color, allocator := context.allocator) -> material {
	data := new(lambertian_data, allocator)
	lambertian_data_init(data, albedo)
	return {lambertian_proc, data}
}

material_make_metallic :: proc(albedo: color, fuzz: f64, allocator := context.allocator) -> material {
	data := new(metallic_data, allocator)
	metallic_data_init(data, albedo, fuzz)
	return {metallic_proc, data}
}

material_make_dielectric :: proc(refractive_index: f64, allocator := context.allocator) -> material {
	data := new(dielectric_data, allocator)
	dielectric_data_init(data, refractive_index)
	return {dielectric_proc, data}
}

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

ray_cast :: proc(r: ^ray, bounces : i64, spheres : []sphere) -> color {
	// < 0 NOT <= 0 or do it in the if below
	// if bounces < 0 do return color{0,0,0}

	// raycast
	closest_hit := hit_record{t=math.F64_MAX}
	material : material
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
				output_color = attenuation * ray_cast(&reflected_ray, bounces-1, spheres)
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
	material : material,
}

camera :: struct {
	position : point3,
	right: vec3,
	up: vec3,
	forward: vec3,
	// orientation : quaternion256,
	aspect_ratio : f64,
	image_size : vec2,
	focus_distance: f64,
	vfov : turns,
	depth_of_field_angle: turns,
	samples_per_pixel : i64,
	max_ray_bounces : i64,
}

render :: proc(str : ^strings.Builder, camera : camera, spheres : []sphere, $print_progress : bool) {
	// TODO(viktor): I don't like the use of vertical fov, the world is mostly horizontal, hfov is more intuitive
	// ![](https://raytracing.github.io/images/fig-1.18-cam-view-geom.jpg|width=200)
	// ![](https://learn.microsoft.com/en-us/windows/uwp/graphics-concepts/images/fovdiag.png|width=200)
	// we position the view plane on the focus plane
	view_plane_half_size: v2
	view_plane_half_size.y = camera.focus_distance * math.tan(turns_to_radians(camera.vfov * 0.5))
	view_plane_half_size.x = view_plane_half_size.y * camera.aspect_ratio

	// view-space
	pixel_deltas_vs := view_plane_half_size * 2 / camera.image_size * {1.0, -1.0} // vertical flip

	// world-space
	pixel_delta_u := camera.right * pixel_deltas_vs.x
	pixel_delta_v := camera.up    * pixel_deltas_vs.y

	// view-space
	view_plane_top_left_pixel_center_vs := view_plane_half_size * {-1.0, 1.0} + pixel_deltas_vs * 0.5

	// world-space
	view_plane_top_left_pixel_center :=
		camera.position +
		camera.forward  * camera.focus_distance +
		camera.right    * view_plane_top_left_pixel_center_vs.x +
		camera.up       * view_plane_top_left_pixel_center_vs.y

	dof_radius := camera.focus_distance * math.tan(turns_to_radians(camera.depth_of_field_angle*0.5))
	dof_disk_u := camera.right * dof_radius
	dof_disk_v := camera.up    * dof_radius

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
					view_plane_top_left_pixel_center +
					(cast(f64)u + offset.x) * pixel_delta_u +
					(cast(f64)v + offset.y) * pixel_delta_v

				ray_origin := camera.position
				if camera.depth_of_field_angle > 0.0 {
					dof_sample := random_point_on_disk()
					ray_origin +=
						(dof_sample.x * dof_disk_u) +
						(dof_sample.y * dof_disk_v)
				}
				ray_direction := pixel_sample - ray_origin
				r := ray{ray_origin, ray_direction}
				pixel_color += ray_cast(&r, camera.max_ray_bounces, spheres)
			}
			write_color(str, pixel_sample_contribution_scale * pixel_color)
		}
	}

	when print_progress do fmt.eprintln("\rDone.                   ")
}

build_dev_scene :: proc(allocator := context.allocator) -> (camera, [dynamic]sphere) {
	spheres := make_dynamic_array([dynamic]sphere, allocator)
	ground     := material_make_lambertian({0.8, 0.8, 0.0}, allocator)
	blue       := material_make_lambertian({0.1, 0.2, 0.5}, allocator)
	// silver     := material_make_metallic({0.8, 0.8, 0.8}, 0.3, allocator)
	glass      := material_make_dielectric(1.5, allocator)
	air_bubble := material_make_dielectric(1.0/1.5, allocator)
	gold       := material_make_metallic({0.8, 0.6, 0.2}, 1.0, allocator)

	append(&spheres, sphere{center={ 0.0, -100.5, -1.0}, radius=100, material=ground})
	append(&spheres, sphere{center={ 0.0,    0.0, -1.2}, radius=0.5, material=blue})
	append(&spheres, sphere{center={-1.0,    0.0, -1.0}, radius=0.5, material=glass})
	append(&spheres, sphere{center={-1.0,    0.0, -1.0}, radius=0.4, material=air_bubble})
	append(&spheres, sphere{center={ 1.0,    0.0, -1.0}, radius=0.5, material=gold})

	// black  := material_make_lambertian(albedo={0.1, 0.1, 0.1})
	// red    := material_make_lambertian(albedo={0.5, 0.2, 0.1})
	// green  := material_make_lambertian(albedo={0.2, 0.5, 0.1})

	// append(&spheres, sphere{center={ 0.0,    0.0,  0.0}, radius=0.1,  material=black})
	// append(&spheres, sphere{center={ 0.0,    0.0,  0.0}, radius=0.05, material=red})
	// append(&spheres, sphere{center={ 0.0,    0.0,  0.0}, radius=0.05, material=green})
	// append(&spheres, sphere{center={ 0.0,    0.0,  0.0}, radius=0.05, material=blue})

	// black_ball   := &spheres[5]
	// right_ball   := &spheres[6]
	// up_ball      := &spheres[7]
	// forward_ball := &spheres[8]

	// right_ball.center, up_ball.center, forward_ball.center =
	// 	lookat(position=black_ball.center, target=spheres[4].center, axis_up={0, 1, 0})
	// right_ball.center *= 0.1
	// up_ball.center *= 0.1
	// forward_ball.center *= 0.1

	camera: camera
	camera.position = {-2, 2, 1}
	camera.right, camera.up, camera.forward = lookat(position=camera.position, target={0, 0, -1}, axis_up={0, 1, 0})
	camera.aspect_ratio = 16.0/9.0
	camera.image_size.x = 400
	camera.image_size.y = camera.image_size.x/camera.aspect_ratio
	camera.focus_distance = 3.4
	camera.vfov = 20.0/360.0
	camera.depth_of_field_angle = 10.0/360.0
	camera.samples_per_pixel = 100
	camera.max_ray_bounces = 50

	return camera, spheres
}

main :: proc () {
	rand.reset(1)

	// Scene
	camera, spheres := build_dev_scene(context.allocator)
	defer for sphere in spheres { free(sphere.material.data) }
	defer delete_dynamic_array(spheres)

	str: strings.Builder
	header := fmt.aprintfln("P3\n%v %v\n255", cast(int)camera.image_size.x, cast(int)camera.image_size.y)
	defer delete_string(header)
	strings.builder_init(&str, 0, len(header) + cast(int)camera.image_size.x * cast(int)camera.image_size.y * 3 * 4)
	fmt.sbprint(&str, header)
	defer strings.builder_destroy(&str)

	render(&str, camera, spheres[:], true)

	fmt.fprint(os.stdout, strings.to_string(str))
}
