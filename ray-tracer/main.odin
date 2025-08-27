package main

import "core:fmt"
import "core:math"
import "core:math/rand"
import "core:strings"


/// types

v3    :: [3]f64
v2    :: [2]f64
turns :: distinct f64 // [0, 1]
// radians :: distinct f64 // [0, math.TAU]
turns_to_radians :: proc "contextless" (θ: turns) -> f64 { return cast(f64)θ * math.TAU }


/// linalg

dot :: proc "contextless" (a, b: v3) -> f64 {
	return a.x*b.x + a.y*b.y + a.z*b.z
}

cross :: proc "contextless" (a, b: v3) -> v3 {
	result: v3
	result[0] = a[1]*b[2] - b[1]*a[2]
	result[1] = a[2]*b[0] - b[2]*a[0]
	result[2] = a[0]*b[1] - b[0]*a[1]
	return result
}

magnitude_squared :: proc "contextless" (vec: v3) -> f64 {
	return dot(vec, vec)
}

magnitude :: proc "contextless" (vec: v3) -> f64 {
	return math.sqrt(magnitude_squared(vec))
}

normalize :: proc "contextless" (vec: v3) -> v3 {
	return vec / magnitude(vec)
}

is_near_zero :: proc "contextless" (vec: v3) -> bool {
	EPSILON :: 1e-8
	return (math.abs(vec.x) <= EPSILON &&
	        math.abs(vec.y) <= EPSILON &&
	        math.abs(vec.z) <= EPSILON)
}

is_normalized :: proc "contextless" (vec: v3) -> bool {
	EPSILON :: 1e-10
	mag2 := magnitude_squared(vec)
	return 1-EPSILON <= mag2 && mag2 <= 1+EPSILON
}

lookat :: proc "contextless" (position: v3 = {}, target: v3 = {0, 0, -1}, axis_up: v3 = {0, 1, 0}) -> (right, up, forward: v3) {
	// view cartesian coordinate system
	forward = normalize(target - position)
	right   = normalize(cross(forward, axis_up))
	up      = cross(right, forward)
	return
}


/// ray tracing

// ![reflection](https://raytracing.github.io/images/fig-1.15-reflection.jpg|width=200)
reflect :: proc "contextless" (vec, normal: v3) -> v3 {
	// vec & normal don't need to be normalized
	return vec-2*dot(vec, normal)*normal
}

// cite: Ray Tracing Gems II (pp.109-114): Majercik, Zander. (2021). The Schlick Fresnel Approximation. 10.1007/978-1-4842-7185-8_9.
// ref: https://www.researchgate.net/publication/354065225_The_Schlick_Fresnel_Approximation
reflectance_fresnel :: proc "contextless" (cos_i, sin_i, src_refractive_index, dst_refractive_index: f64) -> f64 {
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

reflectance_at_normal_incidence :: proc "contextless" (rel_refractive_index: f64) -> f64 {
	r0 := (1.0 - rel_refractive_index) / (1.0 + rel_refractive_index)
	r0 *= r0
	return r0
}

// cite: Schlick, C. An inexpensive BRDF model for physically-based rendering.
// ref: https://onlinelibrary.wiley.com/doi/10.1111/1467-8659.1330233
reflectance_schlick_approximation :: proc "contextless" (cos_i, r0: f64) -> f64 {
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
reflectance_schlick_lazanyi_approximation :: proc "contextless" (cos_i, r0: f64, a: f64, alpha: f64) -> f64 {
	return reflectance_schlick_approximation(r0, cos_i) - a * cos_i * math.pow(1 - cos_i, alpha)
}
// reflectance_schlick_lazanyi_approximation :: proc(r0: v3, cos_i: float, a: v3, alpha: float) -> f64 {
// 	return schlickFresnel(r0, cos_i) - a * cos_i * math.pow(1 - cos_i, alpha);
// }

// cite: Hoffman, N. Fresnel equations considered harmful. In Eurographics Workshop on Material Appearance Modeling, pages 7–11, 2019. DOI: 10.2312/mam.20191305.
// ref: https://diglib.eg.org/server/api/core/bitstreams/726dc384-d7dd-4c0e-8806-eadec0ff3886/content
reflectance_hoffman_approximation :: proc "contextless" (cos_i, r0, h: f64) -> f64 {
	a := 823543 / 46656 * (r0 - h) + 49 / 6 * (1 - r0)
	return reflectance_schlick_lazanyi_approximation(cos_i, r0, a, alpha = 6)
}

// * src_refractive_index: refractive index of the medium the light is exiting, aka. incident refractive index
// * dst_refractive_index: refractive index of the medium the light is entering, aka. transmitted refractive index
refract_with_reference_medium :: proc(vec, normal: v3, src_refractive_index, dst_refractive_index: f64) -> v3 {
	return refract_with_relative_refractive_index(vec, normal, src_refractive_index / dst_refractive_index)
}

// ![refraction](https://raytracing.github.io/images/fig-1.17-refraction.jpg|width=200)
refract_with_relative_refractive_index :: proc(vec, normal: v3, relative_refractive_index: f64) -> v3 {
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

random_v3 :: proc(gen := context.random_generator) -> v3 {
	return {rand.float64(gen), rand.float64(gen), rand.float64(gen)}
}

random_v3_range :: proc(min, max : f64, gen := context.random_generator) -> v3 {
	return {rand.float64_range(min, max, gen), rand.float64_range(min, max, gen), rand.float64_range(min, max, gen)}
}

// Returns the vector to a random point in the [min,min]-[max,max] square.
random_v2_range :: proc(min, max : f64, gen := context.random_generator) -> v3 {
	return {rand.float64_range(min, max, gen), rand.float64_range(min, max, gen), 0}
}

random_unit_vector :: proc(gen := context.random_generator) -> v3 {
	// rejection sampling to ensure normal distribution
	for {
		p := random_v3_range(-1, 1 + math.F64_EPSILON, gen) // excl. max
		length_squared := magnitude_squared(p)
		if 0 < length_squared && length_squared <= 1 {
			return p / math.sqrt(length_squared)
		}
	}
}

random_point_on_hemisphere :: proc(normal: v3, gen := context.random_generator) -> v3 {
	on_unit_sphere := random_unit_vector()
	if dot(on_unit_sphere, normal) > 0.0 { // In the same hemisphere as the normal
		return on_unit_sphere
	} else {
		return -on_unit_sphere
	}
}

random_point_on_disk :: proc(gen := context.random_generator) -> v3 {
	// rejection sampling to ensure normal distribution
	for {
		p := random_v2_range(-1, 1 + math.F64_EPSILON) // excl. max
		length_squared := magnitude_squared(p)
		if length_squared <= 1 {
			return p
		}
	}
}


/// raycast

ray :: struct {
	origin : v3,
	direction : v3,
}

hit_record :: struct {
	p : v3,
	normal : v3,
	front_face : bool,
}

// returns value in range [t_min, t_max) or t_max if no hit
ray_sphere_intersection :: #force_inline proc "contextless" (
	r: ray, sphere_center: v3, sphere_radius: f64, t_min: f64 = 0, t_max: f64 = math.F64_MAX) -> (t: f64) {

	t = t_max

	o_c := sphere_center - r.origin
	a := magnitude_squared(r.direction)
	h := dot(r.direction, o_c)
	c := magnitude_squared(o_c) - sphere_radius*sphere_radius
	discriminant := h*h - a*c
	if discriminant < 0 do return

	sqrtd := math.sqrt(discriminant)

	// Find the nearest root that lies in the acceptable range
	root := (h - sqrtd) / a
	if root < t_min || t_max <= root {
		root = (h + sqrtd) / a
		if root < t_min || t_max <= root do return
	}

	t = root
	return
}


/// materials

material_type :: enum i64 {
	lambertian,
	metallic,
	dielectric,
}

material_data :: struct {
	albedo: v3,
	param1: f64, // fuzz for metal, refractive_index for dielectric
}

make_metallic_data :: proc(albedo: v3, fuzz: f64) -> material_data {
	assert(0 <= fuzz && fuzz <= 1.0)
	// data^ = {albedo, math.clamp(fuzz, 0.0, 1.0)}
	return {albedo=albedo, param1=fuzz}
}

material :: struct {
	// TODO: Maybe get rid of the type at this place? it isn't used inside the actual scatter procs, so its just bloating up the materials array, it's only relevant info at branch-time
	type: material_type,
	data: material_data,
}

// ![](https://raytracing.github.io/images/fig-1.14-rand-unitvec.jpg|width=200)
lambertian_proc :: proc(data: material_data, ray_in: ^ray, hit: ^hit_record) ->
                       (ray_out: ray, attenuation: v3, ok: bool) {
	// Lambertian (diffuse) reflectance can either
	//  * always scatter and attenuate light according to its reflectance R,
	//  * or it can sometimes scatter (with probability 1−R) with no attenuation (where a ray that isn't scattered is just absorbed into the material).
	// It could also be a mixture of both those strategies. We will choose to always scatter

	output_direction := hit.normal + random_unit_vector()
	if is_near_zero(output_direction) {
		output_direction = hit.normal
	}

	ray_out = ray{hit.p, output_direction}

	attenuation = data.albedo
	// Note the third option: we could scatter with some fixed probability p and have attenuation be albedo/p.
	// p := fixed probability
	// attenuation = albedo/p
	ok = true
	return
}

// ![](https://raytracing.github.io/images/fig-1.16-reflect-fuzzy.jpg|width=200)
metallic_proc :: proc(data: material_data, ray_in: ^ray, hit: ^hit_record) ->
					 (ray_out: ray, attenuation: v3, ok: bool) {
	albedo := data.albedo
	fuzz := data.param1

	output_direction := normalize(reflect(ray_in.direction, hit.normal))
	output_direction += fuzz * random_unit_vector()
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
			r0 := reflectance_at_normal_incidence(1.0/HARDCODED_REFRACTIVE_INDEX)
			reflection_factor := reflectance_schlick_approximation(cos_theta, r0)
		} else when METAL_FRESNEL_KIND == 1 {
			H :: 0.5
			r0 := reflectance_at_normal_incidence(1.0/HARDCODED_REFRACTIVE_INDEX)
			reflection_factor := reflectance_hoffman_approximation(cos_theta, r0, H)
		} else {
			sin_theta_squared := 1.0-cos_theta*cos_theta
			sin_theta := math.sqrt(sin_theta_squared)
			reflection_factor := reflectance_fresnel(cos_theta, sin_theta, 1.0, HARDCODED_REFRACTIVE_INDEX)
		}

		// refraction_color*refraction_factor + reflection_color*reflection_factor
		attenuation = math.lerp(albedo, {1.0, 1.0, 1.0}, reflection_factor)
	} else {
		attenuation = albedo
	}
	ok = dot(output_direction, hit.normal) > 0
	return
}

// ![](https://raytracing.github.io/images/fig-1.17-refraction.jpg|width=200)
dielectric_proc :: proc(data: material_data, ray_in: ^ray, hit: ^hit_record) -> (ray_out: ray, attenuation: v3, ok: bool) {
	refractive_index := data.param1

	// materials with a refractive_index < 1 (air/vaccuum) are treated as air/vaccuum materials
	// inside a denser material with the effective refractive index of 1/refractive_index
	src_refractive_index := 1.0              if refractive_index >= 1.0 else 1.0 / refractive_index
	dst_refractive_index := refractive_index if refractive_index >= 1.0 else 1.0

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
	output_direction: v3

	if must_reflect || reflectance_schlick_approximation(cos_theta, reflectance_at_normal_incidence(rel_refractive_index)) > rand.float64() {
		output_direction = reflect(unit_direction, hit.normal)
	} else {
		output_direction = refract(unit_direction, hit.normal, rel_refractive_index)
	}
	ray_out = ray{hit.p, output_direction}

	attenuation = {1,1,1}
	ok = true
	return
}


material_scatter :: proc(material: material, ray_in: ^ray, hit: ^hit_record) -> (ray_out: ray, attenuation: v3, ok: bool) {
	switch material.type {
	case .lambertian: return lambertian_proc(material.data, ray_in, hit)
	case .metallic:   return   metallic_proc(material.data, ray_in, hit)
	case .dielectric: return dielectric_proc(material.data, ray_in, hit)
	}
	unreachable()
}

background_color :: proc "contextless" (r: ^ray) -> v3 {
	// linear gradient between a and b
	a := v3{1.0, 1.0, 1.0}
	b := v3{0.5, 0.7, 1.0}
	t := 0.5 * (r.direction.y + 1.0)
	return math.lerp(a, b, t)
}

image :: struct {
	width, height: i64,
	data: []v3,
}

camera :: struct {
	position : v3,
	right: v3,
	up: v3,
	forward: v3,
	// orientation : quaternion256,
	aspect_ratio : f64,
	image_size : v2,
	focus_distance: f64,
	vfov : turns, // using Hor+ scaling
	depth_of_field_angle: turns,
	samples_per_pixel : i64,
	max_ray_bounces : i64,
}

sphere :: struct {
	center : v3,
	radius : f64,
	material_index: i64,
}

world :: struct {
	materials: [dynamic]material,
	spheres: [dynamic]sphere,
}

world_destroy :: proc(world: ^world) {
	delete(world.spheres)
	delete(world.materials)
}

render :: proc(image: image, camera: camera, world: ^world, $print_progress : bool) {
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

	// TODO(viktor): test behavior when sensor size != render target resolution
	image_width  := cast(int)image.width
	image_height := cast(int)image.height
	pixel_sample_contribution_factor := 1.0 / cast(f64)camera.samples_per_pixel
	for v in 0..<image_height {
		when print_progress do fmt.eprintf("\rScanlines remaining: %v ", image_height - v)
		for u in 0..<image_width {
			pixel_color: v3
			for _ in 0..<camera.samples_per_pixel {
				// if we include max we may sample the max borders twice with adjacent pixels
				offset := random_v2_range(-0.5, 0.5/* +math.F64_EPSILON */)
				pixel_sample_position :=
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
				ray_direction := pixel_sample_position - ray_origin
				ray := ray{ray_origin, ray_direction}
				sample_color := v3{1, 1, 1}
				ray_cast: for _ in 0..=camera.max_ray_bounces {
					spheres := world.spheres
					closest_t := math.F64_MAX
					hit_sphere_index := -1
					SHADOW_ACNE_RAY_OFFSET :: 0.001
					for i in 0..<len(spheres) {
						t := ray_sphere_intersection(ray, spheres[i].center, spheres[i].radius, SHADOW_ACNE_RAY_OFFSET, closest_t)
						if t < closest_t {
							closest_t = t
							hit_sphere_index = i
						}
					}

					if hit_sphere_index >= 0 {
						sphere := &spheres[hit_sphere_index]
						material := &world.materials[sphere.material_index]

						// TODO: maybe move this all the way into the concrete scatter function so this stuff is only calculated if actually needed
						closest_hit: hit_record
						closest_hit.p = ray.origin + closest_t*ray.direction
						closest_hit.normal = (closest_hit.p - sphere.center) / sphere.radius
						closest_hit.front_face = dot(ray.direction, closest_hit.normal) < 0.0
						closest_hit.normal = closest_hit.front_face ? closest_hit.normal : -closest_hit.normal

						// TODO(viktor): @perf: return info whether the ray is inside a sphere and check that sphere first in the next iteration (a bvh would also achieve this)
						// TODO(viktor): instead of including material types in the materials array encode the info in sphere.material_index via bit-shifting
						if scattered_ray, attenuation, ok := material_scatter(material^, &ray, &closest_hit); ok {
							ray = scattered_ray
							sample_color *= attenuation
						} else {
							sample_color = {0,0,0}
							break ray_cast
						}
					} else {
						// NOTE(viktor): only needs to be normalized for the background gradient (dielectric_proc also does it though)
						ray.direction = normalize(ray.direction)
						sample_color *= background_color(&ray)
						break ray_cast
					}
				}
				pixel_color += sample_color
			}
			pixel_color *= pixel_sample_contribution_factor

			// TODO: move this out in front of serialization as a filter pass?
			// linear to gamma2 color correction
			pixel_color.r = math.sqrt(pixel_color.r)
			pixel_color.g = math.sqrt(pixel_color.g)
			pixel_color.b = math.sqrt(pixel_color.b)

			index := u + v * image_width
			image.data[index] = pixel_color
		}
	}

	when print_progress do fmt.eprintln("\rDone.                   ")
}

build_dev_scene :: proc(allocator := context.allocator) -> (camera: camera, world: world) {
	materials := make(type_of(world.materials), allocator)
	/* ground */     append(&materials, material{.lambertian, {albedo={0.8, 0.8, 0.0}}})
	/* blue */       append(&materials, material{.lambertian, {albedo={0.1, 0.2, 0.5}}})
	/* glass */      append(&materials, material{.dielectric, {param1=1.5}})
	/* air_bubble */ append(&materials, material{.dielectric, {param1=1.0/1.5}})
	/* gold */       append(&materials, material{.metallic,   make_metallic_data(albedo={0.8, 0.6, 0.2}, fuzz=1.0)})
	// /* silver */     append(&materials, material{.metallic,   make_metallic_data(albedo={0.8, 0.8, 0.8}, fuzz=0.3)})

	// TODO: an associative relation between spheres & materials would be nicer to ensure that changes to the materials order don't affect sphere representation
	spheres := make(type_of(world.spheres), allocator)
	append(&spheres, sphere{center={ 0.0, -100.5, -1.0}, radius=100, material_index=0/* ground */})
	append(&spheres, sphere{center={ 0.0,    0.0, -1.2}, radius=0.5, material_index=1/* blue */})
	append(&spheres, sphere{center={-1.0,    0.0, -1.0}, radius=0.5, material_index=2/* glass */})
	append(&spheres, sphere{center={-1.0,    0.0, -1.0}, radius=0.4, material_index=3/* air_bubble */})
	append(&spheres, sphere{center={ 1.0,    0.0, -1.0}, radius=0.5, material_index=4/* gold */})

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

	return camera, {materials, spheres}
}

build_final_scene :: proc(allocator := context.allocator) -> (camera: camera, world: world) {
	materials := make(type_of(world.materials), allocator)
	append(&materials, material{.lambertian, {albedo={0.5, 0.5, 0.5}}})
	append(&materials, material{.dielectric, {param1=1.5}})
	append(&materials, material{.lambertian, {albedo={0.4, 0.2, 0.1}}})
	append(&materials, material{.metallic,   {albedo={0.7, 0.6, 0.5}, param1=0.0}})

	spheres := make(type_of(world.spheres), allocator)
	append(&spheres, sphere{center={ 0.0, -1000, 0}, radius=1000, material_index=0})
	append(&spheres, sphere{center={   0,     1, 0}, radius= 1.0, material_index=1})
	append(&spheres, sphere{center={  -4,     1, 0}, radius= 1.0, material_index=2})
	append(&spheres, sphere{center={   4,     1, 0}, radius= 1.0, material_index=3})

	for a in -11..<11 {
		for b in -11..<11 {
			center := v3{cast(f64)a + 0.9 * rand.float64(), 0.2, cast(f64)b + 0.9 * rand.float64()}

			if magnitude(center - v3{4, 0.2, 0}) > 0.9 {
				choose_mat := rand.float64()
				if choose_mat < 0.8 {
					// diffuse
					albedo := random_v3()*random_v3()
					append(&materials, material{.lambertian, {albedo=albedo}})
					append(&spheres, sphere{center=center, radius=0.2, material_index=cast(i64)len(materials)-1})
				} else if choose_mat < 0.95 {
					// metal
					albedo := random_v3_range(0.5, 1)
					fuzz := rand.float64_range(0, 0.5)
					append(&materials, material{.metallic, make_metallic_data(albedo, fuzz)})
					append(&spheres, sphere{center=center, radius=0.2, material_index=cast(i64)len(materials)-1})
				} else {
					append(&materials, material{.dielectric, {param1=1.5}})
					append(&spheres, sphere{center=center, radius=0.2, material_index=cast(i64)len(materials)-1})
				}
			}
		}
	}

	camera.position = {13, 2, 3}
	camera.right, camera.up, camera.forward = lookat(position=camera.position, target={0, 0, 0}, axis_up={0, 1, 0})
	camera.aspect_ratio = 16.0/9.0
	camera.image_size.x = 1200
	camera.image_size.y = camera.image_size.x/camera.aspect_ratio
	camera.focus_distance = 10
	camera.vfov = 20.0/360.0
	camera.depth_of_field_angle = 0.6/360.0
	camera.samples_per_pixel = 500
	camera.max_ray_bounces = 50

	return camera, {materials, spheres}
}

serialize_ppm :: proc(str: ^strings.Builder, image: image) -> string {
	fmt.sbprintfln(str, "P3\n%v %v\n255", image.width, image.height)
	header_size := cast(i64)len(str.buf)
	CHARS_PER_CHANNEL :: 4
	non_zero_resize_dynamic_array(&str.buf, header_size + image.width * image.height * CHARS_PER_CHANNEL * len(v3))

	serialize_channel :: proc(dst: []byte, u: u8, separator: byte) {
		dst[0] = '0' + ((u / 100) % 10)
		dst[1] = '0' + ((u /  10) % 10)
		dst[2] = '0' + ( u        % 10)
		dst[3] = separator

		// replace leading zeroes
		dst[0] = dst[0]=='0'                ? ' ' : dst[0]
		dst[1] = dst[0]==' ' && dst[1]=='0' ? ' ' : dst[1]
	}

	output_stride := CHARS_PER_CHANNEL * len(v3)
	output_slice := str.buf[header_size:]
	for pixel, i in image.data {
		// quantize [0.0,1.0] float values to [0,255] byte range.
		r := u8(256 * math.min(pixel.r, 0.999))
		g := u8(256 * math.min(pixel.g, 0.999))
		b := u8(256 * math.min(pixel.b, 0.999))

		offset := i * output_stride
		output_pixel := output_slice[offset:offset+output_stride]
		serialize_channel(output_pixel[0: 4], r,  ' ')
		serialize_channel(output_pixel[4: 8], g,  ' ')
		serialize_channel(output_pixel[8:12], b, '\n')
	}

	return strings.to_string(str^)
}

main :: proc () {
	rand.reset(1)

	camera, world := build_dev_scene(context.allocator)
	defer world_destroy(&world)

	image: image
	image.width = cast(i64)camera.image_size.x
	image.height = cast(i64)camera.image_size.y
	image.data   = make([]v3, image.width * image.height, context.allocator)
	defer delete(image.data, context.allocator)

	render(image, camera, &world, true)

	str: strings.Builder
	strings.builder_init(&str, context.allocator)
	defer strings.builder_destroy(&str)
	serialized := serialize_ppm(&str, image)

	fmt.print(serialized)
}
