package main

import "core:os"
import "core:fmt"
import "core:strings"
import "core:math"
import "core:math/rand"
import "core:testing"
import "core:time"

render_unoptimized :: proc(camera : camera, spheres : []sphere) {
	// rasterization

	// Calculate the horizontal and vertical delta vectors from pixel to pixel.
	pixel_deltas : vec3 = {
		camera.viewport.x / f64(camera.image_size.x), // 3.556/1280=0.002778
		-camera.viewport.y / f64(camera.image_size.y), // -2/720=-0.002778
		0,
	}

	// Calculate the location of the top left pixel.
	/*
		NOTE(viktor): this is used to calculate the ray_direction which substracts the camera.position away again
		so we might as well leave it out altogether
	*/
	viewport_top_left := /* camera.position */ - vec3{camera.viewport.x*0.5, -camera.viewport.y*0.5, camera.focal_length}
	pixel00_center_in_3d := viewport_top_left + 0.5 * pixel_deltas

	fmt.printfln("P3\n%v %v\n255", camera.image_size.x, camera.image_size.y)
	pixel_samples_scale := 1.0 / f64(camera.samples_per_pixel)
	for j in 0..<camera.image_size.y {
		fmt.eprintf("\rScanlines remaining: %v ", camera.image_size.y - j)
		for i in 0..<camera.image_size.x {
			pixel_color : color
			for _ in 0..<camera.samples_per_pixel {
				// get ray

				offset := sample_square()
				pixel_sample := pixel00_center_in_3d + pixel_deltas * (vec3{f64(i), f64(j), 0} + offset)
				ray_direction := pixel_sample /* - camera.position */
				r := ray{camera.position, ray_direction}

				// ray_color

				// NOTE(viktor): only needs to be normalized for the background gradient
				r.direction = normalize(r.direction)
				sample_color := ray_color(&r) // NOTE(viktor): background color (gradient)
				closest_t := math.inf_f64(0)
				for &sphere in spheres {
					if record, ok := hit_sphere_ranged(&sphere.center, sphere.radius, &r, {0, closest_t}); ok {
						if record.t < closest_t {
							closest_t = record.t
							sample_color = 0.5 * color(record.normal + 1)
						}
					}
				}
				pixel_color += sample_color
			}
			write_color(os.stdout, pixel_samples_scale * pixel_color)
		}
	}

	fmt.eprintln("\rDone.                 ")
}

render_normalize_conditionally :: proc(camera : camera, spheres : []sphere) {
	// rasterization

	// Calculate the horizontal and vertical delta vectors from pixel to pixel.
	pixel_deltas : vec3 = {
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

	fmt.printfln("P3\n%v %v\n255", camera.image_size.x, camera.image_size.y)
	pixel_samples_scale := 1.0 / f64(camera.samples_per_pixel)
	for j in 0..<camera.image_size.y {
		fmt.eprintf("\rScanlines remaining: %v ", camera.image_size.y - j)
		for i in 0..<camera.image_size.x {
			pixel_color : color
			for _ in 0..<camera.samples_per_pixel {
				// get ray

				offset := sample_square()
				pixel_sample := pixel00_center_in_3d + pixel_deltas * (vec3{f64(i), f64(j), 0} + offset)
				ray_direction := pixel_sample /* - camera.position */
				r := ray{camera.position, ray_direction}

				// ray_color

				closest_t := math.inf_f64(0)
				rec : hit_record
				for &sphere in spheres {
					if record, ok := hit_sphere_ranged(&sphere.center, sphere.radius, &r, {0, closest_t}); ok {
						if record.t < closest_t {
							closest_t = record.t
							rec = record
						}
					}
				}

				sample_color : color
				if closest_t < math.inf_f64(0) {
					sample_color = 0.5 * color(rec.normal + 1)
				}
				else {
					// NOTE(viktor): only needs to be normalized for the background gradient
					r.direction = normalize(r.direction)
					sample_color = ray_color(&r) // NOTE(viktor): background color (gradient)
				}
				pixel_color += sample_color
			}
			write_color(os.stdout, pixel_samples_scale * pixel_color)
		}
	}

	fmt.eprintln("\rDone.                 ")
}

render_set_color_conditionally :: proc(camera : camera, spheres : []sphere) {
	// rasterization

	// Calculate the horizontal and vertical delta vectors from pixel to pixel.
	pixel_deltas : vec3 = {
		camera.viewport.x / f64(camera.image_size.x), // 3.556/1280=0.002778
		-camera.viewport.y / f64(camera.image_size.y), // -2/720=-0.002778
		0,
	}

	// Calculate the location of the top left pixel.
	/*
		NOTE(viktor): this is used to calculate the ray_direction which substracts the camera.position away again
		so we might as well leave it out altogether
	*/
	viewport_top_left := /* camera.position */ - vec3{camera.viewport.x*0.5, -camera.viewport.y*0.5, camera.focal_length}
	pixel00_center_in_3d := viewport_top_left + 0.5 * pixel_deltas

	fmt.printfln("P3\n%v %v\n255", camera.image_size.x, camera.image_size.y)
	pixel_samples_scale := 1.0 / f64(camera.samples_per_pixel)
	for j in 0..<camera.image_size.y {
		fmt.eprintf("\rScanlines remaining: %v ", camera.image_size.y - j)
		for i in 0..<camera.image_size.x {
			pixel_color : color
			for _ in 0..<camera.samples_per_pixel {
				// get ray

				offset := sample_square()
				pixel_sample := pixel00_center_in_3d + pixel_deltas * (vec3{f64(i), f64(j), 0} + offset)
				ray_direction := pixel_sample /* - camera.position */
				r := ray{camera.position, ray_direction}

				// ray_color

				// NOTE(viktor): only needs to be normalized for the background gradient
				r.direction = normalize(r.direction)
				sample_color := ray_color(&r) // NOTE(viktor): background color (gradient)
				closest_t := math.inf_f64(0)
				rec : hit_record
				for &sphere in spheres {
					if record, ok := hit_sphere_ranged(&sphere.center, sphere.radius, &r, {0, closest_t}); ok {
						if record.t < closest_t {
							closest_t = record.t
							rec = record
						}
					}
				}
				if closest_t < math.inf_f64(0) {
					sample_color = 0.5 * color(rec.normal + 1)
				}
				pixel_color += sample_color
			}
			write_color(os.stdout, pixel_samples_scale * pixel_color)
		}
	}

	fmt.eprintln("\rDone.                 ")
}

bench_render_normalize_conditionally ::
proc(options: ^time.Benchmark_Options, allocator := context.allocator) -> (err: time.Benchmark_Error) {

	rand.create(0)

	camera : camera
	camera_init(&camera, target_aspect_ratio=16.0/9.0, image_width=400, samples_per_pixel=100)

	// Scene
	spheres :: []sphere{
		{center={0, 0, -1}, radius=0.5},
		// {center={0.2, 0.5, -1}, radius=0.2},
		// {center={0.3, -0.1, -0.7}, radius=0.2},
		// {center={-0.5, -0.25, -0.5}, radius=0.2},
		// {center={0.65, 0.6, -0.8}, radius=0.3},
		{center={0, -100.5, -1}, radius=100},
	}

	for _/* round */ in 1..=options.rounds {
		// fmt.eprintf("\rRound %v/%v", round, options.rounds)
		render_normalize_conditionally(camera, spheres)
	}
	// fmt.eprintln("\rDone.       ")

	options.count     = options.rounds
	options.processed = options.rounds * options.bytes
	return nil
}

bench_render_set_color_conditionally ::
proc(options: ^time.Benchmark_Options, allocator := context.allocator) -> (err: time.Benchmark_Error) {

	rand.create(0)

	camera : camera
	camera_init(&camera, target_aspect_ratio=16.0/9.0, image_width=400, samples_per_pixel=100)

	// Scene
	spheres :: []sphere{
		{center={0, 0, -1}, radius=0.5},
		// {center={0.2, 0.5, -1}, radius=0.2},
		// {center={0.3, -0.1, -0.7}, radius=0.2},
		// {center={-0.5, -0.25, -0.5}, radius=0.2},
		// {center={0.65, 0.6, -0.8}, radius=0.3},
		{center={0, -100.5, -1}, radius=100},
	}

	for _/* round */ in 1..=options.rounds {
		// fmt.eprintf("\rRound %v/%v", round, options.rounds)
		render_set_color_conditionally(camera, spheres)
	}
	// fmt.eprintln("\rDone.       ")

	options.count     = options.rounds
	options.processed = options.rounds * options.bytes
	return nil
}

bench_render_unoptimized ::
proc(options: ^time.Benchmark_Options, allocator := context.allocator) -> (err: time.Benchmark_Error) {

	rand.create(0)

	camera : camera
	camera_init(&camera, target_aspect_ratio=16.0/9.0, image_width=400, samples_per_pixel=100)

	// Scene
	spheres :: []sphere{
		{center={0, 0, -1}, radius=0.5},
		// {center={0.2, 0.5, -1}, radius=0.2},
		// {center={0.3, -0.1, -0.7}, radius=0.2},
		// {center={-0.5, -0.25, -0.5}, radius=0.2},
		// {center={0.65, 0.6, -0.8}, radius=0.3},
		{center={0, -100.5, -1}, radius=100},
	}

	for _/* round */ in 1..=options.rounds {
		// fmt.eprintf("\rRound %v/%v", round, options.rounds)
		render_unoptimized(camera, spheres)
	}
	// fmt.eprintln("\rDone.       ")

	options.count     = options.rounds
	options.processed = options.rounds * options.bytes
	return nil
}

@(test)
benchmark_render :: proc(t: ^testing.T)
{
	str: strings.Builder
	strings.builder_init(&str, context.allocator)
	defer {
		fmt.println(strings.to_string(str))
		strings.builder_destroy(&str)
	}

	{
		name := "render_normalize_conditionally"
		options := &time.Benchmark_Options{
			rounds   = 100,
			setup    = nil,
			bench    = bench_render_normalize_conditionally,
			teardown = nil,
		}
		time.benchmark(options, context.allocator)
		benchmark_print(&str, name, options)
		// 100 rounds in 29370445300 ns
		// 3.405 rounds/s, 0.000 MiB/s
	}

	{
		name := "render_set_color_conditionally"
		options := &time.Benchmark_Options{
			rounds   = 100,
			setup    = nil,
			bench    = bench_render_set_color_conditionally,
			teardown = nil,
		}
		time.benchmark(options, context.allocator)
		benchmark_print(&str, name, options)
		// 100 rounds in 36314637200 ns
		// 2.754 rounds/s, 0.000 MiB/s
	}

	{
		name := "render_unoptimized"
		options := &time.Benchmark_Options{
			rounds   = 100,
			setup    = nil,
			bench    = bench_render_unoptimized,
			teardown = nil,
		}
		time.benchmark(options, context.allocator)
		benchmark_print(&str, name, options)
		// 100 rounds, 0 bytes processed in 66535458600 ns
		// 1.503 rounds/s, 0.000 MiB/s
	}

}

benchmark_print :: proc(str: ^strings.Builder, name: string, options: ^time.Benchmark_Options, loc := #caller_location) {
	fmt.sbprintfln(str, "[%v] %v rounds, %v bytes processed in %v ns\n\t\t%5.3f rounds/s, %5.3f MiB/s\n",
		name,
		options.rounds,
		options.processed,
		time.duration_nanoseconds(options.duration),
		options.rounds_per_second,
		options.megabytes_per_second,
	)
}
