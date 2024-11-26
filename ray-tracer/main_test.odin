package main

import "core:os"
import "core:fmt"
import "core:strings"
import "core:math/rand"
import "core:testing"
import "core:time"

render_materials :: proc(str : ^strings.Builder, camera : camera, spheres : []sphere) {
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

	fmt.sbprintfln(str, "P3\n%v %v\n255", camera.image_size.x, camera.image_size.y)
	pixel_samples_scale := 1.0 / f64(camera.samples_per_pixel)
	for j in 0..<camera.image_size.y {
		// fmt.eprintf("\rScanlines remaining: %v ", camera.image_size.y - j)
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

	// fmt.eprintln("\rDone.                 ")
}

bench_materials ::
proc(options: ^time.Benchmark_Options, allocator := context.allocator) -> (err: time.Benchmark_Error) {

	rand.create(0)

	camera : camera
	camera_init(&camera, target_aspect_ratio=16.0/9.0, image_width=400, samples_per_pixel=100, max_ray_bounces=50)

	// Scene
	ground := material{procedure=lambertian_proc, data=&lambertian_data{albedo={0.8, 0.8, 0.0}}}
	blue   := material{procedure=lambertian_proc, data=&lambertian_data{albedo={0.1, 0.2, 0.5}}}
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
	defer {
		if fd, ferr := os.open("out/test.ppm", os.O_CREATE | os.O_RDWR); ferr != nil {
			err = .Allocation_Error
			// panic("couldn't open test.ppm")
		} else {
			fmt.fprintln(fd, strings.to_string(str))
			os.close(fd)
		}
		strings.builder_destroy(&str)
	}

	for _/* round */ in 1..=options.rounds {
		// fmt.eprintf("\rRound %v/%v", round, options.rounds)
		strings.builder_reset(&str)
		render_materials(&str, camera, spheres)
	}
	// fmt.eprintln("\rDone.       ")

	options.count     = options.rounds
	options.processed = options.rounds * options.bytes
	return err
}

@(test)
benchmark_materials :: proc(t: ^testing.T) {
	str: strings.Builder
	strings.builder_init(&str, context.allocator)
	defer {
		fmt.println(strings.to_string(str))
		strings.builder_destroy(&str)
	}

	{
		name := "materials"
		options := &time.Benchmark_Options{
			rounds   = 100,
			setup    = nil,
			bench    = bench_materials,
			teardown = nil,
		}
		time.benchmark(options, context.allocator)
		benchmark_print(&str, name, options)
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
