package main

import "core:os"
import "core:fmt"
import "core:strings"
import "core:math/rand"
import "core:testing"
import "core:time"

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

	str: strings.Builder
	header := fmt.aprintfln("P3\n%v %v\n255", cast(int)camera.image_size.x, cast(int)camera.image_size.y)
	defer delete_string(header)
	strings.builder_init(&str, 0, len(header) + cast(int)camera.image_size.x * cast(int)camera.image_size.y * 3 * 4)
	defer {
		if fd, ferr := os.open("out/test.ppm", os.O_CREATE | os.O_RDWR); ferr != nil {
			err = .Allocation_Error
			// panic("couldn't open test.ppm")
		} else {
			fmt.fprint(fd, strings.to_string(str))
			os.close(fd)
		}
		strings.builder_destroy(&str)
	}

	for _/* round */ in 1..=options.rounds {
		// fmt.eprintf("\rRound %v/%v", round, options.rounds)
		strings.builder_reset(&str)
		fmt.sbprint(&str, header)
		render(&str, camera, spheres, false)
	}
	// fmt.eprintln("\rDone.       ")

	options.count     = options.rounds
	options.processed = options.rounds * options.bytes
	return err
}

@(test)
benchmark_materials_once :: proc(t: ^testing.T) {
	str: strings.Builder
	strings.builder_init(&str)
	defer {
		fmt.println(strings.to_string(str))
		strings.builder_destroy(&str)
	}

	{
		name := "materials once"
		options := &time.Benchmark_Options{
			rounds   = 1,
			setup    = nil,
			bench    = bench_materials,
			teardown = nil,
		}
		time.benchmark(options)
		benchmark_print(&str, name, options)
	}

}

@(test)
benchmark_materials :: proc(t: ^testing.T) {
	str: strings.Builder
	strings.builder_init(&str)
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
		time.benchmark(options)
		benchmark_print(&str, name, options)
	}

}

benchmark_print :: proc(str: ^strings.Builder, name: string, options: ^time.Benchmark_Options, loc := #caller_location) {
	fmt.sbprintfln(str, "[%v] %v rounds, %v bytes processed in %v ns, %5.3f rounds/s, %5.3f MiB/s\n",
		name,
		options.rounds,
		options.processed,
		time.duration_nanoseconds(options.duration),
		options.rounds_per_second,
		options.megabytes_per_second,
	)
}
