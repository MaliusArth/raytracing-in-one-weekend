package main

import "core:fmt"
import "core:strings"
import "core:testing"
import "core:time"

bench_materials ::
proc(options: ^time.Benchmark_Options, allocator := context.allocator) -> (err: time.Benchmark_Error) {
	// fmt.eprintfln("Rounds: %v", options.rounds)
	// log.infof("Rounds: %v", options.rounds)

	camera, materials, spheres := build_dev_scene(allocator)
	defer delete(materials)
	defer delete(spheres)

	image: image
	image.fourcc = "PPM3"
	image.width = cast(i64)camera.image_size.x
	image.height = cast(i64)camera.image_size.y
	image.data   = make([]color, image.width * image.height, allocator)
	defer delete(image.data, allocator)

	options.bytes = len(image.data)
	for round in 1..=options.rounds {
		fmt.eprintf("\rRound %v/%v", round, options.rounds)
			render(image, camera, materials[:], spheres[:], false)
	}
	fmt.eprintln("\rDone.        ")

	str: strings.Builder
	strings.builder_init(&str)
	defer strings.builder_destroy(&str)

	serialized := serialize(&str, image)

	fmt.print(serialized)

	options.count     = options.rounds
	options.processed = options.rounds * options.bytes
	return err
}

@(test)
benchmark_materials :: proc(t: ^testing.T) {
	str: strings.Builder
	strings.builder_init(&str)
	defer {
		fmt.eprintln(strings.to_string(str))
		// log.info(strings.to_string(str))
		strings.builder_destroy(&str)
	}

	{
		name := "materials"
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
