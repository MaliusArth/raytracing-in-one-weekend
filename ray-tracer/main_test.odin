package main

// import "core:log"
import "core:fmt"
import "core:strings"
// import "core:math/rand"
import "core:testing"
import "core:time"
// import "core:slice"

// build_test_scene :: proc(allocator := context.allocator) -> (camera, [dynamic]sphere) {
// 	spheres := make_dynamic_array([dynamic]sphere, allocator)
// 	ground := material_make(color{0.8, 0.8, 0.0}, allocator)
// 	blue   := material_make(color{0.1, 0.2, 0.5}, allocator)
// 	silver := material_make(color{0.8, 0.8, 0.8}, allocator)
// 	gold   := material_make(color{0.8, 0.6, 0.2}, allocator)
// 	append(&spheres, sphere{center={ 0.0, -100.5, -1.0}, radius=100, material=ground})
// 	append(&spheres, sphere{center={ 0.0,    0.0, -1.2}, radius=0.5, material=blue})
// 	append(&spheres, sphere{center={-1.0,    0.0, -1.0}, radius=0.5, material=silver})
// 	append(&spheres, sphere{center={ 1.0,    0.0, -1.0}, radius=0.5, material=gold})

// 	camera: camera
// 	camera_init(&camera)

// 	return camera, spheres
// }

bench_materials ::
proc(options: ^time.Benchmark_Options, allocator := context.allocator) -> (err: time.Benchmark_Error) {
	// fmt.eprintfln("Rounds: %v", options.rounds)
	// log.infof("Rounds: %v", options.rounds)
	// rand.reset(seed=0)

	camera, spheres := build_dev_scene(allocator)
	defer for sphere in spheres { free(sphere.material.data, allocator) }
	defer delete(spheres)

	image: image
	image.fourcc = "PPM3"
	image.resolution.x = cast(i64)camera.image_size.x
	image.resolution.y = cast(i64)camera.image_size.y
	image.data   = make([]color, image.resolution.x * image.resolution.y, allocator)
	defer delete(image.data, allocator)

	options.bytes = len(image.data)
	for round in 1..=options.rounds {
		fmt.eprintf("\rRound %v/%v", round, options.rounds)
		render(image, camera, spheres[:], false)
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
