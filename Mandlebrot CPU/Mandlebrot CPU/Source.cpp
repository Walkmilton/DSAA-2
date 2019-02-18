// Mandelbrot set example
// Adam Sampson <a.sampson@abertay.ac.uk>


/*
2126, 3969
2113, 3933
2124, 3933
*/

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <complex>
#include <fstream>
#include <iostream>
#include <thread>

// Import things we need from the standard library
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::complex;
using std::cout;
using std::endl;
using std::ofstream;
using std::thread;

// Define the alias "the_clock" for the clock type we're going to use.
typedef std::chrono::steady_clock the_clock;


// The size of the image to generate.
const int WIDTH = 640;
const int HEIGHT = 480;

// The number of times to iterate before we assume that a point isn't in the
// Mandelbrot set.
// (You may need to turn this up if you zoom further into the set.)
const int MAX_ITERATIONS = 500;

// The image data.
// Each pixel is represented as 0xRRGGBB.
uint32_t image[HEIGHT][WIDTH];


// Write the image to a TGA file with the given name.
// Format specification: http://www.gamers.org/dEngine/quake3/TGA.txt
void write_tga(const char *filename)
{
	ofstream outfile(filename, ofstream::binary);

	uint8_t header[18] = {
		0, // no image ID
		0, // no colour map
		2, // uncompressed 24-bit image
		0, 0, 0, 0, 0, // empty colour map specification
		0, 0, // X origin
		0, 0, // Y origin
		WIDTH & 0xFF, (WIDTH >> 8) & 0xFF, // width
		HEIGHT & 0xFF, (HEIGHT >> 8) & 0xFF, // height
		24, // bits per pixel
		0, // image descriptor
	};
	outfile.write((const char *)header, 18);

	for (int y = 0; y < HEIGHT; ++y)
	{
		for (int x = 0; x < WIDTH; ++x)
		{
			uint8_t pixel[3] = {
				image[y][x] & 0xFF, // blue channel
				(image[y][x] >> 8) & 0xFF, // green channel
				(image[y][x] >> 16) & 0xFF, // red channel
			};
			outfile.write((const char *)pixel, 3);
		}
	}

	outfile.close();
	if (!outfile)
	{
		// An error has occurred at some point since we opened the file.
		cout << "Error writing to " << filename << endl;
		exit(1);
	}
}


// Render the Mandelbrot set into the image array.
// The parameters specify the region on the complex plane to plot.
void compute_mandelbrot(double left, double right, double top, double bottom, double YStart, double YEnd)
{
	for (int y = YStart; y < YEnd; ++y)
	{
		for (int x = 0; x < WIDTH; ++x)
		{
			// Work out the point in the complex plane that
			// corresponds to this pixel in the output image.
			complex<double> c(left + (x * (right - left) / WIDTH),
				top + (y * (bottom - top) / HEIGHT));

			// Start off z at (0, 0).
			complex<double> z(0.0, 0.0);

			// Iterate z = z^2 + c until z moves more than 2 units
			// away from (0, 0), or we've iterated too many times.
			int iterations = 0;
			while (abs(z) < 2.0 && iterations < MAX_ITERATIONS)
			{
				z = (z * z) + c;

				++iterations;
			}

			if (iterations < 16) {
				image[y][x] = 0xFF0000;
			}

			if (iterations < 8) {
				image[y][x] = 0x00FF00;
			}

			if (iterations == MAX_ITERATIONS)
			{
				// z didn't escape from the circle.
				// This point is in the Mandelbrot set.
				image[y][x] = MAX_ITERATIONS; // black
			}
		}
	}
}


int main(int argc, char *argv[])
{
	cout << "Please wait..." << endl;

	const int threads = 8;
	int z = 480 / threads;

	thread compute_thread[threads];



	// This shows the whole set.
	int j = 0;
	int x = 0;
	for (int i = z; i <= 480, j <= (480 - z); i += z, j += z) {

		the_clock::time_point start = the_clock::now();

		//compute_mandelbrot(-2.0, 1.0, 1.125, -1.125, j, i);

		compute_thread[x] = thread(compute_mandelbrot, -2.0, 1.0, 1.125, -1.125, j, i);


		the_clock::time_point end = the_clock::now();

		auto time_taken = duration_cast<milliseconds>(end - start).count();
		cout << "Computing the Mandelbrot set took " << time_taken << " ms." << endl;
		x++;
	}

	for (int i = 0; i < threads; i++)
	{
		compute_thread[i].join();
	}

	write_tga("output.tga");

	return 0;
}
