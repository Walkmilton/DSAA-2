//DSAA 2 Coursework
//Declan Doyle
//Code taken and modified from Adam Sampson, with help from lab worksheets. (Labs 2 & 8)

#include <iostream>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <complex>
#include <fstream>
#include <thread>
#include <amp.h>
#include <amp_math.h>

using std::cout;
using std::cin;
using std::endl;
using std::chrono::duration_cast;
using std::chrono::microseconds;
using std::complex;
using std::ofstream;
using std::thread;

using namespace concurrency;

//Defining clock
typedef std::chrono::steady_clock the_clock;


const int WIDTH = 2560;
const int HEIGHT = 1920;
const int MAX_ITERATIONS = 500;

uint32_t image[HEIGHT][WIDTH];


struct Complex1
{
	float x;
	float y;
};


void menuChoice(int menu);
void CPU();
void compute_mandelbrot_CPU(double left, double right, double top, double bottom, double YStart, double YEnd);

void GPU();
Complex1 c_add(Complex1 c1, Complex1 c2) restrict(cpu, amp);
float c_abs(Complex1 c) restrict(cpu, amp);
Complex1 c_mul(Complex1 c1, Complex1 c2) restrict(cpu, amp);
void compute_mandelbrot_GPU(double left, double right, double top, double bottom);


void write_tga(const char *filename);

int main()
{


	//menu variiable
  int menu = 0;

  //Displaying the menu and inputting a menu choice
	while (menu != -1)
	{
		cout << endl << "Declan Doyle Data Structures and Algorithms 2 Coursework Program" << endl;
		cout << "Enter -1 to exit any menu screen" << endl;
		cout << "This program calculates and generates the mandlebrot set on either the CPU, or GPU" << endl;
		cout << endl << "Please choose how the mandlebrot will be calculated:" << endl;
		cout << "1. CPU" << endl;
		cout << "2. GPU" << endl;
		cin >> menu;
		cout << endl << endl;
		menuChoice(menu);

	}


}

void menuChoice(int menu)
{
  if (menu == 1)
  {
    CPU();
  }
  else if (menu == 2)
  {
    GPU();
  }
  else if (menu > 2)
  {
      cout << "Invalid Choice" << endl;
  }

}

void CPU()
{

	ofstream my_file("MandelbrotCPU.csv");

	int threads = 0;

	cout << endl << endl << "Please enter the number of threads to use to generate the mandlebrot. Please choose 1-16:" << endl;
	cin >> threads;

	while (threads > 16 || threads < 1)
	{
		cout << endl << "Please enter a thread count between 1 and 16" << endl;
		cin >> threads;
	}

	cout << endl << "Please wait..." << endl;

	int z = HEIGHT / threads;

	thread *compute_thread = new thread[threads];

	// This shows the whole set.
	int j = 0;
	int x = 0;
	for (int i = z; i <= HEIGHT, j <= (HEIGHT - z); i += z, j += z) {

		the_clock::time_point start = the_clock::now();

		compute_thread[x] = thread(compute_mandelbrot_CPU, -2.0, 1.0, 1.125, -1.125, j, i);

		the_clock::time_point end = the_clock::now();

		auto time_taken = duration_cast<microseconds>(end - start).count();
		
		my_file << time_taken << ", " << x << endl;
		
		//cout << "Computing the Mandelbrot set took " << time_taken << " microseconds." << endl << endl;
		x++;
	}

	for (int i = 0; i < threads; i++)
	{
		compute_thread[i].join();
	}

	write_tga("output.tga");

}

void compute_mandelbrot_CPU(double left, double right, double top, double bottom, double YStart, double YEnd)
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

			if (iterations == MAX_ITERATIONS)
			{
				// z didn't escape from the circle.
				// This point is in the Mandelbrot set.
				image[y][x] = 0x000000; // black
			}
			else
			{
				// z escaped within less than MAX_ITERATIONS
				// iterations. This point isn't in the set.
				image[y][x] = 0xFFFFFF; // white
			}
		}
	}
}

void GPU()
{
	ofstream my_file("MandelbrotGPU.csv");

  cout << "Please wait..." << endl;

	// Start timing
	the_clock::time_point start = the_clock::now();

	// This shows the whole set.
	compute_mandelbrot_GPU(-2.0, 1.0, 1.125, -1.125);

	// This zooms in on an interesting bit of detail.
	//compute_mandelbrot(-0.751085, -0.734975, 0.118378, 0.134488);

	// Stop timing
	the_clock::time_point end = the_clock::now();

	// Compute the difference between the two times in milliseconds
	auto time_taken = duration_cast<microseconds>(end - start).count();
	
	my_file << time_taken << endl;

	//cout << "Computing the Mandelbrot set took " << time_taken << " microseconds." << endl << endl;

	write_tga("output.tga");
}

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

Complex1 c_add(Complex1 c1, Complex1 c2) restrict(cpu, amp) // restrict keyword - able to execute this function on the GPU and CPU
{
	Complex1 tmp;
	float a = c1.x;
	float b = c1.y;
	float c = c2.x;
	float d = c2.y;
	tmp.x = a + c;
	tmp.y = b + d;
	return tmp;
} // c_add

float c_abs(Complex1 c) restrict(cpu, amp)
{
	return concurrency::fast_math::sqrt(c.x*c.x + c.y*c.y);
} //c_abs

Complex1 c_mul(Complex1 c1, Complex1 c2) restrict(cpu, amp)
{
		Complex1 tmp;
		float a = c1.x;
		float b = c1.y;
		float c = c2.x;
		float d = c2.y;
		tmp.x = a*c - b*d;
		tmp.y = b*c + a*d;
		return tmp;
} // c_mul

void compute_mandelbrot_GPU(double left, double right, double top, double bottom)
{

	uint32_t *pImage = &(image[0][0]);

	array_view<uint32_t, 2> a(HEIGHT, WIDTH, pImage);

	a.discard_data();

	parallel_for_each(a.extent, [=](concurrency::index<2> idx) restrict(amp) {

		int x = idx[1];
		int y = idx[0];

		// Work out the point in the complex plane that
		// corresponds to this pixel in the output image.
		Complex1 c;
		c.x = (left + (x * (right - left) / WIDTH));
		c.y = (top + (y * (bottom - top) / HEIGHT));


		// Start off z at (0, 0).
		Complex1 z;
		z.x = 0.0;
		z.y = 0.0;

		// Iterate z = z^2 + c until z moves more than 2 units
		// away from (0, 0), or we've iterated too many times.
		int iterations = 0;
		while (c_abs(z) < 2.0 && iterations < MAX_ITERATIONS)
		{
			z = c_add( (c_mul(z, z)), c);

			++iterations;
		}

		if (iterations == MAX_ITERATIONS)
		{
			// z didn't escape from the circle.
			// This point is in the Mandelbrot set.
			a[y][x] = 0x000000; // black
		}
		else
		{
			// z escaped within less than MAX_ITERATIONS
			// iterations. This point isn't in the set.
			a[y][x] = 0xFFFFFF; // white
		}


	});


}
