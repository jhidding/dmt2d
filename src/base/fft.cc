#include "fft.hh"
#include "misc.hh"
#include <vector>

using namespace Fourier;

Transform::Transform(std::vector<int> const &shape):
	size(System::product(shape)), in(size), out(size)
{
	std::vector<int> ishape(shape.begin(), shape.end());

	d_plan_fwd = fftw_plan_dft(ishape.size(), ishape.data(),
		reinterpret_cast<fftw_complex *>(in.data()), 
		reinterpret_cast<fftw_complex *>(out.data()), FFTW_FORWARD, FFTW_ESTIMATE);

	d_plan_bwd = fftw_plan_dft(ishape.size(), ishape.data(),
		reinterpret_cast<fftw_complex *>(in.data()),
		reinterpret_cast<fftw_complex *>(out.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
}

Transform::~Transform()
{
	fftw_destroy_plan(d_plan_fwd);
	fftw_destroy_plan(d_plan_bwd);
}

void Transform::forward() { fftw_execute(d_plan_fwd); }
void Transform::backward() { fftw_execute(d_plan_bwd); }

