#pragma once
#include <functional>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <memory>

#include "mvector.h"
#include "spaces/ndspace.h"
#include "spaces/ndsymspace.h"



template <typename T, int rank_>
class Array
{
	T	*data;
	mVector<unsigned, rank_> dim;
	mVector<size_t, rank_> M;
	size_t d_size;

	std::shared_ptr<T> dyn;

//	std::function<size_t (mVector<unsigned, rank_> const &)> metric;

	public:
		enum { rank = rank_ };
		typedef T Type;
		typedef mVector<unsigned, rank> Index;
//		typedef std::function<size_t (Index const &)> Metric;

		Array() {}

		Array(mVector<unsigned, rank> const &dim_):
			dim(dim_)
		{
			d_size = std::accumulate(dim.begin(), dim.end(),
				size_t(1), std::multiplies<size_t>());

			size_t rtm = 1;
			for (unsigned k = 0; k < rank; ++k)
			{ M[k] = rtm; rtm *= dim[k]; }

			dyn = std::shared_ptr<T>(new T[d_size]);
			data = dyn.get();
		}


		Array(T *data_, mVector<unsigned, rank_> const &dim_):
			data(data_), dim(dim_)
		{
			size_t rtm = 1;
			for (unsigned k = 0; k < rank; ++k)
			{ M[k] = rtm; rtm *= dim[k]; }	

			d_size = std::accumulate(dim.begin(), dim.end(),
				size_t(1), std::multiplies<size_t>());

//			metric = [&M] (Index const &X) -> size_t
//			{
//				return dot(M, X);
//			};
		}

		void set_data(T *d) { data = d; }
//		void set_metric(Metric &m, Index const &dim_) 
//		{
//			dim = dim_;
//			metric = m; 
//
//			size = metric(dim);	
//		}

		void set_dim(Index const &dim_)
		{
			dim = dim_;
			size_t rtm = 1;
			for (unsigned k = 0; k < rank; ++k)
			{ M[k] = rtm; rtm *= dim[k]; }	
			//std::cout << dim << " | " << M << std::endl;

			d_size = std::accumulate(dim.begin(), dim.end(),
				size_t(1), std::multiplies<size_t>());
		}

		size_t size() const { return d_size; }

		size_t metric(Index const &X)
		{
			return dot(X, M);
		}

		T &operator[](Index const &X)
		{
			size_t offset = metric(X);
		//	std::cout << offset << std::endl;
			return data[offset];
		}

		template <typename F>
		void map_indices(F f)
		{
			Spaces::NdSpace<rank> G(dim);
			std::vector<mVector<unsigned, rank>> D(d_size);
			std::copy(G.begin(), G.end(), D.begin());

			T *loc = data;

			#pragma omp parallel for
			for (size_t i = 0; i < d_size; ++i)
			{
				loc[i] = f(D[i]);
			}
		}

		template <typename F>
		void map_sym_indices(F f)
		{
			Spaces::NdSymSpace<rank> G(dim);
			std::vector<mVector<unsigned, rank>> D(d_size);
			std::copy(G.begin(), G.end(), D.begin());

			unsigned s = (dim[0] * (dim[0] + 1)) / 2;
			#pragma omp parallel for
			for (size_t i = 0; i < s; ++i)
			{
				T v = f(D[i]); (*this)[D[i]] = v;
				while (std::prev_permutation(D[i].begin(), D[i].end()))	
					(*this)[D[i]] = v;
			}
		}

		T *begin() { return data; }
		T const *begin() const { return data; }
		T *end() { return data + d_size; }
		T const *end() const { return data + d_size; }
};

