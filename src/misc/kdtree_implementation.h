#include "kdtree.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iterator>
#include <iostream>

namespace kdTree {

template <typename Point>
inline std::function<bool (Point const &, Point const &)> compare_element(
	unsigned i,
	std::function<double (Point const &, int)> &index)
{
	return [i, &index] (Point const &A, Point const &B) -> bool 
		{ return index(A, i) < index(B, i); };
}

/*
 * BoundingBox ---------------------------------------------------
 */

template <typename Point, int R>
BoundingBox<Point, R>::BoundingBox(
	System::mVector<double, R> const &X1_, 
	System::mVector<double, R> const &X2_,
	std::function<double (Point const &, int)> index_)
	: X1(X1_), X2(X2_), index(index_) {}

template <typename Point, int R>
template <typename Iter>
BoundingBox<Point, R>::BoundingBox(Iter begin, Iter end,
	std::function<double (Point const &, int)> index_):
	index(index_)
{
	for (unsigned i = 0; i < R; ++i)
	{
		X1[i] = index(*std::min_element(begin, end, compare_element<Point>(i, index)), i);
		X2[i] = index(*std::max_element(begin, end, compare_element<Point>(i, index)), i);
	}
}

template <typename Point, int R>
double BoundingBox<Point, R>::half(unsigned i)
{
	return (X2[i] + X1[i]) / 2;
}

/*
 * Tree implementation -------------------------------------------------
 */

template <typename Point, int R>
template <typename Iter>
Tree<Point, R>::Tree(Iter begin, Iter end, 
	std::function<double (Point const &, int)> idx, 
	int dim, int depth)
	: Branch<Point, R>(begin, end, idx)
{
	double boundary = BoundingBox<Point, R>::half(dim);
	Iter mid = std::partition(begin, end,
		[&] (Point const &p)
	{
		return idx(p, dim) < boundary;
	});
	
	if ((end - mid) < 32 or depth > 128)
		B.second = Branch_ptr(dynamic_cast<Branch<Point, R> *>(
			new Leaf<Point, R, Iter>(mid, end, idx)));
	else
		B.second = Branch_ptr(dynamic_cast<Branch<Point, R> *>(
			new Tree<Point, R>(mid, end, idx, (dim + 1) % R, depth+1)));

	if ((mid - begin) < 32 or depth > 128)
		B.first = Branch_ptr(dynamic_cast<Branch<Point, R> *>(
			new Leaf<Point, R, Iter>(begin, mid, idx)));
	else
		B.first = Branch_ptr(dynamic_cast<Branch<Point, R> *>(
			new Tree<Point, R>(begin, mid, idx, (dim + 1) % R, depth+1)));
}


template <typename Point, int R>
size_t Tree<Point, R>::count_if(Predicate<Point, R> const &pred) const
{
	if (pred(*this))
		return  B.first->count_if(pred) +
			B.second->count_if(pred);
	else
		return 0;
}

template <typename Point, int R>
void Tree<Point, R>::traverse(Visitor<Point> &visit, Predicate<Point, R> const &pred) const
{
	if (not pred(*this)) return;

	B.first->traverse(visit, pred);
	B.second->traverse(visit, pred);
}

template <typename Point, int R>
Point const &Tree<Point, R>::nearest_neighbour(Distance<Point, R> const &dist) const
{
	if (dist(*B.second) > dist(*B.first))
	{
		Point const &U = B.first->nearest_neighbour(dist);
		if (dist(*B.second) > dist(U))
			return U;
		
		Point const &V = B.second->nearest_neighbour(dist);
		return (dist(V) > dist(U) ? U : V);
	}
	else
	{
		Point const &U = B.second->nearest_neighbour(dist);
		if (dist(*B.first) > dist(U))
			return U;
		
		Point const &V = B.first->nearest_neighbour(dist);
		return (dist(V) > dist(U) ? U : V);
	}
}

/*
 * Leaf implementation ---------------------------------------------
 */

template <typename Point, typename F>
std::function<bool (Point const &)> Pack(F const &pred)
{
	return [&pred] (Point const &A) -> bool { return pred(A); };
}

template <typename Point, int R, typename Iter>
Leaf<Point, R, Iter>::Leaf(Iter begin, Iter end,
	std::function<double (Point const &, int)> idx)
	: Branch<Point, R>(begin, end, idx), I(begin, end)
{}

template <typename Point, int R, typename Iter>
size_t Leaf<Point, R, Iter>::count_if(Predicate<Point, R> const &pred) const
{
	if (pred(*this))
	{
		return std::count_if(I.first, I.second, Pack<Point>(pred));
	}
	else
	{
		return 0;
	}
}

template <typename Point, int R, typename Iter>
void Leaf<Point, R, Iter>::traverse(Visitor<Point> &visit, Predicate<Point, R> const &pred) const
{
	if (not pred(*this)) return;
	Iter pos = I.first;
	while ((pos = std::find_if(pos, I.second, Pack<Point>(pred))) != I.second)
	{
		visit(*pos);
		++pos;
	}
}

template <typename Point, int R, typename Iter>
Point const &Leaf<Point, R, Iter>::nearest_neighbour(Distance<Point, R> const &dist) const
{
	Point const *P = &(*(I.first));
	double d = dist(*P);

	std::for_each(I.first, I.second,
		[&] (Point const &X)
	{
		double a= dist(X);
		if (a < d)
		{
			d = a;
			P = &X;
		}
	});

	return *P;
}

} // namespace KdTree

